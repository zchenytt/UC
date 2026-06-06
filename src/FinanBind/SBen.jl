"""
pureBenCut: theta > pi'x + Cn
Cn := obj_ast - pi_ast'x_che
pi'x - theta < -Cn
"""
module SBen
import Gurobi, JuMP
import ..Settings, ..Case2383, ..WindGen, ..Static, ..Models
const SB_CUT_COT = 0.05

function _1(ch, s, Θs, n, p)
    (; Cd, x_che, N) = n # (x_che, Θs) here is local-stable
    pi_hat = _sb_core!(n.Xl, n, x_che)
    ObjVal = Settings.getmodeldblattr(n, "ObjVal")
    sec =Settings.getmodeldblattr(n, "Runtime")
    vio = ObjVal - Θs
    @ccall(printf("s=%d, vio=%.3e, sec=%.1e\n"::Cstring; s::Cint, vio::Cdouble, sec::Cdouble)::Cint)
    Cd[1] = Inf # is non-violating
    if vio > SB_CUT_COT
        Cd[end] = Cn = ObjVal - x_che'pi_hat
        Cd[end-1] = -1        # to θ
        Cd[1:end-2] .= pi_hat # to x
    end
    put!(ch, s)
end
_spawn1s(s,Θs,n,Xl,CanSpawn,ch,p) = (n.x_che .= Xl; CanSpawn[s]=false; Threads.@spawn(_1(ch, s, Θs, n, p)))
function trainsb_shortly(mst, sub, pas, NP, lbv, mstreopttime, subtermigap, timensinalgo3, timensinalgo3mst; Times = 20.1)
    (; o, S, N, Θ, Xl, Bv), ch = mst, Channel{Int}();
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 1+length(Θ), N, Bv)
    _111(mst)
    iFloat64, CanSpawn, iUnBn = (Times)S, trues(S), 0
    for s=1:NP _spawn1s(s, Θ[s], sub[s], Xl, CanSpawn, ch, nothing) end; iUnBn += NP
    tabs0 = 1e-9time_ns()
    while true
        if isready(ch) === false
            Settings.opt_ass_opt(mst);          load_che(mst)
            lb = Settings.getmodeldblattr(mst,  "ObjBound")
            mstruntime = Settings.getmodeldblattr(mst, "Runtime")
            push!(timensinalgo3mst, 1e-9time_ns() - tabs0)
            push!(lbv, lb)
            push!(mstreopttime, mstruntime)
            @ccall(printf("i=%6d, lb=%.5e, sec=%.1e\n"::Cstring; iUnBn::Cint, lb::Cdouble, mstruntime::Cdouble)::Cint)
            iUnBn > iFloat64 && break
        end; s = take!(ch)
        # push!(timensinalgo3, 1e-9time_ns() - tabs0)
        # push!(subtermigap, Settings.getmodeldblattr(sub[s], "MIPGap"))
        n, CanSpawn[s] = sub[s], true; Cd = n.Cd
        Cd[1]==Inf || Gurobi.GRBaddconstr(o,N+1,n.Ci,Cd,Cchar('<'),-Cd[end],C_NULL)
        while true # took one, so you must spawn a new one
            iUnBn += 1 # make sure that every scene has re-optimize opportunity
            s = Settings.to_1S(iUnBn, S)
            CanSpawn[s] && (_spawn1s(s, Θ[s], sub[s], Xl, CanSpawn, ch, nothing);break)
        end
    end
    for _=1:NP
        isready(ch) || Gurobi.GRBoptimize(o)
        n=sub[take!(ch)]; Cd = n.Cd
        Cd[1]==Inf || Gurobi.GRBaddconstr(o,N+1,n.Ci,Cd,Cchar('<'),-Cd[end],C_NULL)
    end
    _111(mst)
end

function _m1(ch, s, sub, a...)
    Models.sub!(sub, s, a...)
    put!(ch, s)
end
function prepare_models(S, t, T, NP)
    WD = WindGen.get_case2383(S, t);
    CaD = Case2383.get_Case_Dict();
    F = Case2383.get_PTDF(CaD);
    LD = Case2383.get_load(CaD);
    GD = Case2383.get_gen(CaD);
    Line = Case2383.branch_by_zone(CaD);
    genv, e_v = Settings.Env(), Settings.Env(S)
    xH,uH,pAH = Static.get_History_2383(t, F, Line, WD, LD, GD, genv);
    LineˈP = Line.P;
    mst = Models.mst!(genv, t, T, WD, GD, uH, xH, pAH)
    ch = Channel{Int}()
    sub = Vector{Models.subTy}(undef, S)
    for s=1:NP
        e_v_s = e_v[s]
        Threads.@spawn(_m1(ch, s, sub, T, e_v_s, t, LineˈP, F, WD, LD, GD, uH, xH, pAH))
    end
    i = NP
    while true
        i === S && break
        scene = take!(ch)
        rand() < 0.005 && @ccall(printf("s=%5d\n"::Cstring; scene::Cint)::Cint)
        s = i += 1
        e_v_s = e_v[s]
        Threads.@spawn(_m1(ch, s, sub, T, e_v_s, t, LineˈP, F, WD, LD, GD, uH, xH, pAH))
    end
    for _=1:NP
        scene = take!(ch)
        rand() < 0.005 && @ccall(printf("s=%5d\n"::Cstring; scene::Cint)::Cint)
    end
    pas = nothing
    mst, sub, pas
end

turn_on_mst_CG(mst) = config_mst_CG(mst, 1., -1., 'C')
function config_mst_CG(mst, csumRHS, xMatCoeff, xVType)
    (; o, N, Θ, Bv, Cd, CCi, VCi) = mst
    Settings.setxdblattrelement(mst, N, "RHS", csumRHS)
    Cd .= xMatCoeff; Bv .= xVType
    Gurobi.GRBchgcoeffs(o, N, CCi, VCi, Cd)
    Gurobi.GRBsetcharattrarray(o, "VType", 1+length(Θ), N, Bv)
end


function _sb_core!(pi_hat, n, x_che)
    solve_subLP_no_bias(n, x_che)
    Gurobi.GRBgetdblattrarray(n.o, "RC", 1, n.N, pi_hat) # mutated!
    pi_hat
end
function solve_subLP_no_bias(n, x_che)
    _no_bias_obj_and_fix_at_che(n, x_che)
    t = Settings.opt_and_ter(n)
    t == 2 || (@error("solve_subLP_no_bias> termination = $t"); error())
end
function _no_bias_obj_and_fix_at_che(n, x_che)
    (; Xl2, N, o) = n
    Gurobi.GRBsetdblattrarray(o, "LB", 1, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 1, N, x_che)
end

function ini_lb!(mst, sub, pas)
    (; o, S, N, Xl, Θ), ch = mst, Channel{Int}()
    Settings.opt_ass_opt(mst)
    Gurobi.GRBgetdblattrarray(o, "X", 1+S, N, Xl)
    for (s,n)=enumerate(sub) Threads.@spawn(_ini_1(ch, s, n, Xl, nothing)) end
    for _ = sub
        isready(ch) || Gurobi.GRBoptimize(o)
        n = sub[take!(ch)]
        Gurobi.GRBaddconstr(o,N+1,n.Ci,n.Cd,Cchar('<'),-n.Cd[end],C_NULL) # true
    end
    Θ .= 1/S; Gurobi.GRBsetdblattrarray(o, "Obj", 1, S, Θ) # do only once initially
    _111(mst)
end
function _ini_1(ch, s, n, x_che, p)
    (; N, Cd, Ci) = n
    pi_hat = _sb_core!(n.Xl, n, x_che)
    ObjVal = Settings.getmodeldblattr(n, "ObjVal")
    Cd[end] = Cn = ObjVal - x_che'pi_hat
    Cd[N+1] = -1        # to θ
    Cd[1:N] .= pi_hat   # to x
    put!(ch, s)
end
function _111(mst)
    Settings.opt_ass_opt(mst)
    lb, _ = Settings.getmodeldblattr(mst, "ObjBound"), load_che(mst)
    @ccall(printf("lb = %.5e\n"::Cstring; lb::Cdouble)::Cint)
end
function load_che(mst)
    (; o, S, N, Θ, Xl) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 1, S, Θ)
    Gurobi.GRBgetdblattrarray(o, "X", 1+S, N, Xl)
end

"evaluation is very different from adding cuts---It must generate Binary Trial"
function _eval1(s, n, x_che, ch)
    (; o, N) = n
    Gurobi.GRBsetdblattrarray(o, "LB", 1, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 1, N, x_che)
    t = Settings.opt_and_ter(n)
    if t != 2
        @error("s = $s, subproblem terminate with $t")
        error()
    end
    put!(ch, s)
end
function _eval(mst, sub; ndisplay = 2)
    (; o, S, N, Θ, Xl, Bv, Xn) = mst
    genv, ch = Gurobi.GRBgetenv(o), Channel{Int}()
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", 1+length(Θ), N, Bv)
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(mst.o), "OutputFlag", 1)
    JuMP.optimize!(mst.m)
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(mst.o), "OutputFlag", 0)
    lb, Nsol = Settings.getmodeldblattr(mst, "ObjBound"), Settings.getmodelintattr(mst,"SolCount")
    Rt = Settings.getmodeldblattr(mst, "Runtime")
    printstyled("eval_MIP_Runtime = $Rt\n"; color = :magenta)
    for (i, _) = zip(1:Nsol, 1:ndisplay)
        Gurobi.GRBsetintparam(genv, "SolutionNumber", i-1)
        com, v = Settings.getxdblattrelement(mst, 0, "PoolNX"), Xn[i]
        a, _ = Threads.Atomic{Float64}(com), Gurobi.GRBgetdblattrarray(o, "PoolNX", 1+S, N, v)
        @. v = round(v);        for (s, n)=enumerate(sub) Threads.@spawn(_eval1(s, n, v, ch)) end
        for _=sub Threads.atomic_add!(a, Settings.getmodeldblattr(sub[take!(ch)], "ObjVal")/S) end
        ub = a.value; agap = ub - lb; rgap = agap / ub
        printstyled("Sol=$i: lb = $lb, agap = $agap, rgap = $rgap\n"; color = :magenta)
    end
end


end
