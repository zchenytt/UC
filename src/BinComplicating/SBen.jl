"""
Sometimes it's sufficient to add SBen cuts
"""
module SBen
import Gurobi, JuMP
import ..Settings, ..Case2383, ..WindGen, ..Static, ..Models
const SB_CUT_COT = 0.05

function _1(ch, s, Θs, n, p)
    (; Cd, x_che, N) = n # (x_che, Θs) here is local-stable
    pi_hat = _sb_core!(n.Xl, n, x_che)
    Obn = Settings.getmodeldblattr(n, "ObjBound")
    vio = (pi_hat'x_che + Obn) - Θs # (x_che, Θs) here is local-stable
    rand() < 0.005 && @ccall(printf("s=$s, vio=$vio\n"::Cstring; s::Cint, vio::Cdouble)::Cint)
    Cd[1] = Inf # is non-violating
    if vio > SB_CUT_COT
        Cd[end] = Obn         # height
        Cd[end-1] = -1        # to θ
        Cd[1:end-2] .= pi_hat # to x
    end
    put!(ch, s)
end
_spawn1s(s,Θs,n,Xl,CanSpawn,ch,p) = (n.x_che .= Xl; CanSpawn[s]=false; Threads.@spawn(_1(ch, s, Θs, n, p)))
function trainsb_shortly(mst, sub, pas, NP; Times = 20.1)
    (; o, S, N, Θ, Xl), ch = mst, Channel{Int}(); Threads.@threads(for n = sub
        local_e = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(local_e, "MIPGap", 9e-4)    
        Gurobi.GRBsetdblparam(local_e, "TimeLimit", 15.)    
    end); config_mst_CG(mst, 0., 0., 'C'); _111(mst)
    iFloat64, CanSpawn, iUnBn = (Times)S, trues(S), 0
    for s=1:NP _spawn1s(s, Θ[s], sub[s], Xl, CanSpawn, ch, nothing) end; iUnBn += NP
    while true
        if isready(ch) === false
            Settings.opt_ass_opt(mst);          load_che(mst)
            lb = Settings.getmodeldblattr(mst,  "ObjBound")
            rand() < 0.005 && @ccall(printf("i=%6d, lb=%.5e\n"::Cstring; iUnBn::Cint, lb::Cdouble)::Cint)
            iUnBn > iFloat64 && break
        end; s = take!(ch)
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
    mst = Models.mst!(genv, WD, GD, uH, xH, pAH)
    ch = Channel{Int}()
    sub = Vector{Models.subTy}(undef, S)
    for s=1:NP Threads.@spawn(begin
        Models.sub!(sub, s, T, e_v[s], t, LineˈP, F, WD, LD, GD, uH, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=130, ζCost=150)
        put!(ch, s)
    end) end
    i = NP
    while true
        local s
        i === S && break
        scene = take!(ch)
        rand() < 0.005 && @ccall(printf("s=%5d\n"::Cstring; scene::Cint)::Cint)
        s = i += 1; Threads.@spawn(begin
            Models.sub!(sub, s, T, e_v[s], t, LineˈP, F, WD, LD, GD, uH, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=130, ζCost=150)
            put!(ch, s)
        end)
    end
    for _=1:NP
        scene = take!(ch)
        rand() < 0.005 && @ccall(printf("s=%5d\n"::Cstring; scene::Cint)::Cint)
    end
    # pas = Vector{Models.pasTy}(undef, S); Threads.@threads for s=1:S
    #     Models.pas!(pas, s, mst.N, e_v[s])
    # end
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
    solve_subMIP(n, pi_hat)
    pi_hat
end
function solve_subMIP(n, pi_hat)
    (; Bv, Xl2, N, o) = n
    (@. Xl2 = -1.0 * pi_hat; Gurobi.GRBsetdblattrarray(o, "Obj", 1, N, Xl2))
    Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "LB", 1, N, Xl2)
    Xl2 .= 1; Gurobi.GRBsetdblattrarray(o, "UB", 1, N, Xl2)
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", 1, length(Bv), Bv)
    Settings.opt_ass_time(n)
end
function solve_subLP_no_bias(n, x_che)
    _no_bias_obj_and_fix_at_che(n, x_che)
    (; Bv, o) = n
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 1, length(Bv), Bv)
    t = Settings.opt_and_ter(n)
    t == 2 || (@error("solve_subLP_no_bias> termination = $t"); error())
end
function solve_subMIP_no_bias(n, x_che)
    _no_bias_obj_and_fix_at_che(n, x_che)
    (; Bv, o, N) = n
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", 1+N, length(Bv)-N, Bv)
    Settings.opt_ass_time(n)
end
function _no_bias_obj_and_fix_at_che(n, x_che)
    (; Xl2, N, o) = n
    Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "Obj", 1, N, Xl2)
    Gurobi.GRBsetdblattrarray(o, "LB", 1, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 1, N, x_che)
end

"This is mundane parallel"
function ini_lb!(mst, sub, pas)
    (; o, S, N, Xl, Θ), ch = mst, Channel{Int}()
    Settings.opt_ass_opt(mst); Threads.@threads(for n = sub
        local_e = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(local_e, "MIPGap", 9e-4)
        Gurobi.GRBsetdblparam(local_e, "TimeLimit", 15.)
    end);      Gurobi.GRBgetdblattrarray(o, "X", 1+S, N, Xl)
    for (s,n)=enumerate(sub) Threads.@spawn(_ini_1(ch, s, n, Xl, nothing)) end
    for _ = sub
        isready(ch) || Gurobi.GRBoptimize(o)
        n = sub[take!(ch)]
        Gurobi.GRBaddconstr(o,N+1,n.Ci,n.Cd,Cchar('<'),-n.Cd[end],C_NULL)
    end
    Θ .= 1/S; Gurobi.GRBsetdblattrarray(o, "Obj", 1, S, Θ) # do only once initially
    _111(mst)
end
function _ini_1(ch, s, n, x_che, p)
    (; N, Cd, Ci) = n
    pi_hat = n.Xl
    _sb_core!(pi_hat, n, x_che)
    Cd[end] = Settings.getmodeldblattr(n, "ObjBound")
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
_eval1(s, n, x_che, ch) = (solve_subMIP_no_bias(n, x_che); put!(ch, s))
function _eval(mst, sub; ndisplay = 2)
    (; o, S, N, Θ, Xl, Xn) = mst; Threads.@threads(for n = sub
                local_e = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(local_e, "MIPGap", 0.)    
        Gurobi.GRBsetdblparam(local_e, "TimeLimit", 60.)    
    end);       genv, ch = Gurobi.GRBgetenv(o), Channel{Int}()
    config_mst_CG(mst, 0., 0., 'B');
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(mst.o), "OutputFlag", 1)
    JuMP.optimize!(mst.m) # so that you may interrupt
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
