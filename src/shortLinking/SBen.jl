"""
Sometimes it's sufficient to add SBen cuts
"""
module SBen
import Gurobi, JuMP, Statistics
import ..Settings, ..Case2383, ..WindGen, ..Static, ..Models
const SB_CUT_COT = 0.05

function _1(ch, s, Θs, n)
    (; Cd, x_che, N) = n # (x_che, Θs) here is local-stable
    pi_hat = _sb_core!(n.Xl, n, x_che)
    Obn = Settings.getmodeldblattr(n, "ObjBound")
    vio = (pi_hat'x_che + Obn) - Θs # (x_che, Θs) here is local-stable
    rand() < 1e-4 && @ccall(printf("s=%d, vio=%.3e\n"::Cstring; s::Cint, vio::Cdouble)::Cint)
    Cd[1] = Inf # is non-violating
    if vio > SB_CUT_COT
        Cd[end] = Obn         # height
        Cd[end-1] = -1        # to θ
        Cd[1:end-2] .= pi_hat # to x
    end
    put!(ch, s)
end
_spawn1s(s,Θs,n,Xl,CanSpawn,ch) = (n.x_che .= Xl; CanSpawn[s]=false; Threads.@spawn(_1(ch, s, Θs, n)))
function trainsb_shortly(mst, sub, NP, sampleData; ROT = 1e2)
    (; o, S, N, Θ, Xl, Bv), ch = mst, Channel{Int}(); Threads.@threads(for n = sub
        _e7 = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(_e7, "MIPGap", 9e-4)    
        Gurobi.GRBsetdblparam(_e7, "TimeLimit", 15.)    
    end); Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", S, N, Bv); _111(mst)
    CanSpawn, ı, ȷ, ȷs, rotSum, mstta = trues(S), 0, 0, 0, 0., Threads.@spawn(:interactive, nothing)
    for s=1:NP _spawn1s(s, Θ[s], sub[s], Xl, CanSpawn, ch) end; ı += NP
    (msttxs,mstrotys,mstroiys,mstlbys,mstncys,subtxs,subgapys) = sampleData
    tabs0 = time_ns()
    while true
        if isready(ch) === false && (ȷ > 0 && istaskdone(mstta))
            wait(mstta); ȷs += ȷ
            load_che(mst) # update Θ, Xl at the foreground!
            t9,rot,roi,lb = 1e-9(time_ns()-tabs0), Settings.getmodeldblattr(mst,"Runtime"), Settings.getmodeldblattr(mst,"IterCount"), Settings.getmodeldblattr(mst,"ObjBound")
            push!(msttxs,t9);push!(mstrotys,rot);push!(mstroiys,roi);push!(mstlbys,lb);push!(mstncys,ȷs)
            rotSum += rot; rotSum > ROT && break
            rand() < 5e-4 && @ccall(printf("ȷ=%d, ro=%.1e, rs=%.1e, ri=%.1e\n"::Cstring; ȷ::Cint, rot::Cdouble, rotSum::Cdouble, roi::Cdouble)::Cint)
            mstta, ȷ = Threads.@spawn(:interactive, Settings.opt_ass_opt(mst)), 0
        end
        s = take!(ch); CanSpawn[s]=true; n=sub[s]; Cd = n.Cd
        t9, subgap = 1e-9(time_ns()-tabs0), Settings.getmodeldblattr(n,"MIPGap")
        push!(subtxs,t9);push!(subgapys,subgap)
        Cd[1] === Inf || (Gurobi.GRBaddconstr(o,N+1,n.Ci,Cd,Cchar('<'),-Cd[end],C_NULL); ȷ+=1)
        while true # took one, so you must spawn a new one
            ı += 1; s = Settings.to_1S(ı, S)
            CanSpawn[s] && (_spawn1s(s, Θ[s], sub[s], Xl, CanSpawn, ch);break)
        end
    end
    for _=1:NP
        isready(ch) || Gurobi.GRBoptimize(o)
        n=sub[take!(ch)]; Cd = n.Cd
        Cd[1] === Inf || Gurobi.GRBaddconstr(o,N+1,n.Ci,Cd,Cchar('<'),-Cd[end],C_NULL)
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
    xH,pAH = Static.get_History_2383(t, F, Line, WD, LD, GD, genv);
    LineˈP = Line.P;
    mst = Models.mst!(genv, WD, GD, xH, pAH)
    ch, sub = Channel{Int}(), Vector{Models.subTy}(undef, S)
    for s=1:NP Threads.@spawn(begin
        Models.sub!(sub, s, T, e_v[s], t, LineˈP, F, WD, LD, GD, xH, pAH)
        put!(ch, s)
    end) end
    i = NP
    while true
        local s
        i === S && break
        scene = take!(ch)
        rand() < 5e-4 && @ccall(printf("s=%5d\n"::Cstring; scene::Cint)::Cint)
        s = i += 1; Threads.@spawn(begin
            Models.sub!(sub, s, T, e_v[s], t, LineˈP, F, WD, LD, GD, xH, pAH)
            put!(ch, s)
        end)
    end
    for _=1:NP
        scene = take!(ch)
        rand() < 5e-4 && @ccall(printf("s=%5d\n"::Cstring; scene::Cint)::Cint)
    end
    mst, sub
end

function _sb_core!(pi_hat, n, x_che) # ✅
    solve_subLP_no_bias(n, x_che)
    Gurobi.GRBgetdblattrarray(n.o, "RC", 1, n.N, pi_hat) # mutated!
    solve_subMIP(n, pi_hat)
    pi_hat
end
function solve_subMIP(n, pi_hat) # ✅
    (; Bv, Xl2, N, o) = n
    (@. Xl2 = -1.0 * pi_hat; Gurobi.GRBsetdblattrarray(o, "Obj", 1, N, Xl2))
    Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "LB", 1, N, Xl2)
    Xl2 .= 1; Gurobi.GRBsetdblattrarray(o, "UB", 1, N, Xl2)
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", 1, length(Bv), Bv)
    Settings.opt_ass_time(n)
end
function solve_subLP_no_bias(n, x_che) # ✅
    _no_bias_obj_and_fix_at_che(n, x_che)
    (; Bv, o) = n
    Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 1, length(Bv), Bv)
    Settings.opt_ass_opt(n)
end
function solve_subMIP_no_bias(n, x_che) # ✅
    _no_bias_obj_and_fix_at_che(n, x_che)
    (; Bv, o, N) = n
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", 1+N, length(Bv)-N, Bv)
    Settings.opt_ass_time(n)
end
function _no_bias_obj_and_fix_at_che(n, x_che) # ✅
    (; Xl2, N, o) = n
    Xl2 .= 0; Gurobi.GRBsetdblattrarray(o, "Obj", 1, N, Xl2)
    Gurobi.GRBsetdblattrarray(o, "LB", 1, N, x_che)
    Gurobi.GRBsetdblattrarray(o, "UB", 1, N, x_che)
end

"This is mundane parallel"
function ini_lb!(mst, sub)
    (; o, S, N, Xl, Θ), ch = mst, Channel{Int}()
    Settings.opt_ass_opt(mst); Threads.@threads(for n = sub
        _e7 = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(_e7, "MIPGap", 9e-4)
        Gurobi.GRBsetdblparam(_e7, "TimeLimit", 15.)
    end);      Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
    for (s,n)=enumerate(sub) Threads.@spawn(_ini_1(ch, s, n, Xl, nothing)) end
    for _ = sub
        isready(ch) || Gurobi.GRBoptimize(o)
        n = sub[take!(ch)]
        Gurobi.GRBaddconstr(o,N+1,n.Ci,n.Cd,Cchar('<'),-n.Cd[end],C_NULL)
    end
    Θ .= 1/S; Gurobi.GRBsetdblattrarray(o, "Obj", 0, S, Θ) # do only once initially
    _111(mst)
end
function _ini_1(ch, s, n, x_che, p) # ✅
    (; N, Cd, Ci) = n
    pi_hat = _sb_core!(n.Xl, n, x_che)
    Cd[end] = Settings.getmodeldblattr(n, "ObjBound")
    Cd[N+1] = -1        # to θ
    Cd[1:N] .= pi_hat   # to x
    put!(ch, s)
end
function _111(mst) # ✅
    Settings.opt_ass_opt(mst)
    lb, _ = Settings.getmodeldblattr(mst, "ObjBound"), load_che(mst)
    @ccall(printf("lb = %.5e\n"::Cstring; lb::Cdouble)::Cint)
end
function load_che(mst) # ✅
    (; o, S, N, Θ, Xl) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 0, S,  Θ)
    Gurobi.GRBgetdblattrarray(o, "X", S, N, Xl)
end

"evaluation is very different from adding cuts---It must generate Binary Trial"
_eval1(s, n, x_che, ch) = (solve_subMIP_no_bias(n, x_che); put!(ch, s))
function _eval(mst, sub; ndisplay = 2)
    (; o, S, N, Θ, Xl, Xn, Bv) = mst; Threads.@threads(for n = sub
        _e7 = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(_e7, "MIPGap", 0.)    
        Gurobi.GRBsetdblparam(_e7, "TimeLimit", 60.)    
    end);       genv, ch = Gurobi.GRBgetenv(o), Channel{Int}()
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", S, N, Bv)
    Gurobi.GRBsetintparam(genv, "OutputFlag", 1); Gurobi.GRBsetdblparam(genv, "TimeLimit", 3e2)
    JuMP.optimize!(mst.m) # so that you may interrupt
    Gurobi.GRBsetintparam(genv, "OutputFlag", 0); Gurobi.GRBsetdblparam(genv, "TimeLimit", 1e100)
    lb, Nsol = Settings.getmodeldblattr(mst,"ObjBound"), Settings.getmodelintattr(mst,"SolCount")
    Rt = Settings.getmodeldblattr(mst, "Runtime")
    printstyled("eval_MIP_Runtime = $Rt\n"; color = :magenta)
    for (i, _) = zip(1:Nsol, 1:ndisplay)
        Gurobi.GRBsetintparam(genv, "SolutionNumber", i-1)
        v = Xn[i]; Gurobi.GRBgetdblattrarray(o, "PoolNX", 0, S, Θ)
        Gurobi.GRBgetdblattrarray(o, "PoolNX", S, N, v); @. v = round(v)
        com = Settings.getmodeldblattr(mst, "ObjVal") - Statistics.mean(Θ)
        a = Threads.Atomic{Float64}(com)
        for (s, n)=enumerate(sub) Threads.@spawn(_eval1(s, n, v, ch)) end
        for _=sub Threads.atomic_add!(a, Settings.getmodeldblattr(sub[take!(ch)], "ObjVal")/S) end
        ub = a.value; agap = ub - lb; rgap = agap / ub
        printstyled("Sol=$i: lb = $lb < $ub = ub, agap = $agap, rgap = $rgap\n"; color = :magenta)
    end
end

end