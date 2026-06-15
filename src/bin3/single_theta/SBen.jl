module SBen
import Gurobi, JuMP
import ..Settings, ..Case2383, ..WindGen, ..Static, ..Models
const NphyCore = Sys.CPU_THREADS ÷ 2

_1(ch, s, Θs, n) = (_sb_core!(n.Xl, n, n.x_che); put!(ch, s))
_spawn1s(s,Θs,n,Xl,ch) = (n.x_che .= Xl; Threads.@spawn(_1(ch, s, Θs, n)))
function trainsb_shortly(mst, sub, NP, sampleData)
    (; o, S, N, Xl, Bv, Pi, Ci, Cd), ch = mst, Channel{Int}(); Threads.@threads(for n = sub
        _e7 = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(_e7, "MIPGap", 1e-3)    
        Gurobi.GRBsetdblparam(_e7, "TimeLimit", 15.)    
    end); Bv .= 'C'; Gurobi.GRBsetcharattrarray(o, "VType", 1, N, Bv); _111(mst)
    for ite = 1:240 # 3600s / 15s = 240
        Pi .= ΣobN = 0; Θ_che = Settings.getmodeldblattr(mst, "ObjVal")
        for s=1:NP _spawn1s(s, nothing, sub[s], Xl, ch) end; ı = NP
        while ı < S
            s=take!(ch); n=sub[s]; Πs, ObNs = n.Xl, Settings.getmodeldblattr(n, "ObjBound")
            Pi .+= Πs; ΣobN += ObNs; s = ı += 1; _spawn1s(s, nothing, sub[s], Xl, ch)
        end
        for _=1:NP
            s=take!(ch); n=sub[s]; Πs, ObNs = n.Xl, Settings.getmodeldblattr(n, "ObjBound")
            Pi .+= Πs; ΣobN += ObNs
        end
        vio = (Pi'Xl + ΣobN - (S)Θ_che) / S
        vio > 0 || (printstyled("vio = $vio < 0, thus break"; color = :red, bold = true); break)
        Cd[1:N] .= Pi; Cd[N+1] = -S; Gurobi.GRBaddconstr(o,N+1,Ci,Cd,Cchar('<'),-ΣobN,C_NULL)
        Settings.opt_ass_opt(mst, "master")
        load_che(mst)
        lb = Settings.getmodeldblattr(mst, "ObjBound")
        ite % 4 == 0 && println("ite=", ite, ", vio=", vio, ", lb=", lb)
    end
end

_2345(ch,s,T,sub,e_v,a...) = (Models.sub!(sub,s,T,e_v[s],a...); put!(ch,s))
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
    ch, sub = Channel{Int}(), Vector{Models.subTy}(undef, S); NP = min(NP, NphyCore)
    for s=1:NP Threads.@spawn(_2345(ch,s,T,sub,e_v,t,LineˈP,F,WD,LD,GD,xH,pAH)) end; i = NP
    @time "JuMP Modeling sub vec" while true
        i === S && (for _=1:NP take!(ch) end; break)
        s = take!(ch)
        rand() < 5e-4 && @ccall(printf("s=%5d\n"::Cstring; s::Cint)::Cint)
        s = i += 1
        Threads.@spawn(_2345(ch,s,T,sub,e_v,t,LineˈP,F,WD,LD,GD,xH,pAH))
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
    Settings.opt_ass_opt(n, "subLP_no_bias")
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

function ini_lb!(mst, sub, NP)
    (; o, S, N, Xl, Ci, Cd, Pi), ch = mst, Channel{Int}()
    Settings.opt_ass_opt(mst, "master"); Threads.@threads(for n = sub
        _e7 = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(_e7, "MIPGap", 1e-3)
        Gurobi.GRBsetdblparam(_e7, "TimeLimit", 15.)
    end); load_che(mst);  Pi .= ΣobN = 0.
    for s=1:NP; n = sub[s]; Threads.@spawn(_ini_1(ch, s, n, Xl)) end; i = NP
    while true
        i === S && break
        s = take!(ch)
        n = sub[s]
        Πs, ObNs = n.Xl, Settings.getmodeldblattr(n, "ObjBound")
        Pi .+= Πs; ΣobN += ObNs
        s = i += 1
        n = sub[s]
        Threads.@spawn(_ini_1(ch, s, n, Xl))
    end
    for _=1:NP
        s = take!(ch)
        n = sub[s]
        Πs, ObNs = n.Xl, Settings.getmodeldblattr(n, "ObjBound")
        Pi .+= Πs; ΣobN += ObNs
    end
    Cd[1:N] .= Pi; Cd[N+1] = -S; Gurobi.GRBaddconstr(o,N+1,Ci,Cd,Cchar('<'),-ΣobN,C_NULL)
    Gurobi.GRBsetdblattrelement(o, "Obj", 0, 1.) # do only once initially
    _111(mst)
end
function _ini_1(ch, s, n, x_che) # ✅
    pi_hat = _sb_core!(n.Xl, n, x_che)
    # Cd[end] = Settings.getmodeldblattr(n, "ObjBound")
    # Cd[N+1] = -1        # to θ
    # Cd[1:N] .= pi_hat   # to x
    put!(ch, s)
end

function _111(mst) # ✅
    Settings.opt_ass_opt(mst, "master")
    lb, _ = Settings.getmodeldblattr(mst, "ObjBound"), load_che(mst)
    @ccall(printf("lb = %.5e\n"::Cstring; lb::Cdouble)::Cint)
end
function load_che(mst) # ✅
    (; o, N, Xl) = mst
    Gurobi.GRBgetdblattrarray(o, "X", 1, N, Xl)
    clamp!(Xl, 0., 1.)
end

"evaluation is very different from adding cuts---It must generate Binary Trial"
_eval1(s, n, x_che, ch) = (solve_subMIP_no_bias(n, x_che); put!(ch, s))
function _eval(mst, sub, NP; ndisplay = 1)
    (; o, S, N, Xl, Xn, Bv) = mst; Threads.@threads(for n = sub
        _e7 = Gurobi.GRBgetenv(n.o)
        Gurobi.GRBsetdblparam(_e7, "MIPGap", 1e-4)    
        Gurobi.GRBsetdblparam(_e7, "TimeLimit", 55.)    
    end);       genv, ch = Gurobi.GRBgetenv(o), Channel{Int}(); MuL = 1/S; NP = min(NP, NphyCore)
    Bv .= 'B'; Gurobi.GRBsetcharattrarray(o, "VType", 1, N, Bv)
    Gurobi.GRBsetintparam(genv, "OutputFlag", 1); Gurobi.GRBsetdblparam(genv, "TimeLimit", 3e2)
    JuMP.optimize!(mst.m) # so that you may interrupt
    Gurobi.GRBsetintparam(genv, "OutputFlag", 0); Gurobi.GRBsetdblparam(genv, "TimeLimit", 1e100)
    lb, Nsol = Settings.getmodeldblattr(mst,"ObjBound"), Settings.getmodelintattr(mst,"SolCount")
    Rt = Settings.getmodeldblattr(mst, "Runtime")
    printstyled("eval_MIP_Runtime = $Rt\n"; color = :magenta)
    for (i, _) = zip(1:Nsol, 1:ndisplay)
        Gurobi.GRBsetintparam(genv, "SolutionNumber", i-1)
        v=Xn[i]; Gurobi.GRBgetdblattrarray(o,"PoolNX",1,N,v); @. v = round(v)       
        a = Threads.Atomic{Float64}(0.)
        for s=1:NP
            n = sub[s]; Threads.@spawn(_eval1(s, n, v, ch)) 
        end; k = NP
        while k < S
            _8h2(ch, sub, a, MuL); s = k += 1; n = sub[s]; Threads.@spawn(_eval1(s, n, v, ch))
        end
        for _=1:NP _8h2(ch, sub, a, MuL) end
        ub = a.value; agap = ub - lb; rgap = agap / ub
        printstyled("Sol=$i: lb = $lb < $ub = ub, agap = $agap, rgap = $rgap\n"; color = 30)
    end
end
function _8h2(ch, sub, a, MuL)
    s = take!(ch)
    n = sub[s]
    ObjVal = Settings.getmodeldblattr(n,"ObjVal")
    Gurobi.GRBreset(n.o, 0) # discard Branch-and-Bound Tree to recycle memory
    Threads.atomic_add!(a, MuL * ObjVal)
end

end