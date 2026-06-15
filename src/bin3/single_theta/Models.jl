"""
mst (single-theta) & sub

First-stage have some start-up cost, but the u variable is not the complicating variable,
So in 1st-stage we effectively have 0 cost
"""
module Models
import ..Case2383, ..Settings, JuMP, Gurobi
const SB_CUT_COT = 4.0
subTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, Bv::Vector{Int8}, Xl::Vector{Float64}, Xl2::Vector{Float64}, x_che::Vector{Float64}, N::Int64}
_0() = JuMP.AffExpr()
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end
_variable(m) = JuMP.@variable(m, lower_bound=0, upper_bound=1)
function mst!(genv, WD, GD, xH, pAH) # Now the master problem has Common = 0
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C; S=length(WD.S)
    m=Settings.Model(genv); o,refi,refd=m.moi_backend,Ref{Cint}(),Ref{Cdouble}(); ge=Gurobi.GRBgetenv(o)
    JuMP.@variable(m, #= Single-θ GuIndex0 =# θ); N=0; for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds(begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end); Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on
            if pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT))
                nothing # x1[z,g] = true
            else
                #=x1[z,g] = =#_variable(m); N += 1
            end
        else # was off
            if any(view(xHˈzg, 1:DT)) # cannot turn on
                nothing # x1[z,g] = false
            else
                #=x1[z,g] = =#_variable(m); N += 1
            end
        end
    end; Xl=fill(0.,N); Bv=fill(Cchar('C'),N); nX=2; Pi = fill(0.,N)
    # θ > 1/S * ( sum(pi_s)'x + sum(obN_s) )  _or_   sum(pi_s)'x - S'θ < -sum(obN_s)
    Ci, Cd = Cint[range(1; length=N); 0], Vector{Cdouble}(undef, N+2) # used to finally add cut
    printstyled("master built with N = $N\n"; color = 30)
    _, Xn = Gurobi.GRBsetintparam(ge, "PoolSolutions", nX), [similar(Xl) for _=1:nX]
    Gurobi.GRBsetintparam(ge, "PoolSearchMode", 2)
    (; m, o, refd, refi, S, N, Xl, Xn, Bv, Ci, Cd, Pi)
end

_x2_or_uv() = Dict{Tuple{Int,Int,Int}, Union{Bool,JuMP.VariableRef}}()
function sub!(Ve, s, T, genv, t, LineˈP, F, WD, LD, GD, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=50, ζCost=150)
    ReserveK01v, LoadKref, times,Times, lines = Case2383.Reserve_Curve[reserve_type], Case2383.Load_Curve, t+1:t+T, t:t+T, keys(LineˈP)
    WDˈS, WDˈN, WDˈPMax                              = WD.S, WD.N, WD.PMax
    LDˈi, LDˈg, LDˈn, windMatˈs                      = LD.i, LD.g, LD.n, WDˈS[s]
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈdci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.dci, GD.C
    S,m=length(WDˈS),Settings.Model(genv); o,refi,refd=m.moi_backend,Ref{Cint}(),Ref{Cdouble}(); Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o),"Presolve",2)
    egp =  Dict(i    =>_0() for i=Times);  Qˈy=_0()
    rgp = [Dict(i    =>_0() for i=Times) for z=1:4] # 1≤z≤4
    pfe =  Dict((l,i)=>_0() for i=Times  for l=lines)
    JuMP.@variable(m, qy); x1 = Dict{Tuple{Int,Int}, Union{Bool,JuMP.VariableRef}}(); N = 0
    Cnt_2nd_Bin = 0
    x2 = _x2_or_uv()                    # doesn't contain (the initial) time `t`
    u2, v2 = _x2_or_uv(), _x2_or_uv()   # contain all Times, and u2 induce start-up costs
    #= Allocate a contiguous linking vector =# for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds(begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end); Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on
            u2[z,g,t] = false; if pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT))
                x1[z,g] = true; v2[z,g,t] = false
            else
                x1[z,g] = _variable(m); N+=1
                # v2[z,g,t] = _variable(m); Cnt_2nd_Bin += 1
            end
        else # was off
            v2[z,g,t] = false; if any(view(xHˈzg, 1:DT)) # cannot turn on
                x1[z,g] = u2[z,g,t] = false
            else
                x1[z,g] = _variable(m); N+=1
                # u2[z,g,t] = _variable(m); Cnt_2nd_Bin += 1 # This is 2nd-stage inner, must allocate
            end
        end
    end; Xl, Xl2, x_che = fill(0.,N), fill(0.,N), fill(0.,N) # x_che is locally static
    #= At the initial `t`, the remnant =# for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds(begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end); Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on
            if !(pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT)))
                v2[z,g,t] = v8 = _variable(m); Cnt_2nd_Bin += 1; JuMP.@constraint(m, x1[z,g] + v8 == 1)
            end
        else # was off
            if !any(view(xHˈzg, 1:DT)) # cannot turn on
                u2[z,g,t] = v8 = _variable(m); Cnt_2nd_Bin += 1 # This is 2nd-stage inner, must allocate
                JuMP.@constraint(m, x1[z,g] == v8)
            end
        end
    end
    ### Next allocate (3)Bin variables from t+1:t+T
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds(begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end); Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        for Ofst = 1:T; let i=t+Ofst, xm1=(Ofst>1 ? x2[z,g,i-1] : x1[z,g])
            if xm1 === true # was on
                u2[z,g,i] = false; if any(==(false), view(xHˈzg, 1:UT-Ofst)) || pAHˈzg + Pmin > Sd + Rd * Ofst
                    x2[z,g,i] = true; v2[z,g,i] = false
                else
                    x2[z,g,i] = v8 = _variable(m); v2[z,g,i] = v9 = _variable(m); Cnt_2nd_Bin += 2; JuMP.@constraint(m, v8 + v9 == 1)
                end
            elseif xm1 === false # was off
                v2[z,g,i] = false; if any(view(xHˈzg, 1:DT-Ofst))
                    x2[z,g,i] = u2[z,g,i] = false
                else
                    x2[z,g,i] = u2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1 # here can avoid redundant 2nd-stage inner variable
                end
            else # is decision
                x2[z,g,i],u2[z,g,i],v2[z,g,i]=_variable(m),_variable(m),_variable(m); Cnt_2nd_Bin += 3
            end
        end end
    end; Bv = fill(Cchar('C'), N+Cnt_2nd_Bin)
    ### The rest continuous variables
    JuMP.@variables(m, begin
        0 ≤ p2[z=eachindex(GDˈn), g=eachindex(GDˈn[z])] # ✅ p := (Pmin)x + pA (But here only for pricing purpose, only for x1)
        0 ≤ pe[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times] # effective power of generator, ≤ its actual output, avoid Benders initial over generation
        0 ≤ a2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0           ≤  ϖ[z=eachindex(WDˈN),               Times] ≤ 1
        0           ≤  ζ[z=eachindex(LDˈn),               Times] ≤ 1
        -LineˈP[l]  ≤  pf[l=lines,i=Times]                       ≤ LineˈP[l]
    end);             
    #= Generators =# for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            Gwid, xHˈzg1, x1ˈzg, p2ˈzg = PMax-Pmin, xH[z][g][1], x1[z,g], p2[z,g]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC, sdC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC, GDˈdci[z][g]linC
        end; Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        #= Uvx =# JuMP.@constraints(m, begin
            u2[z,g,t]-v2[z,g,t] == x1ˈzg-xHˈzg1
            u2[z,g,t+1]-v2[z,g,t+1] == x2[z,g,t+1]-x1ˈzg
            [i=t+2:t+T], u2[z,g,i]-v2[z,g,i] == x2[z,g,i]-x2[z,g,i-1]
            u2[z,g,t]         ≤ x1ˈzg
            x1ˈzg + v2[z,g,t] ≤ 1
        end); for i=times; let j = i, e = _0(), x = x2[z,g,i]
            for _ = 1:UT
                j < t && break
                JuMP.add_to_expression!(e, u2[z,g,j])
                j -= 1
            end; JuMP.@constraint(m, e ≤ x); j = i; e = _0()
            for _ = 1:DT
                j < t && break
                JuMP.add_to_expression!(e, v2[z,g,j])
                j -= 1
            end; JuMP.@constraint(m, x + e ≤ 1)
        end end
        #= Qˈy / pe =# JuMP.@constraint(m, p2ˈzg == Pmin * x1ˈzg + a2[z,g,t])
        JuMP.@constraint(m, pe[z,g,t] ≤ p2ˈzg)
        JuMP.add_to_expression!(Qˈy, linC, p2ˈzg); for i=times
            tmp14 = x2[z,g,i]
            tmp14 isa Bool || JuMP.add_to_expression!(Qˈy, linC * Pmin, tmp14)
            JuMP.add_to_expression!(Qˈy, linC, a2[z,g,i])
            JuMP.@constraint(m, pe[z,g,i] <= Pmin * x2[z,g,i] + a2[z,g,i])
        end; for i=Times
            tmp14 = u2[z,g,i]
            tmp14 isa Bool || JuMP.add_to_expression!(Qˈy, suC, tmp14)
            tmp14 = v2[z,g,i]
            tmp14 isa Bool || JuMP.add_to_expression!(Qˈy, sdC, tmp14)
        end
        #= Reserve =# 1≤z≤4 && for i=Times JuMP.add_to_expression!(rgp[z][i], ur[z,g,i]) end
        #= Energy & LineFlow =#
            pe8 = pe[z,g,t]
            JuMP.add_to_expression!(egp[t], pe8)
            for l=lines
                (pfl, Cnl) = (pfe[l,t], F[l, node])
                JuMP.add_to_expression!(pfl, Cnl, pe8)
            end
        for i=times
            pe8 = pe[z,g,i]
            JuMP.add_to_expression!(egp[i], pe8)
            for l=lines
                (pfl, Cnl) = (pfe[l,i], F[l, node])
                JuMP.add_to_expression!(pfl, Cnl, pe8)
            end
        end
        if MmRatio > 1.1
            Ju, Jv = Su - Pmin, Sd - Pmin - Rd; JuMP.@constraints(m, begin
                a2[z,g,t] + ur[z,g,t] - pAHˈzg ≤ Ru * xHˈzg1 + Ju * u2[z,g,t]
                a2[z,g,t+1] + ur[z,g,t+1] - a2[z,g,t] ≤ Ru * x1ˈzg + Ju * u2[z,g,t+1]
                [i=t+2:t+T], a2[z,g,i] + ur[z,g,i] - a2[z,g,i-1] ≤ Ru * x2[z,g,i-1] + Ju * u2[z,g,i]
                [i=t+2:t+T], a2[z,g,i-1] - a2[z,g,i]             ≤ Rd * x2[z,g,i-1] + Jv * v2[z,g,i]
                a2[z,g,t] - a2[z,g,t+1]               ≤ Rd * x1ˈzg + Jv * v2[z,g,t+1]
                pAHˈzg - a2[z,g,t]             ≤ Rd * xHˈzg1 + Jv * v2[z,g,t]
            end)
            if UT > 1
                JuMP.@constraints(m, begin
                               a2[z,g,t] + ur[z,g,t] ≤ Gwid * x1ˈzg     - (PMax - Su) * u2[z,g,t] - (PMax - Sd) * v2[z,g,t+1]
                [i=t+1:t+T-1], a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i] - (PMax - Su) * u2[z,g,i] - (PMax - Sd) * v2[z,g,i+1]
                [i=[t+T]],     a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i] - (PMax - Su) * u2[z,g,i]
                end)
            else
                JuMP.@constraints(m, begin               
                               a2[z,g,t] + ur[z,g,t] ≤ Gwid * x1ˈzg     - max(Sd-Su, 0.) * u2[z,g,t] - (PMax - Sd) * v2[z,g,t+1]
                [i=t+1:t+T-1], a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i] - max(Sd-Su, 0.) * u2[z,g,i] - (PMax - Sd) * v2[z,g,i+1]
                [i=[t+T]],     a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i] - max(Sd-Su, 0.) * u2[z,g,i]
                               a2[z,g,t] + ur[z,g,t] ≤ Gwid * x1ˈzg     - (PMax - Su) * u2[z,g,t] - max(Su-Sd, 0.) * v2[z,g,t+1]
                [i=t+1:t+T-1], a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i] - (PMax - Su) * u2[z,g,i] - max(Su-Sd, 0.) * v2[z,g,i+1]
                [i=[t+T]],     a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i] - (PMax - Su) * u2[z,g,i]
                end)
            end
        else
            JuMP.@constraint(m,             a2[z,g,t] + ur[z,g,t] ≤ Gwid *     x1ˈzg)
            JuMP.@constraint(m, [i=times],  a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i])
        end
    end
    #= Wind/Load =# for i=Times
        egpˈi = egp[i]; for (z,Ø)=enumerate(WDˈN)
            Kˈz, ϖˈz, ϖCostˈz = windMatˈs[i,z], ϖ[z,i], length(Ø)ϖCost
            JuMP.add_to_expression!(Qˈy, ϖCostˈz, ϖˈz)
            for (g,node)=enumerate(Ø)
                pWind = Kˈz * WDˈPMax[z][g]
                JuMP.add_to_expression!(egpˈi, pWind)
                JuMP.add_to_expression!(egpˈi, -pWind, ϖˈz)
                for l=lines
                    (pfl, Cnl) = (pfe[l,i], F[l, node])
                    JuMP.add_to_expression!(pfl, Cnl * pWind)
                    JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈz)
                end
            end
        end
        ReserveK01 = ReserveK01v[i]; for (z,Ø)=enumerate(LDˈn)
            ζCostˈz, ζˈz = length(Ø)ζCost, ζ[z,i] # load shed
            JuMP.add_to_expression!(Qˈy, ζCostˈz, ζˈz)
            1≤z≤4 && (rgpˈz = rgp[z][i])
            for (g,node)=enumerate(Ø)
                Ty, ConstLoad = LDˈi[z][g], LDˈg[z][g]
                load = LoadKref[Ty][i] * ConstLoad
                if 1≤z≤4 # reserve requirements in isolated systems depends on overall demand [10.1109/TSG.2015.2469134]
                    upResDemand = *(ReserveK01, ReserveK2D, ConstLoad)
                    JuMP.add_to_expression!(rgpˈz, -upResDemand)
                    JuMP.add_to_expression!(rgpˈz, upResDemand, ζˈz)
                end
                JuMP.add_to_expression!(egpˈi, -load)
                JuMP.add_to_expression!(egpˈi, load, ζˈz)
                for l=lines
                    (pfl, Cnl) = (pfe[l,i], F[l, node])
                    JuMP.add_to_expression!(pfl, Cnl * -load)
                    JuMP.add_to_expression!(pfl, Cnl * load, ζˈz)
                end
            end
        end
    end
    JuMP.@objective(m, Min, Qˈy); JuMP.@constraints(m, begin
                                Qˈy == qy
        [z=1:4, i=Times],    rgp[z][i] ≥ 0 # Doesn't need relaxation
        [i=t:t+T],             egp[i] == 0
        [l=lines,i=Times],   pfe[l,i] == pf[l,i]
    end)
    Ve[s] = (; m, o, refd, refi, Bv, Xl, Xl2, x_che, N)
end

end