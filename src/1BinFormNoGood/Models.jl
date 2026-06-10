"""
For current-period binary complicating
"""
module Models
import ..Case2383, ..Settings, JuMP, Gurobi
subTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, Ci::Vector{Int32}, Cd::Vector{Float64}, Bv::Vector{Int8}, Xl::Vector{Float64}, Xl2::Vector{Float64}, x_che::Vector{Float64}, N::Int64}
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end
_0() = JuMP.AffExpr(0.)
_variable(m) = JuMP.@variable(m, lower_bound=0, upper_bound=1)
function mst!(genv, WD, GD, xH, pAH)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Common, S, m = _0(), length(WD.S), Settings.Model(genv)
    o,refi,refd=m.moi_backend,Ref{Cint}(),Ref{Cdouble}(); ge=Gurobi.GRBgetenv(o)
    JuMP.@variable(m, θ[1:S]); x1, Xl = Dict{Tuple{Int, Int}, Union{Bool, JuMP.VariableRef}}(), Float64[]
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on, then don't have start-up cost
            if Pmin + pAHˈzg > Sd || any(==(false), view(xHˈzg, 1:UT))
                x1[z,g] = true
            else
                x1[z,g] = _variable(m); push!(Xl, 0)
            end
        else # was off
            if any(view(xHˈzg, 1:DT)) # cannot turn on, then don't have start-up cost
                x1[z,g] = false
            else
                x1[z,g] = also_u = _variable(m); push!(Xl, 0)
                JuMP.add_to_expression!(Common, suC, also_u)
            end
        end
    end; N=length(Xl); Bv=fill(Cchar('C'),N); Θ=fill(-Inf,S); nX=2
    JuMP.@objective(m, Min, Common)
    _, Xn = Gurobi.GRBsetintparam(ge, "PoolSolutions", nX), [similar(Xl) for _=1:nX]
    Gurobi.GRBsetintparam(ge, "PoolSearchMode", 2)
    (; m, o, refd, refi, S, N, Θ, Xl, Xn, Bv)
end

function sub!(Ve, s, T, genv, t, LineˈP, F, WD, LD, GD, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=50, ζCost=150)
    ReserveK01v, LoadKref, times,Times, lines = Case2383.Reserve_Curve[reserve_type], Case2383.Load_Curve, t+1:t+T,t:t+T, keys(LineˈP)
    WDˈS, WDˈN, WDˈPMax                              = WD.S, WD.N, WD.PMax
    LDˈi, LDˈg, LDˈn, windMatˈs                      = LD.i, LD.g, LD.n, WDˈS[s]
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Qˈy, S, m = _0(), length(WDˈS), Settings.Model(genv); o,refi,refd = m.moi_backend,Ref{Cint}(),Ref{Cdouble}(); Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "Presolve", 2)
    egp         =  Dict(i    =>_0() for i=Times)
    rgp         = [Dict(i    =>_0() for i=Times) for z=1:4] # 1≤z≤4
    pfe         =  Dict((l,i)=>_0() for i=Times  for l=lines)
    JuMP.@variable(m, qy)
    x1, Xl = Dict{Tuple{Int, Int}, Union{Bool, JuMP.VariableRef}}(), Float64[]
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on, then don't have start-up cost
            if pAHˈzg + Pmin > Sd || any(==(false), view(xHˈzg, 1:UT))
                x1[z,g] = true
            else
                x1[z,g] = _variable(m); push!(Xl, 0)
            end
        else # was off
            if any(view(xHˈzg, 1:DT)) # cannot turn on, then don't have start-up cost
                x1[z,g] = false
            else
                x1[z,g] = also_u = _variable(m); push!(Xl, 0)
                # JuMP.add_to_expression!(Common, suC, also_u)
            end
        end
    end; N, Xl2, x_che = length(Xl), similar(Xl), similar(Xl) # x_che is locally static
    Ci = Cint[range(S; length=N); s-1] # [READ-ONLY] provide lag-cut to mst, not myself
    Cd = Vector{Cdouble}(undef, N+2) # [General Container, mainly to store lag-cut info]

    Cnt_2nd_Bin, x2 = 0, Dict{Tuple{Int,Int,Int}, Union{Bool,JuMP.VariableRef}}()
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        for Ofst = 1:T
            let i=t+Ofst, xm1=(Ofst>1 ? x2[z,g,i-1] : x1[z,g])
                if xm1 === true # was on
                    if any(==(false), view(xHˈzg, 1:UT-Ofst)) || pAHˈzg + Pmin > Sd + Rd * Ofst
                        x2[z,g,i] = true
                    else
                        x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                    end
                elseif xm1 === false # was off
                    if any(view(xHˈzg, 1:DT-Ofst))
                        x2[z,g,i] = false
                    else
                        x2[z,g,i] = also_u = _variable(m); Cnt_2nd_Bin += 1
                        JuMP.add_to_expression!(Qˈy, suC, also_u)
                    end
                else # is decision
                    x2[z,g,i] = x_2 = _variable(m); u_2 = _variable(m); Cnt_2nd_Bin += 2
                    JuMP.@constraint(m, u_2 >= x_2 - xm1)
                    JuMP.add_to_expression!(Qˈy, suC, u_2)
                end
            end
        end
    end; Bv = fill(Cchar('B'), N+Cnt_2nd_Bin)

    JuMP.@variables(m, begin
        0 ≤ p2[z=eachindex(GDˈn), g=eachindex(GDˈn[z])] # ✅ p := (Pmin)x + pA (But here only for pricing purpose, only for x1)
        0 ≤ a2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0 ≤ pe[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times] # effective power of generator, less than its actual output
        0           ≤  ϖ[z=eachindex(WDˈN),               Times] ≤ 1
        0           ≤  ζ[z=eachindex(LDˈn),               Times] ≤ 1
        -LineˈP[l]  ≤  pf[l=lines,i=Times]                       ≤ LineˈP[l]
    end);             

    #= Generators =# for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            Gwid, xHˈzg1 = PMax-Pmin, xH[z][g][1]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g] * linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
            x1ˈzg, p2ˈzg = x1[z,g], p2[z,g]
        end
        #= Uvx =# JuMP.@constraints(m, begin
            UT * (x1ˈzg - xHˈzg1) <= x1ˈzg + sum(x2[z,g,j] for j=range(t+1;length=UT-1))
            UT * (x2[z,g,t+1] - x1ˈzg) <= sum(x2[z,g,j] for j=range(t+1;length=UT))
            -DT * (x1ˈzg - xHˈzg1) <= 1-x1ˈzg + sum(1-x2[z,g,j] for j=range(t+1;length=DT-1))
            -DT * (x2[z,g,t+1] - x1ˈzg) <= sum(1-x2[z,g,j] for j=range(t+1;length=DT))
        end); for i=t+2:t+T
            kLe = t+T - i + 1
            kU, kD = min(kLe, UT), min(kLe, DT)
            JuMP.@constraints(m, begin
                kU * (x2[z,g,i] - x2[z,g,i-1]) <= sum(x2[z,g,j] for j=range(i;length=kU))
                -kD * (x2[z,g,i] - x2[z,g,i-1]) <= sum(1-x2[z,g,j] for j=range(i;length=kD))
            end)
        end
        #= Qˈy / pe =# JuMP.@constraint(m, p2ˈzg == Pmin * x1ˈzg + a2[z,g,t])
        JuMP.@constraint(m, pe[z,g,t] <= p2ˈzg)
        JuMP.add_to_expression!(Qˈy, linC, p2ˈzg); for i=times
            JuMP.add_to_expression!(Qˈy, linC * Pmin, x2[z,g,i])
            JuMP.add_to_expression!(Qˈy, linC, a2[z,g,i])
            JuMP.@constraint(m, pe[z,g,i] <= Pmin * x2[z,g,i] + a2[z,g,i])
        end
        #= Reserve =# 1≤z≤4 && for i=Times JuMP.add_to_expression!(rgp[z][i], ur[z,g,i]) end
        #= Energy & LineFlow =# for i=Times
            JuMP.add_to_expression!(egp[i], pe[z,g,i])
            for l=lines
                (pfl, Cnl) = (pfe[l,i], F[l, node])
                JuMP.add_to_expression!(pfl, Cnl, pe[z,g,i])
            end
        end
        #= Ramp Up/Dn =# MmRatio > 1.1 && JuMP.@constraints(m, begin
            a2[z,g,t] + ur[z,g,t] - pAHˈzg <=                   Ru * xHˈzg1      + (Su - Pmin) * (1 - xHˈzg1)
            a2[z,g,t+1] + ur[z,g,t+1] - a2[z,g,t] <=            Ru * x1ˈzg       + (Su - Pmin) * (1 - x1ˈzg)
            [i=t+2:t+T], a2[z,g,i] + ur[z,g,i] - a2[z,g,i-1] <= Ru * x2[z,g,i-1] + (Su - Pmin) * (1 - x2[z,g,i-1])
            pAHˈzg - a2[z,g,t]                              <= Rd * x1ˈzg       + (Sd - Pmin) * (1 - x1ˈzg)
            [i=t+1:t+T], a2[z,g,i-1] - a2[z,g,i]            <= Rd * x2[z,g,i]   + (Sd - Pmin) * (1 - x2[z,g,i])
        end)
        JuMP.@constraint(m, #= Intra-period =# a2[z,g,t] + ur[z,g,t] ≤ Gwid *     x1ˈzg)
        JuMP.@constraint(m, [i=times],         a2[z,g,i] + ur[z,g,i] ≤ Gwid * x2[z,g,i])
    end

    #= Wind/Load =# for i=Times
        egpˈi, ReserveK01 = egp[i], ReserveK01v[i]
        for (z,Ø)=enumerate(WDˈN)
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
        for (z,Ø)=enumerate(LDˈn)
            ζCostˈz, ζˈz = length(Ø)ζCost, ζ[z,i] # load shed
            JuMP.add_to_expression!(Qˈy, ζCostˈz, ζˈz)
            1≤z≤4 && (rgpˈz = rgp[z][i])
            for (g,node)=enumerate(Ø)
                Ty, ConstLoad = LDˈi[z][g], LDˈg[z][g]
                load = LoadKref[Ty][i] * ConstLoad
                if 1≤z≤4
                    # reserve requirements in isolated systems depends on overall demand [10.1109/TSG.2015.2469134]
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
    JuMP.@objective(m, Min, Qˈy) # start with this
    JuMP.@constraints(m, begin
        Qˈy == qy
        [z=1:4, i=Times],   rgp[z][i] ≥ 0 # Doesn't need relaxation
        [i=Times],          egp[i] == 0
        [l=lines,i=Times],   pfe[l,i] == pf[l,i]
    end)
    Ve[s] = (; m, o, refd, refi, Ci, Cd, Bv, Xl, Xl2, x_che, N)
end

end
