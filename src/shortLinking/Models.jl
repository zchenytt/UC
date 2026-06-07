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
function mst!(genv, WD, GD, xH, pAH) # It should be reoptimized only if there was separating cut added!
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Common, S, m = _0(), length(WD.S), Settings.Model(genv)
    o,refi,refd=m.moi_backend,Ref{Cint}(),Ref{Cdouble}(); ge=Gurobi.GRBgetenv(o)
    # theta: range(0;length=S); Xl starts at S;
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
            if Pmin + pAHˈzg > Sd || any(==(0), view(xHˈzg, 1:UT))
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
    end; N = length(Xl); Bv = fill(Cchar('B'), N); Θ = fill(-Inf, S); nX = 2
    JuMP.@objective(m, Min, Common)
    _, Xn = Gurobi.GRBsetintparam(ge, "PoolSolutions", nX), [similar(Xl) for _=1:nX]
    Gurobi.GRBsetintparam(ge, "PoolSearchMode", 2)
    (; m, o, refd, refi, S, N, Θ, Xl, Xn, Bv)
end

function sub!(Ve, s, T, genv, t, LineˈP, F, WD, LD, GD, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=130, ζCost=150)
    ReserveK01v, LoadKref, times,Times, lines = Case2383.Reserve_Curve[reserve_type], Case2383.Load_Curve, t+1:t+T,t:t+T, keys(LineˈP)
    WDˈS, WDˈN, WDˈPMax                              = WD.S, WD.N, WD.PMax
    LDˈi, LDˈg, LDˈn, windMatˈs                      = LD.i, LD.g, LD.n, WDˈS[s]
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Qˈy, S, m = _0(), length(WDˈS), Settings.Model(genv); o,refi,refd = m.moi_backend,Ref{Cint}(),Ref{Cdouble}()
    egp         =  Dict(i    =>_0() for i=Times)
    rgp         = [Dict(i    =>_0() for i=Times) for z=1:4] # 1≤z≤4
    pfe         =  Dict((l,i)=>_0() for i=Times  for l=lines)
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "Presolve", 2)
    JuMP.@variable(m, qy) # index 0
    x1, Xl = Dict{Tuple{Int, Int}, Union{Bool, JuMP.VariableRef}}(), Float64[] # Xl == pi_hat
    u1 = Dict{Tuple{Int, Int}, Bool}()
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on, then don't have start-up cost
            u1[z,g] = false
            if pAHˈzg + Pmin > Sd || any(==(0), view(xHˈzg, 1:UT))
                x1[z,g] = true
            else
                x1[z,g] = _variable(m); push!(Xl, 0)
            end
        else # was off
            if any(view(xHˈzg, 1:DT)) # cannot turn on, then don't have start-up cost
                u1[z,g] = x1[z,g] = false
            else
                x1[z,g] = also_u = _variable(m); push!(Xl, 0)
                # JuMP.add_to_expression!(Common, suC, also_u)
            end
        end
    end; N, Xl2, x_che = length(Xl), similar(Xl), similar(Xl) # x_che is locally static
    Ci = Cint[range(S; length=N); s-1] # [READ-ONLY] provide lag-cut to mst, not myself
    Cd = Vector{Cdouble}(undef, N+2) # [General Container, mainly to store lag-cut info]

    Cnt_2nd_Bin, x2 = 0, Dict{Tuple{Int,Int,Int}, Union{Bool,JuMP.VariableRef}}()
    u2 = Dict{Tuple{Int,Int,Int}, Union{Bool,JuMP.VariableRef}}() # needed as 2nd-stage cost includes it
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø)
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g]linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        let Ofst=1, i=t+Ofst
            if x1[z,g] === true # was on
                u2[z,g,i] = false
                if any(==(0), view(xHˈzg, 1:UT-Ofst)) || pAHˈzg + Pmin > Sd + Rd * Ofst
                    x2[z,g,i] = true
                else
                    x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                end
            elseif x1[z,g] === false # was off
                if any(view(xHˈzg, 1:DT-Ofst))
                    u2[z,g,i] = x2[z,g,i] = false
                else
                    u2[z,g,i] = x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                    JuMP.add_to_expression!(Qˈy, suC, u2[z,g,i])
                end
            else # is decision
                u2[z,g,i],x2[z,g,i] = _variable(m),_variable(m); Cnt_2nd_Bin += 2
                JuMP.add_to_expression!(Qˈy, suC, u2[z,g,i])
            end
        end
        for Ofst = 2:T
            i=t+Ofst
            if x2[z,g,i-1] === true # was on
                u2[z,g,i] = false
                if any(==(0), view(xHˈzg, 1:UT-Ofst)) || pAHˈzg + Pmin > Sd + Rd * Ofst
                    x2[z,g,i] = true
                else
                    x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                end
            elseif x2[z,g,i-1] === false # was off
                if any(view(xHˈzg, 1:DT-Ofst))
                    u2[z,g,i] = x2[z,g,i] = false
                else
                    u2[z,g,i] = x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                    JuMP.add_to_expression!(Qˈy, suC, u2[z,g,i])
                end
            else # is decision
                u2[z,g,i],x2[z,g,i] = _variable(m),_variable(m); Cnt_2nd_Bin += 2
                JuMP.add_to_expression!(Qˈy, suC, u2[z,g,i])
            end
        end
    end; Bv = fill(Cchar('B'), N + Cnt_2nd_Bin)

    JuMP.@variables(m, begin
        0 ≤ a2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0 ≤ p2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times] # ✅ p := (Pmin)x + pA (But here only for pricing purpose)
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0           ≤  ϖ[z=eachindex(WDˈN),               Times] ≤ 1
        0           ≤  ζ[z=eachindex(LDˈn),               Times] ≤ 1
        -LineˈP[l]  ≤ pf[l=lines,                         Times] ≤ LineˈP[l]
    end);             

    #= Generators =# for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            Gwid, xHˈzg = PMax-Pmin, xH[z][g]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g] * linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
            x1ˈzg = x1[z,g]; u1ˈzg = haskey(u1, (z,g)) ? false : x1ˈzg
        end
        T12 = Sd-PMax
        T13 = -max(Sd-Su,0.)
        T15 = -max(Su-Sd,0.)
        T16 = Su-PMax
        R1 = Su-Pmin
        R2 = Sd-Pmin-Rd
        begin # Uvx-polytope
            JuMP.@constraint(m,            u2[z,g,t+1]-x2[z,g,t+1]+x1ˈzg >= 0)
            JuMP.@constraint(m, [i=t+2:t+T], u2[z,g,i]-x2[z,g,i]+x2[z,g,i-1] >= 0)
            for i=times
                h = _0()
                for j = range(i; step=-1, length=UT)
                    if j > t
                        JuMP.add_to_expression!(h, u2[z,g,j])
                    elseif t > j
                        k=t-j; JuMP.add_to_expression!(h, max(0, xHˈzg[k] - xHˈzg[k+1]))
                    else
                        JuMP.add_to_expression!(h, u1ˈzg)
                    end
                end
                JuMP.@constraint(m, h <= x2[z,g,i]); h = _0(); j = i
                for _ = 1:DT
                    if j > t
                        JuMP.add_to_expression!(h, u2[z,g,j])
                    elseif t > j
                        k=t-j; JuMP.add_to_expression!(h, max(0, xHˈzg[k] - xHˈzg[k+1]))
                    else
                        JuMP.add_to_expression!(h, u1ˈzg)
                    end
                    j -= 1
                end
                if j > t
                    JuMP.add_to_expression!(h, x2[z,g,j])
                elseif t > j
                    k=t-j; JuMP.add_to_expression!(h, xHˈzg[k])
                else
                    JuMP.add_to_expression!(h, x1ˈzg)
                end
                JuMP.@constraint(m, h <= 1)
            end
        end
        #= global-wide =# for i=times
            1≤z≤4 && JuMP.add_to_expression!(rgp[z][i], ur[z,g,i])
            JuMP.add_to_expression!(egp[i], Pmin, x2[z,g,i])
            JuMP.add_to_expression!(egp[i], a2[z,g,i])
            for l=lines
                (pfl, Cnl) = (pfe[l,i], F[l, node])
                JuMP.add_to_expression!(pfl, Cnl * Pmin, x2[z,g,i])
                JuMP.add_to_expression!(pfl, Cnl,       a2[z,g,i])
            end
            JuMP.@constraint(m, p2[z,g,i] == Pmin * x2[z,g,i] + a2[z,g,i])
            JuMP.add_to_expression!(Qˈy, linC, p2[z,g,i])
        end
        #= global-wide =# let i = t
            1≤z≤4 && JuMP.add_to_expression!(rgp[z][i], ur[z,g,i]) 
            JuMP.add_to_expression!(egp[i], Pmin, x1ˈzg)
            JuMP.add_to_expression!(egp[i], a2[z,g,i])
            for l=lines
                (pfl, Cnl) = (pfe[l,i], F[l, node])
                JuMP.add_to_expression!(pfl, Cnl * Pmin, x1ˈzg)
                JuMP.add_to_expression!(pfl, Cnl,       a2[z,g,i])
            end
            JuMP.@constraint(m, p2[z,g,i] == Pmin * x1ˈzg + a2[z,g,i])
            JuMP.add_to_expression!(Qˈy, linC, p2[z,g,i])
        end
        let i = t
            if MmRatio > 1.1
                JuMP.@constraint(m, ur[z,g,i] + a2[z,g,i] ≤ Gwid * x1ˈzg + T16 * u1ˈzg)
            else
                JuMP.@constraint(m, ur[z,g,i] + a2[z,g,i] ≤ Gwid * x1ˈzg)
            end
        end
        for i=times
            urˈzg, pAˈzg, xˈzg, uˈzg = ur[z,g,i], a2[z,g,i], x2[z,g,i], u2[z,g,i]
            if MmRatio > 1.1
                #= UpperLimit Constr =# if i == t+T
                    JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg + T16 * uˈzg) # This is THE classic
                elseif UT > 1
                    JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg + T12 * (u2[z,g,i+1]-x2[z,g,i+1]+xˈzg) + T16 * uˈzg)
                else
                    JuMP.@constraints(m, begin
                        urˈzg+pAˈzg ≤ Gwid*xˈzg + T12 * (u2[z,g,i+1]-x2[z,g,i+1]+xˈzg) + T13 * uˈzg
                        urˈzg+pAˈzg ≤ Gwid*xˈzg + T15 * (u2[z,g,i+1]-x2[z,g,i+1]+xˈzg) + T16 * uˈzg
                    end)
                end
                #= Normal Ramp Constr =# if i == t+1
                    JuMP.@constraints(m, begin
                        urˈzg+pAˈzg - a2[z,g,i-1] ≤ Ru * x1ˈzg + R1 *  uˈzg
                        a2[z,g,i-1] - pAˈzg       ≤ Rd * x1ˈzg + R2 * (uˈzg - xˈzg + x1ˈzg)
                    end)
                else
                    JuMP.@constraints(m, begin
                        urˈzg+pAˈzg - a2[z,g,i-1] ≤ Ru * x2[z,g,i-1] + R1 *  uˈzg
                        a2[z,g,i-1] - pAˈzg       ≤ Rd * x2[z,g,i-1] + R2 * (uˈzg - xˈzg + x2[z,g,i-1])
                    end)
                end
            else
                JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg)
            end
        end
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
    JuMP.@variable(m, surp[i=Times] >= 0); for i=surp JuMP.add_to_expression!(Qˈy, 3e3, i) end
    JuMP.@objective(m, Min, Qˈy) # start with this
    JuMP.@constraints(m, begin
        Qˈy == qy
        [i=Times],          egp[i] == surp[i]
        [z=1:4, i=Times],   rgp[z][i] ≥ 0
        [l=lines,i=Times],  pfe[l,i] == pf[l,i]
    end)
    Ve[s] = (; m, o, refd, refi, Ci, Cd, Bv, Xl, Xl2, x_che, N)
end

end