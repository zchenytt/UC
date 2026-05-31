module Models
import ..Case2383, ..Uvx, ..Settings, JuMP, Gurobi
subTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, Ci::Vector{Int32}, Cd::Vector{Float64}, Bv::Vector{Int8}, Xl::Vector{Float64}, Xl2::Vector{Float64}, x_che::Vector{Float64}, N::Int}
pasTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refi::Base.RefValue{Int32}, refd::Base.RefValue{Float64}, VCi::Vector{Int32}, Cd::Vector{Float64}, QCi::Vector{Int32}, levelCCi::Vector{Int32}}

_variable(m) = JuMP.@variable(m, lower_bound=0, upper_bound=1)
function mst!(genv, WD, GD, uH, xH, pAH)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Common, S, m = _0(), length(WD.S), Settings.Model(genv)
    JuMP.@variable(m, common); o,refi,refd=m.moi_backend,Ref{Cint}(),Ref{Cdouble}(); ge=Gurobi.GRBgetenv(o)
    JuMP.@variable(m, θ[1:S]) # for `x` we may have FixedData or Decision; Since `u`'s decision is just === x, we only allocate the FixedData part of `u` 
    x1     = Dict{Tuple{Int, Int}, Union{Bool, JuMP.VariableRef}}()
    u1, Xl = Dict{Tuple{Int, Int},       Bool                   }(), Float64[]
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g] * linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on
            u1[z,g] = false
            if Sd < Pmin + pAHˈzg || any(==(0), view(xHˈzg, 1:UT))
                x1[z,g] = true # must keep; Since this is not a decision, it won't incur meaningful cost
            else
                x1[z,g] = _variable(m); push!(Xl, 0)
            end
        else # was off
            if any(view(xHˈzg, 1:DT))
                u1[z,g] = x1[z,g] = false
            else
                x1[z,g] = also_u = _variable(m); push!(Xl, 0)
                JuMP.add_to_expression!(Common, suC, also_u)
            end
        end
    end; N = length(Xl)
    VCi, Bv = Cint.(range(1+S; length=N)), fill(Cchar('B'), N)
    CCi, _, Cd = Cint.(0:N), JuMP.@constraint(m, [1:N+1], 0 == 0), fill(NaN, N+2) # Ci is READ-ONLY
                             JuMP.@constraint(m, Common == common)
    Θ = fill(-Inf, S); JuMP.@objective(m, Min, Common); nX = 2
    _, Xn = Gurobi.GRBsetintparam(ge, "PoolSolutions", nX), [similar(Xl) for _=1:nX]
    Gurobi.GRBsetintparam(ge, "PoolSearchMode", 2)
    (; m, o, refd, refi, S, N, Θ, Xl, Xn, Bv, CCi, VCi, Cd)
end

function sub!(Ve, s, T, genv, t, LineˈP, F, WD, LD, GD, uH, xH, pAH; ReserveK2D=0.04, reserve_type=1, ϖCost=130, ζCost=150)
    ReserveK01v, LoadKref, times,Times, lines = Case2383.Reserve_Curve[reserve_type], Case2383.Load_Curve, t+1:t+T,t:t+T, keys(LineˈP)
    WDˈS, WDˈN, WDˈPMax                              = WD.S, WD.N, WD.PMax
    LDˈi, LDˈg, LDˈn, windMatˈs                      = LD.i, LD.g, LD.n, WDˈS[s]
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Common, Qˈy, S, m = _0(), _0(), length(WDˈS), Settings.Model(genv); o,refi,refd = m.moi_backend,Ref{Cint}(),Ref{Cdouble}()
    egp         =  Dict(i    =>_0() for i=Times)
    rgp         = [Dict(i    =>_0() for i=Times) for z=1:4] # 1≤z≤4
    pfe         =  Dict((l,i)=>_0() for i=Times  for l=lines)
    Gurobi.GRBsetintparam(Gurobi.GRBgetenv(o), "Presolve", 2)
    JuMP.@variable(m, qy) # index 0

    x1     = Dict{Tuple{Int, Int}, Union{Bool, JuMP.VariableRef}}()
    u1, Xl = Dict{Tuple{Int, Int},       Bool                   }(), Float64[] # Xl == pi_hat
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            xHˈzg, Gwid = xH[z][g], PMax-Pmin
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g] * linC
        end
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        if xHˈzg[1] # was on
            u1[z,g] = false
            if Sd < Pmin + pAHˈzg || any(==(0), view(xHˈzg, 1:UT))
                x1[z,g] = true
            else
                x1[z,g] = _variable(m); push!(Xl, 0)
            end
        else # was off
            if any(view(xHˈzg, 1:DT))
                u1[z,g] = x1[z,g] = false
            else 
                x1[z,g] = also_u = _variable(m); push!(Xl, 0)
            end
        end
    end; N, Xl2, x_che = length(Xl), similar(Xl), similar(Xl) # x_che is locally static
    Ci = Cint[range(1+S; length=N); s] # [READ-ONLY] provide lag-cut to mst, not myself
    Cd = Vector{Cdouble}(undef, N+2) # [General Container, mainly to store lag-cut info]

    Cnt_2nd_Bin,x2,u2 = 0,Dict{Tuple{Int, Int, Int}, Union{Bool, JuMP.VariableRef}}(),Dict{Tuple{Int, Int, Int}, Union{Bool, JuMP.VariableRef}}()
    for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        xHˈzg, UT, DT = @inbounds (xH[z][g], GDˈUT[z][g], GDˈDT[z][g])
        let Ofst=1, i=t+Ofst
            if x1[z,g] === true # was on
                u2[z,g,i] = false
                if any(==(0), view(xHˈzg, 1:UT-Ofst))
                    x2[z,g,i] = true
                else
                    x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                end
            elseif x1[z,g] === false # was off
                if any(view(xHˈzg, 1:DT-Ofst))
                    u2[z,g,i] = x2[z,g,i] = false
                else
                    u2[z,g,i] = x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                end
            else # is decision
                u2[z,g,i],x2[z,g,i] = _variable(m),_variable(m); Cnt_2nd_Bin += 2
            end
        end
        for Ofst = 2:T
            i=t+Ofst
            if x2[z,g,i-1] === true # was on
                u2[z,g,i] = false
                if any(==(0), view(xHˈzg, 1:UT-Ofst))
                    x2[z,g,i] = true
                else
                    x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                end
            elseif x2[z,g,i-1] === false # was off
                if any(view(xHˈzg, 1:DT-Ofst))
                    u2[z,g,i] = x2[z,g,i] = false
                else
                    u2[z,g,i] = x2[z,g,i] = _variable(m); Cnt_2nd_Bin += 1
                end
            else # is decision
                u2[z,g,i],x2[z,g,i] = _variable(m),_variable(m); Cnt_2nd_Bin += 2
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
            Gwid, uHˈzg, xHˈzg = PMax-Pmin, uH[z][g], xH[z][g]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s, suC = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g], GDˈsci[z][g] * linC
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
            x1ˈzg = x1[z,g]; u1ˈzg = haskey(u1, (z,g)) ? u1[z,g] : x1ˈzg
        end
        # if UT == 1, use these coefficients
        T12 = Sd-PMax
        T13 = -max(Sd-Su,0.)
        T15 = -max(Su-Sd,0.)
        T16 = Su-PMax
        # Ramp constr, use these
        R1 = Su-Pmin
        R2 = Sd-Pmin-Rd
        for i=times JuMP.@constraints(m, begin
            Uvx.vge0!(_0(), i, t, z, g, x1ˈzg, x2) ≥ -u2[z,g,i]
            Uvx.u!(_0(),UT,i,t,z,g,uHˈzg,u1ˈzg,u2) ≤ x2[z,g,i]
            Uvx.v!(_0(),DT,i,t,z,g,xHˈzg,x1ˈzg,x2,uHˈzg,u1ˈzg,u2) ≤ 1
        end) end # all the logic constraints
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
            JuMP.add_to_expression!(Qˈy, suC,  u2[z,g,i]) # only for times
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
    JuMP.@objective(m, Min, Qˈy) # start with this
    JuMP.@constraints(m, begin
        Qˈy == qy
        [i=Times],          egp[i] == 0
        [z=1:4, i=Times],   rgp[z][i] ≥ 0
        [l=lines,i=Times],  pfe[l,i] == pf[l,i]
    end)
    Ve[s] = (; m, o, refd, refi, Ci, Cd, Bv, Xl, Xl2, x_che, N)
end

_0() = JuMP.AffExpr(0.)
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end

end

