"""
For the finantial-binding variant:

TODO you need to distinguish the t=1 cost vs common cost,
in order to compare with the previous bin_complicating case
"""
module Models
import ..Case2383, ..Uvx, ..Settings, JuMP, Gurobi
subTy = @NamedTuple{m::JuMP.Model, o::Gurobi.Optimizer, refd::Base.RefValue{Float64}, refi::Base.RefValue{Int32}, Ci::Vector{Int32}, Cd::Vector{Float64}, Xl::Vector{Float64}, Xl2::Vector{Float64}, x_che::Vector{Float64}, N::Int64}

_variable(m) = JuMP.@variable(m, lower_bound=0, upper_bound=1)
function mst!(genv, t, T, WD, GD, uH, xH, pAH)
    GDˈn, GDˈUT, GDˈDT, GDˈpmax, GDˈMmRatio, GDˈpmin = GD.n, GD.UT, GD.DT, GD.pmax, GD.MmRatio, GD.pmin
    GDˈrui, GDˈrdi, GDˈr2s, GDˈsci, GDˈC             = GD.rui, GD.rdi, GD.r2s, GD.sci, GD.C
    Common, S, m = _0(), length(WD.S), Settings.Model(genv); o,refi,refd=m.moi_backend,Ref{Cint}(),Ref{Cdouble}()
    JuMP.@variable(m, common);      JuMP.@variable(m, θ[1:S])
    JuMP.@variable(m, 0 <= x[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=t:t+T]); N = length(x)
    JuMP.@variable(m, 0 <= u[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=t:t+T])
    for z=eachindex(GDˈn), g=eachindex(GDˈn[z])
        PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
        Gwid, xHˈzg = PMax-Pmin, xH[z][g]; UT, DT = GDˈUT[z][g], GDˈDT[z][g]
        rui, rdi, r2s = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g]
        Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        let Ofst=0, i=t+Ofst, uzgi=u[z,g,i], xzgi=x[z,g,i]
            if xHˈzg[1] # the unit was ON
                JuMP.set_upper_bound(uzgi, 0)
                if any(==(0), view(xHˈzg, 1:UT-Ofst)) || pAHˈzg+Pmin > Rd*Ofst + Sd
                    JuMP.set_lower_bound(xzgi, 1)
                end
            else # the unit was OFF
                if any(view(xHˈzg, 1:DT-Ofst))
                    JuMP.set_upper_bound(uzgi, 0); JuMP.set_upper_bound(xzgi, 0)
                else
                    JuMP.@constraint(m, uzgi == xzgi)
                end
            end
        end
        for Ofst = 1:T
            let i=t+Ofst, uzgi=u[z,g,i], xzgi=x[z,g,i]
                if JuMP.lower_bound(x[z,g,i-1]) == 1 # the unit was fixed at ON
                    JuMP.set_upper_bound(uzgi, 0)
                    if any(==(0), view(xHˈzg, 1:UT-Ofst)) || pAHˈzg+Pmin > Rd*Ofst + Sd
                        JuMP.set_lower_bound(xzgi, 1)
                    end
                elseif JuMP.has_upper_bound(x[z,g,i-1]) && JuMP.upper_bound(x[z,g,i-1]) == 0
                    if any(view(xHˈzg, 1:DT-Ofst))
                        JuMP.set_upper_bound(uzgi, 0); JuMP.set_upper_bound(xzgi, 0)
                    else
                        JuMP.@constraint(m, uzgi == xzgi)
                    end
                end
            end
        end
    end
    JuMP.@constraint(m, [z=eachindex(GDˈn), g=eachindex(GDˈn[z])], u[z,g,t]-x[z,g,t]+xH[z][g][1] >= 0)
    JuMP.@constraint(m, [z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=t+1:t+T], u[z,g,i]-x[z,g,i]+x[z,g,i-1] >= 0)
    for z=eachindex(GDˈn), g=eachindex(GDˈn[z])
        u1ˈzg, x1ˈzg, suC = u[z,g,t], x[z,g,t], GDˈC[z][g]GDˈsci[z][g]
        uHˈzg, xHˈzg, UT, DT = uH[z][g], xH[z][g], GDˈUT[z][g], GDˈDT[z][g]
        for i=t:t+T
            JuMP.add_to_expression!(Common, suC, u[z,g,i])
            JuMP.@constraint(m, Uvx.u!(_0(),UT,i,t,z,g,uHˈzg,u1ˈzg,u) ≤ x[z,g,i])
            JuMP.@constraint(m, Uvx.v!(_0(),DT,i,t,z,g,xHˈzg,x1ˈzg,x,uHˈzg,u1ˈzg,u) ≤ 1)
        end
    end
    Bv, Xl, Θ = fill(Cchar('B'), N), fill(NaN, N), fill(-Inf, S)
    JuMP.@constraint(m, Common == common)
    JuMP.@objective(m, Min, Common); nX = 2; ge=Gurobi.GRBgetenv(o)
    _, Xn = Gurobi.GRBsetintparam(ge, "PoolSolutions", nX), [similar(Xl) for _=1:nX]
    Gurobi.GRBsetintparam(ge, "PoolSearchMode", 2)
    (; m, o, refd, refi, S, N, Θ, Xl, Xn, Bv)
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
    JuMP.@variable(m, 0 <= x2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]); N = length(x2)
    Xl = fill(NaN, N); Xl2 = similar(Xl); x_che = similar(Xl) # x_che is locally static
    Ci = Cint[range(1+S; length=N); s] # [READ-ONLY] provide lag-cut to mst, not myself
    Cd = Vector{Cdouble}(undef, N+2) # [General Container, mainly to store lag-cut info]
    JuMP.@variables(m, begin
        0 ≤ a2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0 ≤ p2[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times] # ✅ p := (Pmin)x + pA (But here only for pricing purpose)
        0 ≤ ur[z=eachindex(GDˈn), g=eachindex(GDˈn[z]), i=Times]
        0           ≤  ϖ[z=eachindex(WDˈN),               Times] ≤ 1
        0           ≤  ζ[z=eachindex(LDˈn),               Times] ≤ 1
        -LineˈP[l]  ≤ pf[l=lines,                         Times] ≤ LineˈP[l]
    end)
    slk = JuMP.VariableRef[]
    #= Generators =# for (z,Ø)=enumerate(GDˈn), (g,node)=enumerate(Ø) # Put Generator at outer layer, so technic params are calculated once
        @inbounds begin
            PMax, MmRatio, Pmin, pAHˈzg = GDˈpmax[z][g], GDˈMmRatio[z][g], GDˈpmin[z][g], pAH[z][g]
            Gwid, xHˈzg = PMax-Pmin, xH[z][g]
            UT, DT, linC = GDˈUT[z][g], GDˈDT[z][g], GDˈC[z][g]
            rui, rdi, r2s = GDˈrui[z][g], GDˈrdi[z][g], GDˈr2s[z][g]
            Ru, Rd, Su, Sd = _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
        end
        #= global-wide =# for i=Times
            slk1 = JuMP.@variable(m, lower_bound = 0)
            1≤z≤4 && JuMP.add_to_expression!(rgp[z][i], ur[z,g,i])
            JuMP.add_to_expression!(egp[i], Pmin, x2[z,g,i])
            JuMP.add_to_expression!(egp[i], a2[z,g,i])
            JuMP.add_to_expression!(egp[i], -1., slk1)
            for l=lines
                (pfl, Cnl) = (pfe[l,i], F[l, node])
                JuMP.add_to_expression!(pfl, Cnl * Pmin, x2[z,g,i])
                JuMP.add_to_expression!(pfl, Cnl, a2[z,g,i])
                JuMP.add_to_expression!(pfl, -Cnl, slk1)
            end
            JuMP.@constraint(m, p2[z,g,i] == Pmin * x2[z,g,i] + a2[z,g,i])
            push!(slk, slk1)
            JuMP.add_to_expression!(Qˈy, linC, p2[z,g,i])
        end
        for i=Times
            urˈzg, pAˈzg, xˈzg = ur[z,g,i], a2[z,g,i], x2[z,g,i]
            JuMP.@constraint(m, urˈzg+pAˈzg ≤ Gwid*xˈzg) # UpperLimit Constr
            if MmRatio > 1.1
                #= Normal Ramp Constr =# if i == t
                    JuMP.@constraint(m, urˈzg+pAˈzg - pAHˈzg ≤ Ru * xHˈzg[1] + (Su - Pmin) * (1 - xHˈzg[1]))
                    JuMP.@constraint(m, pAHˈzg - pAˈzg ≤ Rd * xˈzg + (Sd - Pmin) * (1 - xˈzg))
                else
                    JuMP.@constraint(m, urˈzg+pAˈzg - a2[z,g,i-1] ≤ Ru * x2[z,g,i-1] + (Su - Pmin) * (1 - x2[z,g,i-1]))
                    JuMP.@constraint(m, a2[z,g,i-1] - pAˈzg ≤ Rd * xˈzg + (Sd - Pmin) * (1 - xˈzg))                    
                end
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
    JuMP.@objective(m, Min, Qˈy + 7000 * sum(slk)) # start with this
    JuMP.@constraints(m, begin
        Qˈy == qy
        [i=Times],          egp[i] == 0
        [z=1:4, i=Times],   rgp[z][i] ≥ 0
        [l=lines,i=Times],  pfe[l,i] == pf[l,i]
    end)
    Ve[s] = (; m, o, refd, refi, Ci, Cd, Xl, Xl2, x_che, N)
end

_0() = JuMP.AffExpr(0.)
function _RURDSUSD(Gwid, rui, rdi, r2s, Pmin, PMax)
    Ru, Rd = Gwid/rui, Gwid/rdi
    Su, Sd = clamp(r2s * Ru, Pmin, PMax), clamp(r2s * Rd, Pmin, PMax)
    Ru, Rd, Su, Sd
end

end

