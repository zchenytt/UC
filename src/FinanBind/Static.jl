module Static
import JuMP, ..Case2383, ..Settings

_0() = JuMP.AffExpr(0.)
function Model2383!(m, t, F, liD, WD, LD, GN, GMax, Gmin;ReserveK2D=0.04, reserve_type=1)
    ReserveK01, LoadKref = Case2383.Reserve_Curve[reserve_type][t], Case2383.Load_Curve
    LDˈi, LDˈg, LDˈn = LD.i, LD.g, LD.n
    lines = keys(liD)
    pfe, egp, rgp = Dict(k=>_0() for k=lines), _0(), [_0() for z=1:4]
    JuMP.@variables(m, begin
        -liD[k] ≤ pf[k=lines] ≤ liD[k]
        x[z=eachindex(GN), g=eachindex(GN[z])], Bin # static problem, no turn-on/off actions
        0 ≤ pA[z=eachindex(GN), g=eachindex(GN[z])]
        0 ≤ ur[z=eachindex(GN), g=eachindex(GN[z])]
        0 ≤ ϖ[z=eachindex(WD.N)] ≤ 0
        0 ≤ ζ[z=eachindex(LDˈn)] ≤ 0
    end); JuMP.@constraint(m, [z=eachindex(GN), g=eachindex(GN[z])], 
        pA[z,g]+ur[z,g] ≤ (GMax[z][g]-Gmin[z][g])x[z,g])
    for (z,v)=enumerate(LDˈg)
        ζˈz = ζ[z] # load shed
        1≤z≤4 && (rˈz = rgp[z]) 
        for (g,p)=enumerate(v)
            Ty, n = LDˈi[z][g], LDˈn[z][g]
            load = LoadKref[Ty][t] * p
            JuMP.add_to_expression!(egp, -load)
            JuMP.add_to_expression!(egp, load, ζˈz)
            if 1≤z≤4
                upResDemand = ReserveK01 * ReserveK2D * p
                JuMP.add_to_expression!(rˈz, -upResDemand)
                JuMP.add_to_expression!(rˈz, upResDemand, ζˈz)
            end
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * -load)
                JuMP.add_to_expression!(pfl, Cnl * load, ζˈz)
            end
        end
    end
    for (z,v)=enumerate(WD.PMax)
        Kˈz, ϖˈz = WD.S[1][1, z], ϖ[z]
        for (g,pMa)=enumerate(v)
            pWind = Kˈz * pMa
            JuMP.add_to_expression!(egp, pWind)
            JuMP.add_to_expression!(egp, -pWind, ϖˈz)
            n = WD.N[z][g]
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * pWind)
                JuMP.add_to_expression!(pfl, Cnl * -pWind, ϖˈz)
            end
        end
    end
    for (z,v)=enumerate(Gmin)
        1≤z≤4 && (rˈz = rgp[z])
        for (g,pmin)=enumerate(v)
            1≤z≤4 && JuMP.add_to_expression!(rˈz, ur[z,g])
            JuMP.add_to_expression!(egp, pmin, x[z,g])
            JuMP.add_to_expression!(egp, pA[z,g])
            n = GN[z][g]
            for l=lines
                Cnl, pfl = F[l,n], pfe[l]
                JuMP.add_to_expression!(pfl, Cnl * pmin, x[z,g])
                JuMP.add_to_expression!(pfl, Cnl, pA[z,g])
            end
        end
    end
    JuMP.@constraints(m, begin
        egp == 0
        rgp .>= 0
        [k=lines], pfe[k] == pf[k]
    end)
    JuMP.@objective(m, Min, sum(rand(5.0:1e-4:8.0)i for i=x) + sum(rand(0.5:1e-4:1.5)i for i=pA))
end

function get_History_2383(t, F, Line, WD, LD, GD, genv)
    GDˈUT = GD.UT
    uL = max(maximum(maximum, GDˈUT),maximum(maximum, GD.DT))-1 # uH is not only for `u` but also for `v`
    uH =  [[falses(uL  ) for _=v] for v=GDˈUT]
    xH =  [[falses(uL+1) for _=v] for v=GDˈUT]
    pAH = [[NaN          for _=v] for v=GDˈUT]
    model = Settings.Model(genv)
    Model2383!(model, t, F, Line.P, WD, LD, GD.n, GD.pmax, GD.pmin)
    m = (m = model, o = model.moi_backend, refd = Ref{Cdouble}(), refi = Ref{Cint}())
    ter = Settings.opt_and_ter(m)
    ter == 2 || error("terminate code $ter")
    pA, x = m.m[:pA], m.m[:x]
    for (z,v)=enumerate(uH), (g,bv)=enumerate(v)
        pAH[z][g] = Settings.getxdblattrelement(m, pA[z,g], "X")
        if round(Bool, Settings.getxdblattrelement(m, x[z,g], "X"))
            if rand() ≤ .9
                j = rand(eachindex(bv))
                xH[z][g][1:j] .= bv[j] = true
            else
                xH[z][g] .= true
            end
        else
            if rand() ≤ .9
                j = rand(eachindex(bv))
                xH[z][g][j+1:end] .= true
            end
        end
    end
    xH,uH,pAH
end

end
