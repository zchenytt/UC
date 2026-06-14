# TODO create 2-master, one monotheta, one multitheta
import Random, Statistics, JLD2
get_free_memory() = Sys.free_memory() / 1024^3;
include("src/Settings.jl");
include("src/Case2383.jl");
include("src/WindGen.jl");
include("src/Static.jl");
include("src/Models.jl");
include("src/SBen.jl");
t, T = 1, 23;
NP = Sys.CPU_THREADS - 1;
S = 7000;
NcutWanted = 50000;

function run_one_test(SEEDi, t, T, NP, S, NcutWanted)
    printstyled("   julia> SEEDi=$SEEDi, S=$S, NP=$NP, memory=$(get_free_memory()) GiB\n"; color=208)
    Settings.printinfo()
    Random.seed!(hash(SEEDi))
    mst, sub = SBen.prepare_models(S, t, T, NP)
    @time "ini_lb!" SBen.ini_lb!(mst, sub, NP)
    # @time "pre_eval" SBen._eval(mst, sub, NP);
    sD = (msttxs, mstroiys, mstlbys, mstrotys, subtxs, subrotys, subgapys, subvioys, subVμys) = (Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])
    @time "Algo3" SBen.trainsb_shortly(mst, sub, NP, sD; NcutWanted = NcutWanted)
    JLD2.save(
        "algo3data$SEEDi.jld2", Dict(
        "msttxs"    => msttxs,
        "mstroiys"  => mstroiys,
        "mstlbys"   => mstlbys,
        "mstrotys"  => mstrotys,
        "subtxs"    => subtxs,
        "subrotys"  => subrotys,
        "subgapys"  => subgapys,
        "subvioys"  => subvioys,
        "subVμys"  => subVμys)
    )
    @time "final_eval" SBen._eval(mst, sub, NP)
    Vv = mst.Vv
    println("ending vioMean=", Statistics.mean(Vv), ", Extrema=", extrema(Vv), ", Median=", Statistics.median(Vv))
    printstyled("ending memory=", get_free_memory(), "GiB\n"; color=27)
end;
main(r, a...) = for SEEDi=r run_one_test(SEEDi, a...) end;
main(0:10, t, T, NP, S, NcutWanted)
