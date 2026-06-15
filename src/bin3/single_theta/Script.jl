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

function run_one_test(SEEDi, t, T, NP, S)
    printstyled("   julia> SEEDi=$SEEDi, S=$S, NP=$NP, memory=$(get_free_memory()) GiB\n"; color=208)
    Settings.printinfo()
    Random.seed!(hash(SEEDi))
    mst, sub = SBen.prepare_models(S, t, T, NP)
    @time "ini_lb!" SBen.ini_lb!(mst, sub, NP)
    # TODO write the main algorithm with single-theta
    @time "Algo3" SBen.trainsb_shortly(mst, sub, NP, nothing)
    @time "final_eval" SBen._eval(mst, sub, NP)
    printstyled("ending memory=", get_free_memory(), "GiB\n"; color=27)
end;
function main(t, T, NP)
    run_one_test(4, t, T, NP, 4000)
    run_one_test(5, t, T, NP, 5000)
    run_one_test(6, t, T, NP, 6000)
    run_one_test(7, t, T, NP, 7000)
end
main(t, T, NP)
