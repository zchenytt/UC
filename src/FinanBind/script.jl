# SEEDi = parse(Int, ARGS[1])
# printstyled("SEEDi = $SEEDi\n"; color=:magenta)

SEEDi = 2
import Random, Gurobi, JuMP, LinearAlgebra
SEED = hash(SEEDi); # rand(1:typemax(Int32));
Random.seed!(SEED);
t, T = 1, 23;
NP, NB = 126, 40; # number per batch, number of batches
S = (NB)NP; # number of scenes

include("src/Settings.jl");
include("src/Case2383.jl");
include("src/WindGen.jl");
include("src/Static.jl");
include("src/Uvx.jl");
include("src/Models.jl");
include("src/SBen.jl");

@time mst, sub, pas = SBen.prepare_models(S, t, T, NP);

inid = @timed SBen.ini_lb!(mst, sub, pas); # This operation causes main memory explosion e.g. 400+G/500G
# 74.628084699

lbv = Float64[];
timensinalgo3mst = Float64[];
timensinalgo3 = nothing;
mstreopttime = Float64[];
subtermigap = nothing;
algo3d = @timed SBen.trainsb_shortly(mst, sub, pas, NP, lbv, mstreopttime, subtermigap, timensinalgo3, timensinalgo3mst; Times = 11.);
# 7145.092644989
# () 446GiB total memory occupation

# eval_MIP_Runtime = (t^opt) = 1302.359650850296
# Sol=1: lb = 15503.470357491518, agap = 14.131738000907717, rgap = 0.0009106908344435977
# Sol=2: lb = 15503.470357491518, agap = 451.3995121021635, rgap = 0.028292271625632453

@time SBen._eval(mst, sub);
Non-default parameters:
Threads  1
PoolSolutions  2
PoolSearchMode  2
Optimize a model with 66745 rows, 11137 columns and 174849416 nonzeros (Min)

map(length, (timensinalgo3mst, lbv, mstreopttime))
