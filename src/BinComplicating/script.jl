SEEDi = parse(Int, ARGS[1])
printstyled("SEEDi = $SEEDi\n"; color=:magenta)
import Random, Gurobi, JuMP, LinearAlgebra
SEED = hash(SEEDi); # rand(1:typemax(Int32));
Random.seed!(SEED);
t, T = 1, 23;
NP, NB = 126, 48; # number per batch, number of batches
S = (NB)NP; # number of scenes

include("Settings.jl");
include("Case2383.jl");
include("WindGen.jl");
include("Static.jl");
include("Uvx.jl");
include("Models.jl");
include("SBen.jl");

@time mst, sub, pas = SBen.prepare_models(S, t, T, NP);
m1 = Sys.free_memory() / 1024^3;
@time SBen.ini_lb!(mst, sub, pas); # This operation causes main memory explosion e.g. 400+G/500G
m2 = Sys.free_memory() / 1024^3;
println("Mem used in ini_lb> $(m1-m2) GiB")
@time SBen._eval(mst, sub; ndisplay = 1);
@time SBen.trainsb_shortly(mst, sub, pas, NP; Times = 8.);
@time SBen._eval(mst, sub);
