# SEEDi = parse(Int, ARGS[1])


import Random, Gurobi, JuMP, LinearAlgebra, JLD2
include("src/Settings.jl");
include("src/Case2383.jl");
include("src/WindGen.jl");
include("src/Static.jl");
include("src/Models.jl");
include("src/Models.jl");
include("src/SBen.jl");

t, T = 1, 23;
NP, NB = 254, 24; # number per batch, number of batches
S = (NB)NP; # number of scenes
ROT = 3e3;

get_free_memory() = Sys.free_memory() / 1024^3;
function run_one_test(SEEDi, t, T, NP, S, ROT)
    printstyled("   julia> SEEDi = $SEEDi, memory = $(get_free_memory()) GiB\n"; color=:cyan)
    Random.seed!(hash(SEEDi))
    m2 = get_free_memory();
    @time "JuMP_Model" mst, sub = SBen.prepare_models(S, t, T, NP);
    m1 = get_free_memory(); printstyled("Mem used in JuMP Model> $(m2-m1) GiB\n"; color = :cyan)
    @time "ini_lb!" SBen.ini_lb!(mst, sub); # This operation causes main memory explosion e.g. 400+G/500G
    m2 = get_free_memory(); printstyled("Mem used in ini_lb> $(m1-m2) GiB\n"; color = :cyan)
    # @time "pre_eval" SBen._eval(mst, sub; ndisplay = 1);
    msttxs,mstrotys,mstroiys,mstlbys,mstncys = Float64[],Float64[],Float64[],Float64[],Int[]
    subtxs,subgapys = Float64[],Float64[]
    sampleData = (msttxs,mstrotys,mstroiys,mstlbys,mstncys,subtxs,subgapys)
    @time "Algo3" SBen.trainsb_shortly(mst, sub, NP, sampleData; ROT = ROT);
    JLD2.save(
        "algo3data$SEEDi.jld2", Dict(
        "msttxs"    => msttxs,
        "mstrotys"  => mstrotys,
        "mstroiys"  => mstroiys,
        "mstlbys"   => mstlbys,
        "mstncys"   => mstncys,
        "subtxs"    => subtxs,
        "subgapys"  => subgapys)
    ); printstyled("algo3data$SEEDi.jld2 saved.\n"; color = 27)
    @time "final_eval" SBen._eval(mst, sub);
    printstyled("ending memory = $(get_free_memory()) GiB\n"; color = :cyan)
end;
main(r, a...) = for SEEDi=r run_one_test(SEEDi, a...) end;

main(0:10, t, T, NP, S, ROT)



# amd@amd:~/julia_projects/uc/Uc2$ julia --project=. --threads=254,2 Script.jl
#    julia> SEEDi = 0, memory = 311.6371383666992 GiB
# s= 1394
# s= 1956
# s= 2759
# s= 3232
# s= 4177
# JuMP_Model: 1310.393737 seconds (85.14 G allocations: 2.463 TiB, 80.55% gc time, 105 lock conflicts, 106.07% compilation time: <1% of which was recompilation)
# Mem used in JuMP Model> -58.3213996887207 GiB
# lb = 3.05379e+03
# ini_lb!: 473.155417 seconds (661.70 k allocations: 33.227 MiB, 63 lock conflicts, 5.65% compilation time)
# Mem used in ini_lb> 251.66358184814453 GiB
# lb = 3.05379e+03
# s=3848, vio=8.500e+02
# s=1748, vio=4.261e+02
# s=2622, vio=1.535e+03
# s=2627, vio=9.142e+02
# s=5971, vio=2.629e+02
# s=3198, vio=9.958e+01
# ȷ=4, ro=1.7e-01, rs=9.7e+02, ri=2.6e+01
# s=3246, vio=4.608e+01
# s=3963, vio=3.384e+01
# s=4336, vio=1.296e+02
# s=4095, vio=6.181e+01
# s=1378, vio=1.356e+01
# ȷ=9, ro=4.0e-01, rs=2.4e+03, ri=2.0e+01
# s=2957, vio=1.984e+00
# s=1093, vio=5.403e-01
# lb = 4.55360e+03
# Algo3: 4458.906601 seconds (1.05 M allocations: 58.625 MiB, 3949 lock conflicts, 0.37% compilation time)
# Set parameter OutputFlag to value 1
# Gurobi Optimizer version 13.0.1 build v13.0.1rc0 (linux64gpu - "Ubuntu 24.04.4 LTS")

# CPU model: AMD EPYC 7763 64-Core Processor, instruction set [SSE2|AVX|AVX2]
# Thread count: 128 physical cores, 256 logical processors, using up to 1 threads

# GPU model: NVIDIA RTX A6000, CUDA compute version 8.6, NVIDIA driver compatible with CUDA version 13

# Non-default parameters:
# Threads  1
# PoolSolutions  2
# PoolSearchMode  2

# Optimize a model with 35426 rows, 6155 columns and 2125505 nonzeros (Min)
# Model fingerprint: 0x82fe7e24
# Model has 6155 linear objective coefficients
# Variable types: 6096 continuous, 59 integer (59 binary)
# Coefficient statistics:
#   Matrix range     [2e-05, 4e+04]
#   Objective range  [2e-04, 1e+02]
#   Bounds range     [1e+00, 1e+00]
#   RHS range        [4e+02, 7e+03]

# Found heuristic solution: objective 5156.5656466
# Found heuristic solution: objective 5156.5656466
# Presolve added 0 rows and 9 columns
# Presolve removed 62 rows and 0 columns
# Presolve time: 2.78s
# Presolved: 35364 rows, 6164 columns, 1750944 nonzeros
# Variable types: 6100 continuous, 64 integer (59 binary)
# Found heuristic solution: objective 5126.6984606
# Root relaxation presolved: 35364 rows, 6164 columns, 1750944 nonzeros


# Root simplex log...

# Iteration    Objective       Primal Inf.    Dual Inf.      Time
#     3779    4.4879926e+03   2.920135e+05   0.000000e+00      5s
#     9264    4.5527980e+03   5.813785e+03   0.000000e+00     10s
#    11893    4.5534968e+03   2.867177e+03   0.000000e+00     15s
#    14273    4.5536031e+03   0.000000e+00   0.000000e+00     20s

# Root relaxation: objective 4.553603e+03, 14273 iterations, 15.93 seconds (25.86 work units)

#     Nodes    |    Current Node    |     Objective Bounds      |     Work
#  Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

#      0     0 4553.60306    0   38 5126.69846 4553.60306  11.2%     -   20s
# H    0     0                    4713.2632254 4553.60306  3.39%     -   21s
# H    0     0                    4553.9966690 4553.60306  0.01%     -   21s
#      0     0 4553.67298    0   41 4553.99667 4553.67298  0.01%     -   31s
#      0     0 4553.69317    0   42 4553.99667 4553.69317  0.01%     -   34s
#      0     0 4553.70159    0   43 4553.99667 4553.70159  0.01%     -   36s
#      0     0 4553.70744    0   43 4553.99667 4553.70744  0.01%     -   38s
#      0     0 4553.70961    0   43 4553.99667 4553.70961  0.01%     -   39s
#      0     0 4553.71170    0   43 4553.99667 4553.71170  0.01%     -   40s
#      0     0 4553.71344    0   43 4553.99667 4553.71344  0.01%     -   41s
#      0     0 4553.71441    0   43 4553.99667 4553.71441  0.01%     -   42s
#      0     0 4553.71544    0   43 4553.99667 4553.71544  0.01%     -   43s
#      0     0 4553.71696    0   43 4553.99667 4553.71696  0.01%     -   44s
#      0     0 4553.71786    0   45 4553.99667 4553.71786  0.01%     -   45s
#      0     0 4553.71878    0   43 4553.99667 4553.71878  0.01%     -   46s
#      0     0 4553.71920    0   43 4553.99667 4553.71920  0.01%     -   46s
#      0     0 4553.72608    0   45 4553.99667 4553.72608  0.01%     -   73s
#      0     0 4553.72796    0   45 4553.99667 4553.72796  0.01%     -   75s
#      0     0 4553.72897    0   45 4553.99667 4553.72897  0.01%     -   76s
#      0     0 4553.72991    0   45 4553.99667 4553.72991  0.01%     -   77s
#      0     0 4553.73067    0   45 4553.99667 4553.73067  0.01%     -   77s
#      0     0 4553.73123    0   45 4553.99667 4553.73123  0.01%     -   78s
# ^C
# Cutting planes:
#   MIR: 1311

# Explored 1 nodes (21656 simplex iterations) in 86.54 seconds (92.30 work units)
# Thread count was 1 (of 256 available processors)

# Solution count 2: 4554 4713.26 

# Solve interrupted
# Best objective 4.553996669013e+03, best bound 4.553731234623e+03, gap 0.0058%

# User-callback calls 26229, time in user-callback 0.01 sec
# eval_MIP_Runtime = 86.53938102722168
# Sol=1: lb = 4553.594916892824 < 4563.413564040575 = ub, agap = 9.818647147751108, rgap = 0.002151601429491611
# Sol=2: lb = 4553.594916892824 < 4563.965374880923 = ub, agap = 10.370457988099588, rgap = 0.0022722472973121886
# final_eval: 2349.517627 seconds (361.65 k allocations: 19.233 MiB, 15 lock conflicts, 0.14% compilation time)
# ending memory = 10.9178466796875 GiB
