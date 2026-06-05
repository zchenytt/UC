include("src/WindGen.jl")
using CairoMakie
fig = Figure(size = (500, 600));
ax11 = Axis(fig[1, 1], xticks=0:2:23, xticklabelsvisible=false, xticksvisible=false, yticksvisible=false, title = "Single Scenario 1,2");
ax12 = Axis(fig[1, 2], xticks=0:2:23, xticklabelsvisible=false, xticksvisible=false, yticksvisible=false, title = "Single Scenario 3,4");
ax21 = Axis(fig[2, 1], xticks=0:2:23, xticksvisible=false, yticksvisible=false, title = "Scenario", titlevisible=false);
ax22 = Axis(fig[2, 2], xticks=0:2:23, xticksvisible=false, yticksvisible=false, title = "Scenario", titlevisible=false);
ax31 = Axis(fig[3, 1:2], xticks=0:1:23, xticklabelsvisible=false, xticksvisible=false, yticksvisible=false, title = "500 samples - Smaller Variance");
ax41 = Axis(fig[4, 1:2], xticks=0:1:23, xticksvisible=false, yticksvisible=false, title = "500 samples - Larger Variance");
S = 4
m = WindGen._S(S, 1);
lines!(ax11, 0:23, m[1][:, 1], color=:red)
lines!(ax11, 0:23, m[1][:, 2], color=:blue)
lines!(ax11, 0:23, m[1][:, 3], color=:green)
lines!(ax12, 0:23, m[2][:, 1], color=:red)
lines!(ax12, 0:23, m[2][:, 2], color=:blue)
lines!(ax12, 0:23, m[2][:, 3], color=:green)
lines!(ax21, 0:23, m[3][:, 1], color=:red)
lines!(ax21, 0:23, m[3][:, 2], color=:blue)
lines!(ax21, 0:23, m[3][:, 3], color=:green)
lines!(ax22, 0:23, m[4][:, 1], color=:red)
lines!(ax22, 0:23, m[4][:, 2], color=:blue)
lines!(ax22, 0:23, m[4][:, 3], color=:green)
S = 500
m = WindGen._S(S, 1);
for i=1:S
    lines!(ax41, 0:23, m[i][:, 1], color=:red, linewidth=0.1)
    lines!(ax41, 0:23, m[i][:, 2], color=:blue, linewidth=0.1)
    lines!(ax41, 0:23, m[i][:, 3], color=:green, linewidth=0.1)
end
fig



using JLD2
# @load "myvec.jld2" timensinalgo3mst lbv mstreopttime timensinalgo3 subtermigap
using CairoMakie
f = Figure(size=(680,400), figure_padding = 0.01);
ax1 = Axis(f[1, 1], ylabel = "Master Lower Bound", ylabelcolor=:blue, yticklabelcolor = :blue, xticksvisible=false, xticklabelsvisible=false);
ax2 = Axis(f[1, 1], ylabelcolor=:green, ylabel = "Master Reoptimize Time (s)", yticklabelcolor = :green, yaxisposition = :right);hidespines!(ax2);hidexdecorations!(ax2);
ax3 = Axis(f[2, 1]; yscale=log10, ylabel = "Subproblem Terminating Gap", xlabel="Runtime of Algorithm 3 (s)", xticksvisible=false);
lines!(ax1, timensinalgo3mst, lbv, color = :blue);
scatter!(ax2, timensinalgo3mst, mstreopttime, color = :green, markersize=4);
scatter!(ax3, timensinalgo3, subtermigap, markersize=3, color = :tomato);
f
