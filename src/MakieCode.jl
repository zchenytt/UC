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


import JLD2
a = JLD2.load("algo3data0.jld2");

subtxs = a["subtxs"];
subgapys = a["subgapys"];
msttxs = a["msttxs"]; #
mstlbys = a["mstlbys"]; #
mstrotys = a["mstrotys"]; #
mstroiys = a["mstroiys"]; #
mstncys = a["mstncys"];
map(length, (subtxs, subgapys, msttxs, mstlbys, mstrotys, mstroiys, mstncys));

using CairoMakie
f = Figure(size=(500,500), figure_padding = 1.);
ax1 = Axis(f[1, 1], 
    ylabel = "Master Lb",
    ylabelcolor=:snow4, 
    yticklabelcolor = :snow4, 
    xticksvisible=false, 
    xticklabelsvisible=false
);
ax2 = Axis(f[1, 1], 
    ylabel = "Master Reoptimize Time (s)",
    ylabelcolor=:green, 
    yticklabelcolor = :green,
    ylabelsize = 11,
    yaxisposition = :right
);hidespines!(ax2);hidexdecorations!(ax2);
ax3 = Axis(f[2, 1];
    ylabel = "Master Reoptimize #Iter",
    ylabelcolor = :coral,
    yticklabelcolor = :coral,
    ylabelsize = 13,
    xticksvisible=false,
    yaxisposition = :right,
    xticklabelsvisible=false
);hidespines!(ax3);hidexdecorations!(ax3);
ax4 = Axis(f[2, 1];
    ylabel = "#Separating Cuts",
    ylabelcolor=:blue, 
    yticklabelcolor = :blue,
    xticksvisible=false, 
    xticklabelsvisible=false
);
ax5 = Axis(f[3, 1];
    ylabel = "Subproblem Terminating Gap",
    ylabelsize = 13,
    yscale = log10,
    xticksvisible = false,
    xlabel = "Runtime of Algorithm 3 (s) (x-axis is shared)",
);
lines!(ax1, msttxs, mstlbys, color = :snow4);
scatter!(ax2, msttxs, mstrotys, color = :green, markersize=1.5);
scatter!(ax3, msttxs, mstroiys, color = :coral, markersize=1.5);
lines!(ax4, msttxs, mstncys, color = :blue);
scatter!(ax5, subtxs, subgapys, color = :goldenrod2, markersize=0.9);
f
save("algo3plot3.pdf", f)
