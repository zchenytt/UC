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
using CairoMakie
using LaTeXStrings

a = JLD2.load("algo3data1.jld2");
msttxs = a["msttxs"];
mstroiys = a["mstroiys"];
mstlbys = a["mstlbys"];
mstrotys = a["mstrotys"];
subtxs = a["subtxs"];
subrotys = a["subrotys"];
subgapys = a["subgapys"];
subvioys = a["subvioys"];
subVμys = a["subVμys"];
map(length, (msttxs, mstroiys, mstlbys, mstrotys, subtxs, subrotys, subgapys, subvioys, subVμys))
for (i,e)=enumerate(subgapys) e < 1e-6 && (subgapys[i] = 1e-6) end

f = Figure(size=(550, 1100));
ax1 = Axis(f[1, 1], 
    ylabel = "Master lb",
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
    yaxisposition = :right,
    xticksvisible=false, 
    xticklabelsvisible=false
);hidespines!(ax2);hidexdecorations!(ax2);
scatter!(ax1, msttxs, mstlbys, color = :snow4, markersize=1.5);
scatter!(ax2, msttxs, mstrotys, color = :green, markersize=1.5);
ax3 = Axis(f[2, 1];
    ylabel = "Master lb",
    ylabelcolor = :snow4,
    yticklabelcolor = :snow4,
    ylabelsize = 13,
    xticksvisible=false,
    xticklabelsvisible=false
);
ax4 = Axis(f[2, 1];
    ylabel = "Master Reoptimize #SimplexIteration",
    ylabelcolor = :coral,
    yticklabelcolor = :coral,
    ylabelsize = 13,
    yaxisposition = :right,
    xticksvisible=false, 
    xticklabelsvisible=false
);hidespines!(ax4);hidexdecorations!(ax4);
scatter!(ax3, msttxs, mstlbys, color = :snow4, markersize=1.5);
ylims!(ax3, 3600., 4015.)
scatter!(ax4, msttxs, mstroiys, color = :coral, markersize=1.5);
ax5 = Axis(f[3, 1];
    ylabel = "Subproblem Terminating Gap",
    ylabelsize = 13,
    yscale = log10,
    ylabelcolor = :goldenrod2,
    yticklabelcolor = :goldenrod2,
    xticksvisible=false, 
    xticklabelsvisible=false
);
ax6 = Axis(f[4, 1];
    ylabel = "Subproblem Reoptimize Time (s)",
    ylabelcolor = :blueviolet,    
    yticklabelcolor = :blueviolet,
    ylabelsize = 13,
    xticksvisible=false, 
    xticklabelsvisible=false
)
ax7 = Axis(f[5,1];
    ylabel = "Vio",
    ylabelcolor = :red,
    yticklabelcolor = :red,
    ylabelsize = 16,
    xticksvisible=false, 
    xticklabelsvisible=false
)
ax8 = Axis(f[5,1];
    ylabel = L"$\bar{v}$",
    ylabelcolor = :blue,
    yticklabelcolor = :blue,
    yaxisposition = :right,
    ylabelsize = 18,
    xticksvisible=false, 
    xticklabelsvisible=false
); hidespines!(ax8); hidexdecorations!(ax8);
ax9 = Axis(f[6, 1];
    ylabel = "Vio",
    ylabelcolor = :red,
    yticklabelcolor = :red,
    ylabelsize = 16,
    xlabel = "Runtime of Algorithm 3 (s) (x-axis is shared)",
);
ax10 = Axis(f[6,1];
    ylabel = L"$\bar{v}$",
    ylabelcolor = :blue,
    yticklabelcolor = :blue,
    yaxisposition = :right,
    ylabelsize = 18,
    xticksvisible=false, 
    xticklabelsvisible=false
); hidespines!(ax10); hidexdecorations!(ax10);
scatter!(ax5, subtxs, subgapys, color = :goldenrod2, markersize=0.9);
scatter!(ax6, subtxs, subrotys, color = :blueviolet, markersize=1.5);
scatter!(ax7, subtxs, subvioys, color = :red, markersize=1.2);
scatter!(ax8, subtxs, subVμys, color = :blue, markersize=0.7);
scatter!(ax9, subtxs, subvioys, color = :red, markersize=1.2);
scatter!(ax10, subtxs, subVμys, color = :blue, markersize=0.9);
ylims!(ax9, 0., 900.);
ylims!(ax10, 0., 900.);
f
CairoMakie.save("longplot.png", f)
