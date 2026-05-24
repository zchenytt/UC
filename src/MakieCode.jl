#######################################################

include("src/WindGen.jl")
using CairoMakie
fig = Figure(size = (1000, 400));
left = fig[1:4, 1] = GridLayout();
axs_left = [Axis(left[i, j]) for i = 1:4, j = 1:2];
axs_left[1, 1].title = "Single Scenario"
axs_left[1, 2].title = "Single Scenario"
right = fig[1:4, 2] = GridLayout();
ax_r1 = Axis(right[1:2, 1]);
ax_r1.title = "200 Scenarios, small variance"
ax_r2 = Axis(right[3:4, 1]);
ax_r2.title = "200 Scenarios, large variance"
colgap!(fig.layout, 1, 20)  # gap between left and right blocks
rowgap!(left, 10)
rowgap!(right, 10)
fig
include("src/WindGen.jl")
S = 200
m = WindGen._S(S, 1);
for i = 1:S
    lines!(ax_r2, 1:24, m[i][:, 1], color=:red)
    lines!(ax_r2, 1:24, m[i][:, 2], color=:blue)
    lines!(ax_r2, 1:24, m[i][:, 3], color=:green)
end
fig
S = 8;
m = WindGen._S(S, 1);
i = 0;
for ax = axs_left
    i += 1
    lines!(ax, 1:24, m[i][:, 1], color=:red)
    lines!(ax, 1:24, m[i][:, 2], color=:blue)
    lines!(ax, 1:24, m[i][:, 3], color=:green)
end
ax_r2.xlabel = "time (h)"
axs_left[4, 1].xlabel = "time (h)"
axs_left[4, 2].xlabel = "time (h)"

#######################################################
