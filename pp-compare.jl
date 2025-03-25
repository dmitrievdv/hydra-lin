using Makie
using GLMakie

include("ampsolvers.jl")

z_arr = [5:0.1:13;]
n_z = length(z_arr)
T_l = 7e3
n_H = 1e12
T_s = 4e3
W = 0.01

v = 1e7

n_levels = 15

n_stat = zeros((n_levels, n_z))
b_stat = zeros((n_levels, n_z))

pp_n_stat = zeros((n_levels, n_z))
pp_b_stat = zeros((n_levels, n_z))

for i_z = 1:n_z
    z = 10^z_arr[i_z]
    println(z)
    
    n_stat[:, i_z] = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 100, ε = 1e-6,  plane_parallel = false)
    b_stat[:, i_z] = [n_stat[lev, i_z] / lteni(lev, n_H - sum(n_stat[:,i_z]), T_l) for lev in 1:n_levels]

    pp_n_stat[:, i_z] = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 100, ε = 1e-6,  plane_parallel = true)
    pp_b_stat[:, i_z] = [pp_n_stat[lev, i_z] / lteni(lev, n_H - sum(pp_n_stat[:,i_z]), T_l) for lev in 1:n_levels]
end

begin
    fig = Figure()
    ax = Axis(fig[1,1])

    plt_lev = 5

    for lev in 1:plt_lev
        lines!(ax, z_arr, log10.(b_stat[lev, :]))
    end
    # fig

    for lev in 1:plt_lev
        lines!(ax, z_arr, log10.(pp_b_stat[lev, :]), linestyle = :dash)
    end
    fig
end