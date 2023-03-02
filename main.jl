using Plots

include("matrices.jl")

lg_z_arr = [3:0.1:14;]
n_z = length(lg_z_arr)
n_levels = 15
v = 3e7
W = 1e-2
T_l = 11e3
T_s = 4e3
n_H = 1e11

n_matrix = zeros(n_z, n_levels)
b_matrix = zeros(n_z, n_levels)

levels = [1:n_levels;]

for i = 1:n_z
    z = 10^lg_z_arr[i]
    println(z)
    n = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 1)
    n_e = n_H - sum(n)
    n_matrix[i, :] .= n
    b_matrix[i, :] .= n ./ lteni.(levels, n_e, T_l)
end

b_matrix
plt = plot(lg_z_arr, log10.(n_matrix[:, 1]))
plot!(plt, lg_z_arr, log10.(n_matrix[:, 2]))
plot!(plt, lg_z_arr, log10.(n_matrix[:, 3]))
