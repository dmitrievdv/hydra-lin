using Plots

include("matrices.jl")

lg_z_arr = [3:0.1:14;]
n_z = length(lg_z_arr)
n_levels = 30
v = 3e7
W = 1e-2
T_l = 8e3
T_s = 5e3
n_H = 1e11

n_matrix = zeros(n_z, n_levels)
b_matrix = zeros(n_z, n_levels)

levels = [1:n_levels;]

# for i = 1:n_z
#     z = 10^lg_z_arr[i]
#     println(z)
#     n = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 1)
#     n_e = n_H - sum(n)
#     n_matrix[i, :] .= n
#     b_matrix[i, :] .= n ./ lteni.(levels, n_e, T_l)
# end

# b_matrix
# plt = plot(lg_z_arr, log10.(n_matrix[:, 1]))
# plot!(plt, lg_z_arr, log10.(n_matrix[:, 2]))
# plot!(plt, lg_z_arr, log10.(n_matrix[:, 3]))

z = 1e12
k = 1e-10

function getampandphaseshift(n_H, T_l, W, T_s, z, v, k)
    n = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 100, ε = 1e-6)
    α = calculatelinearizedescapepropability(v, z, n)
    β = calculateescapepropabilities(v, z, n)
    n_e = n_H - sum(n)

    j = 3
    ϵ = 1e-1
    δnH = 1
    δTl = 1

    δσ_ij = linearisedpopulationsmatrix(n, α, β, n_e, T_l, T_s, W)
    δσ_inH = linearizednhsources(n, n_e, β, W, T_l, T_s)
    δσ_iTl = linearizedTesources(n, n_e, β, W, T_l, T_s)

    L = zeros(Complex, (n_levels, n_levels))
    V = zeros(Complex, n_levels)
    L .= δσ_ij
    V = (-δσ_inH*δnH - δσ_iTl*δTl) ./ n

    for i = 1:n_levels
        L[i, :] = (L[i, :] - δσ_inH) .* n / n[i]
    end

    δf = L\V
    A_0 = abs.(δf)
    ϕ_0 = @. imag(log(δf / A_0))

    n_nostat_levels = n_levels
    for i = 1:n_nostat_levels
        L[i,i] += k*v*1im
    end

    δf = L\V
    A = abs.(δf)
    ϕ = @. imag(log(δf / A))

    Δϕ = ϕ - ϕ_0
    for i = 1:n_levels
       if Δϕ[i] > π
            Δϕ[i] -= 2π
        elseif Δϕ[i] < -π
            Δϕ[i] += 2π
        end
    end

    return n, A_0 ./ A, Δϕ

    # 
end

T_ls = [7e3:1e1:12e3;]
n_Tl = length(T_ls)
ns = zeros(n_levels, n_Tl)
ΔAs = zeros(n_levels, n_Tl)
Δϕs = zeros(n_levels, n_Tl)

for j = 1:n_Tl
    T_l = T_ls[j]
    n, ΔA, Δϕ = getampandphaseshift(n_H, T_l, W, T_s, z, v, k)
    ns[:, j] = n
    ΔAs[:, j] = ΔA
    Δϕs[:, j] = Δϕ
end

plt = plot(T_ls, log10.(ΔAs[1, :]))
for i = 2:5
    plot!(plt, T_ls, log10.(ΔAs[i, :]))
end
plt
# n_e = n_H - sum(n)
# M, R = populationsmatrix(β, n_e, T_l, T_s, W)
# σ1 = M*n + R

# n_e = n_H + ϵ*n_H - sum(n)
# # β = calculateescapepropabilities(v, z, n)
# M, R = populationsmatrix(β, n_e, T_l, T_s, W)
# σ2 = M*n + R
# δσ = (σ2 - σ1)/(ϵ*n_H)
# println((δσ - δσ_inH)./σ1)

# n2 = deepcopy(n)
# n2[j] += ϵ*n[j]
# n_e2 = n_H - sum(n2)
# β2 = calculateescapepropabilities(v, z, n)
# M, R = populationsmatrix(β2, n_e2, T_l, T_s, W)
# σ2 = M*n2 + R

# δσ = (σ2 - σ1)/(ϵ*n[j])
# println((δσ - (δσ_ij[:,j] - δσ_inH))./σ1)