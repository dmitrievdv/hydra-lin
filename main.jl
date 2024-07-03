using Plots
using LinearSolve

include("matrices.jl")

lg_z_arr = [3:0.1:14;]
n_z = length(lg_z_arr)
n_levels = 15
v = 3e7
W = 1e-2
T_l = 8e3
T_s = 5e3
n_H = 1e12

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

z = 1e11
k = 1e-10


T_ls = [1e3:1e1:20e3;]
n_Tl = length(T_ls)
ns = zeros(n_levels, n_Tl)
ΔAs = zeros(n_levels, n_Tl)
Δϕs = zeros(n_levels, n_Tl)

A0s = zeros(n_levels, n_Tl)
ϕ0s = zeros(n_levels, n_Tl)

Ls = zeros(n_levels, n_levels, n_Tl)
Vs = zeros(n_levels, n_Tl)
Δfs = zeros(Complex, n_levels, n_Tl)
f0s = zeros(Complex, n_levels, n_Tl)
fs = zeros(Complex, n_levels, n_Tl)

f1s = zeros(Complex, n_levels, n_Tl)
Δf1s = zeros(Complex, n_levels, n_Tl)
ΔA1s = zeros(n_levels, n_Tl)
Δϕ1s = zeros(n_levels, n_Tl)

f0s_res = zeros(Complex, n_levels, n_Tl)
fs_res = zeros(Complex, n_levels, n_Tl)
Δf_res = zeros(Complex, n_levels, n_Tl)

for j = 1:n_Tl
    T_l = T_ls[j]
    n = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 100, ε = 1e-6)
    α = calculatelinearizedescapepropability(v, z, n)
    β = calculateescapepropabilities(v, z, n)
    n_e = n_H - sum(n)

    # j = 3
    ϵ = 1e-4
    δnH = 1
    δTl = 1

    δσ_ij = linearisedpopulationsmatrix(n, α, β, n_e, T_l, T_s, W)
    δσ_inH = linearizednhsources(n, n_e, β, W, T_l, T_s)
    δσ_iTl = linearizedTesources(n, n_e, β, W, T_l, T_s)

    L = zeros(ComplexF64, (n_levels, n_levels))
    V = zeros(ComplexF64, n_levels)
    V .= (-δσ_inH*δnH - δσ_iTl*δTl) ./ n
    ones = fill(1.0 + 0.0im, n_levels)
    # V .= reverse(V)
    # println(V)

    for i = 1:n_levels
        L[i, :] = (δσ_ij[i,:] .- δσ_inH[i]) .* n / n[i] #/ V[i]
    end

    problem = LinearProblem(L, V)
    solution = solve(problem)
    δf = solution.u
    # δf = L \ V
    f0s[:, j] = δf
    f0s_res[:, j] = L*δf - V
    A_0 = abs.(δf)
    ϕ_0 = @. imag(log(δf / A_0 * (1 + 0.0im)))

    V_nonstat = zeros(Complex{Float64}, n_levels)

    L[1,1] += k*v*1im
    V_nonstat[1] = -k*v*1im*δf[1]

    problem = LinearProblem(L, V_nonstat)
    solution = solve(problem)
    Δf1 = solution.u
    δf1 = δf + Δf1
    A1 = abs.((δf1))
    ϕ1 = @. imag(log(δf1 / abs(δf1)))

    Δϕ1 = ϕ1 - ϕ_0
    for i = 1:n_levels
       if Δϕ1[i] > π
            Δϕ1[i] -= 2π
        elseif Δϕ1[i] < -π
            Δϕ1[i] += 2π
        end
    end

    ΔA1 = A_0 ./ A1

    ΔA1s[:, j] = ΔA1
    Δϕ1s[:, j] = Δϕ1
    Δf1s[:, j] = Δf1
    f1s[:, j] = δf1

    n_nostat_levels = n_levels
    for i = 2:n_nostat_levels
        L[i,i] += k*v*1im# / V[i]
        V_nonstat[i] = -k*v*1im*δf[i]
    end

    # V_nonstat = fill(-k*v*1im) .* δf# ./ V
    # V_nonstat .= reverse(V_nonstat)
    # Δf = L \ V_nonstat

    problem = LinearProblem(L, V_nonstat)
    solution = solve(problem)
    Δf = solution.u

    problem = LinearProblem(L, V)
    solution = solve(problem)
    # δf = solution.u
    fs_res[:, j] = L*δf - V
    # δf = L \ V
    δf = δf + Δf
    A = abs.((δf))
    ϕ = @. imag(log(δf / abs(δf)))

    Δϕ = ϕ - ϕ_0
    for i = 1:n_levels
       if Δϕ[i] > π
            Δϕ[i] -= 2π
        elseif Δϕ[i] < -π
            Δϕ[i] += 2π
        end
    end

    ΔA = A_0 ./ A

    


    ns[:, j] = n
    ΔAs[:, j] = ΔA
    Δϕs[:, j] = Δϕ
    A0s[:, j] = A_0
    ϕ0s[:, j] = ϕ_0
    Ls[:,:,j] .= abs.(L)
    Vs[:, j] .= abs.(V)
    Δfs[:, j] = Δf
    fs[:, j] = δf

end

begin
lev = 3
plt = plot(T_ls, log10.(Ls[lev,lev, :]), lc = :black) 
for i = 1:10
    plot!(plt, T_ls, log10.(Ls[lev,lev+i, :]), lc = i) 
    plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
end
plt
end

begin
    lev = 3
    plt = plot(T_ls, log10.(abs.(real.(fs[1,:])) ./ abs.(f0s[1,:])), lc = 1, label = "1, real") 
    plot!(plt, T_ls, log10.(abs.(imag.(fs[1,:])) ./ abs.(f0s[1,:])), lc = 1, ls = :dash, label = "1, imag") 
    for i = 2:lev
        plot!(plt, T_ls, log10.(abs.(real.(fs[i,:])) ./ abs.(f0s[i,:])), lc = i, label = "$i, real") 
        plot!(plt, T_ls, log10.(abs.(imag.(fs[i,:])) ./ abs.(f0s[i,:])), lc = i, ls = :dash, label = "$i, imag") 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(T_ls, log10.(A0s[1,:]), lc = 1) 
    plot!(plt, T_ls, log10.(A0s[1,:] ./ ΔAs[1,:]), lc = 1, ls = :dash) 
    for i = 2:lev
        plot!(plt, T_ls, log10.(A0s[i,:]), lc = i) 
        plot!(plt, T_ls, log10.(A0s[i,:] ./ ΔAs[i,:]), lc = i, ls = :dash) 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(T_ls, log10.(ΔAs[1,:]) , lc = 1) 
    plot!(plt, T_ls, log10.(ΔA1s[1,:]), lc = 1, ls = :dash) 
    for i = 2:lev
        plot!(plt, T_ls, log10.(ΔAs[i,:]), lc = i) 
        plot!(plt, T_ls, log10.(ΔA1s[i,:]), lc = i, ls = :dash) 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(T_ls, (Δϕs[1,:]), lc = 1) 
    for i = 2:lev
        plot!(plt, T_ls, (Δϕs[i,:]), lc = i) 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 1
    
    plt = plot(real.(Δfs[lev,:] ./ real.(f0s[lev,:])), imag.(Δfs[lev,:] ./ real.(f0s[lev,:])), line_z = T_ls)
    xlims!(plt, -2, 2)
    ylims!(plt, -2, 2)
end

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
# println((δσ - (δσ_ij[:,j] - δσ_inH))./σ1):log10.(abs.(δfs_stat[1,:,i_W])) - log10.(abs.(δfs_nonstat[1,:,i_W]))