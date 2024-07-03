using Plots
using LinearSolve

include("matrices.jl")

# solvestationary(n_H, T_l, W, T_s, z, v; n_levels = n_levels, max_it = 100, ε = 1e-6)

function stationarycomplexamp(n_stat, n_H, T_l, W, T_s, z, v; varpars = ("T_l",))
    n_levels = length(n_stat)
    
    δnH = 0
    δTl = 0
    if "T_l" in varpars
        δTl = 1
    end; if "n_H" in varpars
        δnH = 1
    end

    n_e = n_H - sum(n_stat)
    α = calculatelinearizedescapepropability(v, z, n_stat)
    β = calculateescapepropabilities(v, z, n_stat)
    δσ_ij = linearisedpopulationsmatrix(n_stat, α, β, n_e, T_l, T_s, W)
    δσ_inH = linearizednhsources(n_stat, n_e, β, W, T_l, T_s)
    δσ_iTl = linearizedTesources(n_stat, n_e, β, W, T_l, T_s)



    L = zeros(ComplexF64, (n_levels, n_levels))
    V = zeros(ComplexF64, n_levels)
    V .= (-δσ_inH*δnH - δσ_iTl*δTl) ./ n_stat
    for i = 1:n_levels
        L[i, :] = (δσ_ij[i,:] .- δσ_inH[i]) .* n_stat / n_stat[i]
    end

    problem = LinearProblem(L, V)
    solution = solve(problem)
    return solution.u
end

function nonstationarycomplexamp(n_stat, statcomplexamps, n_H, T_l, W, T_s, z, v; nonstat_levels = 0)
    n_levels = length(n_stat)

    n_e = n_H - sum(n_stat)
    α = calculatelinearizedescapepropability(v, z, n_stat)
    β = calculateescapepropabilities(v, z, n_stat)
    δσ_ij = linearisedpopulationsmatrix(n_stat, α, β, n_e, T_l, T_s, W)
    δσ_inH = linearizednhsources(n_stat, n_e, β, W, T_l, T_s)

    L = zeros(ComplexF64, (n_levels, n_levels))
    V = zeros(ComplexF64, n_levels)
    # V .= (-δσ_inH*δnH - δσ_iTl*δTl) ./ n_stat
    for i = 1:n_levels
        L[i, :] = (δσ_ij[i,:] .- δσ_inH[i]) .* n_stat / n_stat[i]
    end

    n_nonstat_levels = if nonstat_levels == 0
        n_levels
    else 
        nonstat_levels
    end

    for i = 1:n_nonstat_levels
        L[i,i] += k*v*1im# / V[i]
        V[i] = -k*v*1im*statcomplexamps[i]
    end

    problem = LinearProblem(L, V)
    solution = solve(problem)
    return statcomplexamps + solution.u
end

T_ls = [7e3:2e1:8e3;]#[3e3:5e2:1e4;]
n_Tl = length(T_ls)
n_H = 1e11
n_levels = 15
v = 3e7
lgWs = [-2.75:0.02:-1.5;]; n_W = length(lgWs)
T_l = 8e3
T_s = 5.5e3
n_H = 1e12
z = 1e11
k = 1e-10

n_stats = zeros(n_levels, n_Tl, n_W)
b_stats = zeros(n_levels, n_Tl, n_W)
δfs_stat = zeros(Complex, (n_levels, n_Tl, n_W))
δfs_nonstat1 = zeros(Complex, (n_levels, n_Tl, n_W))
δfs_nonstat = zeros(Complex, (n_levels, n_Tl, n_W))

for i = 1:n_Tl, j = 1:n_W
    T_l = T_ls[i]; W = 10^lgWs[j]
    n_stat = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 100, ε = 1e-6)
    b_stat = [n_stat[lev] / lteni(lev, n_H - sum(n_stat), T_l) for lev in 1:n_levels]
    δf_stat = stationarycomplexamp(n_stat, n_H, T_l, W, T_s, z, v; varpars = ("T_l",))
    δf_nonstat1 = nonstationarycomplexamp(n_stat, δf_stat, n_H, T_l, W, T_s, z, v; nonstat_levels = 1)
    δf_nonstat = nonstationarycomplexamp(n_stat, δf_stat, n_H, T_l, W, T_s, z, v)

    n_stats[:, i, j] =  n_stat
    b_stats[:, i, j] =  b_stat
    δfs_stat[:, i, j] = δf_stat
    δfs_nonstat1[:, i, j] = δf_nonstat1
    δfs_nonstat[:, i, j] = δf_nonstat
end

W_samp = 5e-3
i_W = findmin(abs.(lgWs .- log10(W_samp)))[2]
Tl_samp = 7500
i_Tl = findmin(abs.(T_ls .- Tl_samp))[2]

begin
    lev = 5
    plt = plot(T_ls, log10.(abs.(δfs_stat[1,:,i_W])), lc = 1) 
    plot!(plt, T_ls, log10.(abs.(δfs_nonstat[1,:,i_W])), lc = 1, ls = :dash) 
    scatter!(plt, T_ls, log10.(abs.(δfs_nonstat1[1,:,i_W])), mc = 1, ms = 1, msw = 0) 
    for i = 2:lev
        plot!(plt, T_ls, log10.(abs.(δfs_stat[i,:,i_W])), lc = i) 
        plot!(plt, T_ls, log10.(abs.(δfs_nonstat[i,:,i_W])), lc = i, ls = :dash) 
        scatter!(plt, T_ls, log10.(abs.(δfs_nonstat1[i,:,i_W])), mc = i, ms = 1, msw = 0) 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(lgWs, log10.(abs.(δfs_stat[1,i_Tl,:])), lc = 1, label = "1, stat") 
    plot!(plt, lgWs, log10.(abs.(δfs_nonstat[1,i_Tl,:])), lc = 1, ls = :dash, label = "1, nonstat") 
    # scatter!(plt, lgWs, log10.(abs.(δfs_nonstat1[1,i_Tl,:])), mc = 1, ms = 1, msw = 0) 
    for i = 2:lev
        plot!(plt, lgWs, log10.(abs.(δfs_stat[i,i_Tl,:])), lc = i, label = "$i, stat") 
        plot!(plt, lgWs, log10.(abs.(δfs_nonstat[i,i_Tl,:])), lc = i, ls = :dash, label = "$i, nonstat") 
        # scatter!(plt, lgWs, log10.(abs.(δfs_nonstat1[i,i_Tl,:])), mc = i, ms = 1, msw = 0) 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(T_ls, log10.(abs.(δfs_stat[1,:,i_W])) - log10.(abs.(δfs_nonstat[1,:,i_W])), lc = 1, label = "1, stat/nonstat") 
    for i = 2:lev
        plot!(plt, T_ls, log10.(abs.(δfs_stat[i,:,i_W])) - log10.(abs.(δfs_nonstat[i,:,i_W])), lc = i, label = "$i, stat/nonstat") 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(lgWs, log10.(abs.(δfs_stat[1,i_Tl,:])) - log10.(abs.(δfs_nonstat[1,i_Tl,:])), lc = 1, label = "1, stat/nonstat") 
    for i = 2:lev
        plot!(plt, lgWs, log10.(abs.(δfs_stat[i,i_Tl,:])) - log10.(abs.(δfs_nonstat[i,i_Tl,:])), lc = i, label = "$i, stat/nonstat") 
        # plot!(plt, T_ls, log10.(Ls[lev+i,lev, :]), lc = i, ls = :dash) 
    end
    plt
end

begin
    lev = 5
    plt = plot(T_ls, log10.(b_stats[1,:,i_W]), lc = 1, label = "1")
    for i = 2:lev
        plot!(plt, T_ls, log10.(b_stats[i,:,i_W]), lc = i, label = "$i")
    endS
    plt
end

begin
    lev = 5
    plt = plot(lgWs, log10.(b_stats[1,i_Tl,:]), lc = 1, label = "1")
    for i = 2:lev
        plot!(plt, lgWs, log10.(b_stats[i,i_Tl,:]), lc = i, label = "$i")
    end
    plt
end

heatmap(lgWs, T_ls, log10.(abs.(δfs_stat[1,:,:])) - log10.(abs.(δfs_nonstat[1,:,:])))
