using Plots
# using Distributed
# n_procs = 4
# addprocs(n_procs)

# @everywhere using SharedArrays




include("ampsolvers.jl")

# solvestationary(n_H, T_l, W, T_s, z, v; n_levels = n_levels, max_it = 100, ε = 1e-6)

begin 
    T_ls = [1e3:5e2:15e3;]
    n_Tl = length(T_ls)
    n_H = 1e11
    n_levels = 15
    v = 1e7
    # lgWs = [-3:0.1:0;]; n_W = length(lgWs)  
    W = 1e-2
    T_l = 8e3
    T_s = 4e3
    lg_nHs = [9:0.1:13;]
    n_nH = length(lg_nHs)
    n_H = 1e11
    z = 1e11
    k = 1e-11
    n_levels = 15

    ks = zeros(n_Tl, n_nH)

    for i=1:n_Tl, j=1:n_nH
        n_H = 10^lg_nHs[j]
        T_l = T_ls[i]
        ks[i,j] = findstatk(n_H, T_l, W, T_s, z, v, 2)

    end
    heatmap(lg_nHs, T_ls, log10.(ks), clims = (-12,-7))
end



begin 
T_ls = [1e3:5e2:15e3;]
n_Tl = length(T_ls)
n_H = 1e11
n_levels = 15
v = 1e7
lgWs = [-3:0.1:0;]; n_W = length(lgWs)
T_l = 8e3
T_s = 4e3
n_H = 1e11
z = 1e11
k = 1e-11
n_levels = 15

n_stats = zeros((n_levels, n_Tl, n_W))
b_stats = zeros((n_levels, n_Tl, n_W))
δfs_stat = zeros(ComplexF64, (n_levels, n_Tl, n_W))
δfs_nonstat1 = zeros(ComplexF64, (n_levels, n_Tl, n_W))
δfs_nonstat = zeros(ComplexF64, (n_levels, n_Tl, n_W))

for i = 1:n_Tl; for j = 1:n_W
    T_l = T_ls[i]; W = 10^lgWs[j]
    n_stat = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, max_it = 100, ε = 1e-6)
    b_stat = [n_stat[lev] / lteni(lev, n_H - sum(n_stat), T_l) for lev in 1:n_levels]
    δf_stat = stationarycomplexamp(n_stat, n_H, T_l, W, T_s, z, v; varpars = ("T_l",))
    δf_nonstat1 = nonstationarycomplexamp(n_stat, δf_stat, n_H, T_l, W, T_s, z, v, k; nonstat_levels = 1)
    δf_nonstat = nonstationarycomplexamp(n_stat, δf_stat, n_H, T_l, W, T_s, z, v, k)

    n_stats[:, i, j] .=  n_stat
    b_stats[:, i, j] .=  b_stat
    δfs_stat[:, i, j] .= δf_stat
    δfs_nonstat1[:, i, j] .= δf_nonstat1
    δfs_nonstat[:, i, j] .= δf_nonstat
end; end

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
    end
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
end