using LinearSolve

include("matrices.jl")

# solvestationary(n_H, T_l, W, T_s, z, v; n_levels = n_levels, max_it = 100, ε = 1e-6)

function findstatk(n_H, T_l, W, T_s, z, v, target; n_levels = 15, kwargs...)
    n_stat = solvestationary(n_H, T_l, W, T_s, z, v, n_levels = n_levels, kwargs...)
    lg_k_low = -20.0
    lg_k_high = 0.0
    function ionnonstat(k)
        δf_stat = stationarycomplexamp(n_stat, n_H, T_l, W, T_s, z, v; varpars = ("T_l",))
        δf_nonstat = nonstationarycomplexamp(n_stat, δf_stat, n_H, T_l, W, T_s, z, v, k)
        return abs(δf_stat[1])/abs(δf_nonstat[1]) - target
    end

    while lg_k_high - lg_k_low > 1e-1
        lg_k_mid = (lg_k_low + lg_k_high)/2
        f_high = ionnonstat(10^lg_k_high)
        f_low = ionnonstat(10^lg_k_low)
        if f_high*f_low > 0
            return NaN
        end
        f_mid = ionnonstat(10^lg_k_mid)
        if f_high*f_mid ≤ 0
            lg_k_low = lg_k_mid
        else
            lg_k_high = lg_k_mid
        end
    end

    return 10^lg_k_high
end

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

function nonstationarycomplexamp(n_stat, statcomplexamps, n_H, T_l, W, T_s, z, v, k; nonstat_levels = 0)
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