using LinearAlgebra

include("coef.jl")

"""
    ltene(n_H, T; levels = 5)

Compute electron concentration in hydrogen gas with temperature `T` and hydroen concentration `n_H`. Assuming `levels` level atom.
"""
function ltene(n_H, T; levels = 5)
    sum_neutral = 0.0
    for i ∈ 1:levels
        sum_neutral = sum_neutral + i^2*menzel_cst/T^1.5*exp(χ1/i^2/T)
    end
    return (-1 + √(1+4sum_neutral*n_H))/2sum_neutral
end


"""
    lteni(i :: Int, n_e, T)

Compute i-th level population in hydrogen gas with temperature `T` and electron concentration `n_e`. 
"""
function lteni(i :: Int, n_e, T)
    return i^2*menzel_cst/T^1.5*exp(χ1/i^2/T)*n_e^2
end

"""
    lteni(i :: Int, n_H, T; levels = 5)

Compute i-th level population in hydrogen gas with temperature `T` and hydrogen concentration `n_H`. 
Calls `n_e = ltene(n_H, T, levels = levels)` and then `lteni(i, n_e, T)`.
"""
function lteni(i :: Int, n_H, T, levels)
    n_e = ltene(n_H, T; levels = levels)
    lteni(i, n_e, T)
end


"""
    menzelne(n, n_H, T)

Compute `n_e` from menzel relations `b` if hydogen concentration is `n_H` and temperature is `T`.
"""
function menzelne(b, n_H, T)
    sum_neutral = 0.0
    levels = length(b)
    for i ∈ 1:levels
        sum_neutral = sum_neutral + i^2*menzel_cst/T^1.5*exp(χ1/i^2/T)*b[i]
    end
    return (-1 + √(1+4sum_neutral*n_H))/2sum_neutral
end

function calc_simple_ne(b1, n_H, T)
    sum_neutral = 0.0
    sum_neutral = sum_neutral + menzel_cst/T^1.5*exp(χ1/T)*b1
    return (-1 + √(1+4sum_neutral*n_H))/2sum_neutral
end

function calc_upper_levels_ne(n_1, b_upper, n_H, T)
    sum_neutral = 0.0
    levels = length(b_upper) + 1
    for i ∈ 2:levels
        sum_neutral = sum_neutral + i^2*menzel_cst/T^1.5*exp(χ1/i^2/T)*b_upper[i-1]
    end
    return (-1 + √(1+4sum_neutral*(n_H-n_1)))/2sum_neutral
end

"""
    exitpropability(τ)

Compute exitpropability `β = (1 - exp(-τ))/τ`.
"""
function exitpropability(τ)
    if abs(τ) < 1e-1
        return 1 - τ/2 + τ^2/6
    else
        return (1 - exp(-τ))/τ
    end
end

"""
    calculateescapepropabilities(v, z, n)

Calculate escape propabilities from 1d medium with given hydrogen level populations `n` using Sobolev's approximation. 
Returns a matrix β, where β[i,j] = β[j,i] are line escape propabilities and β[i,i] are continuum escappe propabilities.

The velocity gradient is `v/z`, where `v` is expansion velocity and `z` is length of the medium. 
The continnum `τ_c = bflimitabsorbtioncoefficient(i)*n[i]*z`.
"""
function calculateescapepropabilities(v, z, n)
    n_levels = length(n)
    β = zeros(n_levels, n_levels)
    for i = 1: n_levels
        τ = bflimitabsorbtioncoefficient(i)*n[i]*z
        β[i,i] = exitpropability(τ)
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            ν = linefrequencyul(i, j)
            τ = c^3/8π*A_ij/ν^3*i^2/j^2*abs(z/v)*n[j]*abs(1- n[i]/n[j]*j^2/i^2)
            β[i,j] = exitpropability(τ)
            β[j,i] = β[i,j]
        end
    end
    return β
end

function calculate_plane_parallel_escape_propability(v, z, n)
    n_levels = length(n)
    β = zeros(n_levels, n_levels)
    for i = 1: n_levels
        τ = bflimitabsorbtioncoefficient(i)*n[i]*z
        β[i,i] = exitpropability(τ)
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            ν = linefrequencyul(i, j)
            τ = c^3/8π*A_ij/ν^3*i^2/j^2*abs(z/v)*n[j]*abs(1- n[i]/n[j]*j^2/i^2)
            β[i,j] = plane_parallel_escape_propability(τ)
            β[j,i] = β[i,j]
        end
    end
    return β
end

function calculate_plane_parallel_linearized_escape_propability(v, z, n)
    n_levels = length(n)
    α = zeros(n_levels, n_levels)
    # E₂ = expint(2, τ)
    for i = 1: n_levels
        τ = bflimitabsorbtioncoefficient(i)*n[i]*z
        E₂ = expint(2, τ)
        α[i,i] = E₂ - plane_parallel_escape_propability(τ)
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            ν = linefrequencyul(i, j)
            τ = c^3/8π*A_ij/ν^3*i^2/j^2*abs(z/v)*n[j]*abs(1- n[i]/n[j]*j^2/i^2)
            # β = plane_parallel_escape_propability(τ)
            E₂ = expint(2, τ)
            α[i, j] = E₂ - plane_parallel_escape_propability(τ)
        end
    end
    return α
end

function plane_parallel_escape_propability(τ)
    return (1/2 - expint(3, τ))/τ
end


"""
    calculatelinearizedescapepropability(v, z, n)

Calculate `α = ∂β/∂τ*τ` (first taylor expansion term for `β`) where β is escape propability (see `calculateescapepropabilities`).
"""
function calculatelinearizedescapepropability(v, z, n)
    n_levels = length(n)
    α = zeros(n_levels, n_levels)
    for i = 1: n_levels
        τ = bflimitabsorbtioncoefficient(i)*n[i]*z
        β = exitpropability(τ)
        α[i,i] = 1 - β - β*τ
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            ν = linefrequencyul(i, j)
            τ = c^3/8π*A_ij/ν^3*i^2/j^2*abs(z/v)*n[j]*abs(1- n[i]/n[j]*j^2/i^2)
            β = exitpropability(τ)
            α[i,j] = 1 - β - β*τ
            α[j,i] = α[i,j]
        end
    end
    return α
end

"""
    populationsmatrix(β, n_e, T_l, T_s, W)

Computes matrix `M` and vector `R` from escape propabilities `β` and electron concentration `n_e`,
such that hydrogen level sources `σ = M*n + R`, where n is level populations.
`T_l`, `T_s`, `W` are local temperature, radiation source temperature and geometric dilution factor.
Returns a tuple `(M, R)`.
"""
function populationsmatrix(β, n_e, T_l, T_s, W)
    n_levels = size(β)[1]
    M = zeros(n_levels, n_levels)
    R = zeros(n_levels)
    for i = 1:n_levels
        A_ci = fbspontaneous(i, T_l)
        B_ci = fbradiative(i, T_s, T_l)
        B_ic = bfradiative(i, T_s)
        q_ic = bfcollision(i, T_l)
        β_ic = β[i,i]
        M[i,i] += -β_ic*W*B_ic - n_e*q_ic
        R[i] += β_ic*n_e^2*(A_ci + β_ic*W*B_ci) + n_e^3*i^2*h^3/(2π*mₑ*kB*T_l)^(3/2)*exp(χ1/i^2/T_l)*q_ic
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            J_ij = blackbodybracket(i, j, T_s)
            β_ji = β[j,i]
            q_ij = bbulcollision(i, j, T_l)
            M[i,i] += -A_ij*β_ji*(1 + W*J_ij) - n_e*q_ij
            M[i,j] += i^2/j^2*(A_ij*W*β_ji*J_ij + exp(χ1/T_l*(1/i^2 - 1/j^2))*q_ij*n_e)
        end
        for k = i+1:n_levels
            A_ki = bbulspontaneous(k, i)
            J_ki = blackbodybracket(k, i, T_s)
            q_ik = bbulcollision(k, i, T_l)
            β_ik = β[i,k]
            M[i,i] += -k^2/i^2*A_ki*W*β_ik*J_ki - n_e*q_ik
            M[i,k] += A_ki*β_ik*(1 + W*J_ki) + i^2/k^2*exp(χ1/T_l*(1/i^2 - 1/k^2))*q_ik*n_e
        end
    end
    return M, R
end

"""
    menzelmatrix(β, n_e, T_l, T_s, W)

Computes matrix `M` and vector `R` from escape propabilities `β` and electron concentration `n_e`,
such that hydrogen level sources `σ = M*(n ./ n_LTE) + R`, where n is level populations vector and 
n_LTE is LTE level populations vector. `T_l`, `T_s`, `W` are local temperature, 
radiation source temperature and geometric dilution factor.
Returns a tuple `(M, R)`.
"""
function menzelmatrix(β, n_e, T_l, T_s, W)
    n_levels = size(β)[1]
    M = zeros(n_levels, n_levels)
    R = zeros(n_levels)
    for i = 1:n_levels
        A_ci = fbspontaneous(i, T_l)
        B_ci = fbradiative(i, T_s, T_l)
        B_ic = bfradiative(i, T_s)
        q_ic = bfcollision(i, T_l)
        β_ic = β[i,i]
        M[i,i] += -β_ic*W*B_ic - n_e*q_ic
        R[i] += (2π*mₑ*kB*T_l)^(3/2)/(i^2*h^3)*exp(-χ1/i^2/T_l)*β_ic*(A_ci + β_ic*W*B_ci) + n_e*q_ic
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            J_ij = blackbodybracket(i, j, T_s)
            β_ji = β[j,i]
            q_ij = bbulcollision(i, j, T_l)
            M[i,i] += -A_ij*β_ji*(1 + W*J_ij) - n_e*q_ij
            M[i,j] += exp(χ1/T_l*(1/j^2 - 1/i^2))*A_ij*W*β_ji*J_ij + q_ij*n_e
        end
        for k = i+1:n_levels
            A_ki = bbulspontaneous(k, i)
            J_ki = blackbodybracket(k, i, T_s)
            q_ik = bblucollision(i, k, T_l)
            β_ik = β[i,k]
            M[i,i] += -k^2/i^2*A_ki*W*β_ik*J_ki - n_e*q_ik
            M[i,k] += k^2/i^2*exp(χ1/T_l*(1/k^2-1/i^2))*A_ki*β_ik*(1 + W*J_ki) + q_ik*n_e
        end
    end
    return M, R
end

function linearisedpopulationsmatrix(n, α, β, n_e, T_l, T_s, W)
    n_levels = size(β)[1]
    M = zeros(n_levels, n_levels)
    for i = 1:n_levels
        A_ci = fbspontaneous(i, T_l)
        B_ci = fbradiative(i, T_s, T_l)
        B_ic = bfradiative(i, T_s)
        q_ic = bfcollision(i, T_l)
        β_ic = β[i,i]
        α_ic = α[i,i]
        X_i = χ1/i^2/T_l
        M[i,i] += -β_ic*W*B_ic - n_e*q_ic + α_ic*(n_e^2/n[i]*(A_ci + 2β_ic*W*B_ci) - W*B_ic)
        for j = 1:i-1
            A_ij = bbulspontaneous(i, j)
            J_ij = blackbodybracket(i, j, T_s)
            S_ij = 1/abs(n[j]/n[i]*i^2/j^2 - 1)
            β_ji = β[j,i]
            α_ji = α[j,i]
            q_ij = bbulcollision(i, j, T_l)
            X_j = χ1/j^2/T_l
            M[i,i] += A_ij*(α_ji*(S_ij - W*J_ij) - W*β_ji*J_ij) - A_ij*β_ji - n_e*q_ij
            M[i,j] = A_ij*i^2/j^2*W*J_ij*(α_ji + β_ji) - A_ij*n[i]*i^2*α_ji/(n[j]*i^2-n[i]*j^2) + n_e*i^2/j^2*exp(X_i - X_j)*q_ij
            # println("$i $j $(A_ij*(α_ji*(S_ij - W*J_ij) - β_ji)) $(-A_ij*β_ji) $(n_e*q_ij)")
        end
        for k = i+1:n_levels
            A_ik = bbluspontaneous(i, k)
            J_ik = blackbodybracket(k, i, T_s)
            q_ik = bblucollision(i, k, T_l)
            S_ik = 1/abs(n[i]/n[k]*k^2/i^2 - 1)
            β_ki = β[i,k]
            α_ki = α[i,k]
            X_k = χ1/k^2/T_l
            M[i,i] += A_ik*(α_ki*(S_ik - W*J_ik) - W*β_ki*J_ik) - n_e*q_ik
            M[i,k] = A_ik*i^2/k^2*(1 + W*J_ik)*(α_ki + β_ki) - A_ik*n[i]*i^2*α_ki/(n[i]*k^2-n[k]*i^2) + n_e*i^2/k^2*exp(X_i - X_k)*q_ik
            # println("$i $k $(A_ik*(α_ki*(S_ik - W*J_ik) - β_ki)) $(n_e*q_ik)")
        end
        # println("$i $(M[i,i]) $(-β_ic*W*B_ic) $(-n_e*q_ic) $(α_ic*(n_e^2/n[i]*(A_ci + 2β_ic*W*B_ci) - W*B_ic))")
    end
    return M
end

function linearizednhsources(n, n_e, β, W, T_l, T_s)
    n_levels = length(n)
    δσ = zeros(n_levels)
    for i = 1:n_levels
        A_ci = fbspontaneous(i, T_l)
        B_ci = fbradiative(i, T_s, T_l)
        β_ic = β[i,i]
        X_i = χ1/T_l/i^2
        C_i = i^2*h^3/(2π*mₑ*kB*T_l)^(3/2)*exp(X_i)
        q_ic = bfcollision(i, T_l)
        δσ[i] += 2*n_e*β_ic*(A_ci + β_ic*W*B_ci)
        δσ[i] += q_ic*(3*n_e^2*C_i - n[i])
        for j = 1:i-1
            X_j = χ1/T_l/j^2
            q_ij = bbulcollision(i, j, T_l)
            δσ[i] += q_ij*(n[j]*i^2/j^2*exp(X_i - X_j) - n[i])
        end
        for k = i+1:n_levels
            X_k = χ1/T_l/k^2
            q_ik = bblucollision(i, k, T_l)
            δσ[i] += q_ik*(n[k]*i^2/k^2*exp(X_i - X_k) - n[i])
        end
    end
    return δσ
end

function linearizedTesources(n, n_e, β, W, T_l, T_s; ϵ = 1e-8)
    n_levels = length(n)
    δσ = zeros(n_levels)
    dT_l = ϵ*T_l
    T_l2 = T_l*(1+ϵ) 
    for i = 1:n_levels
        q_ic = bfcollision(i, T_l)
        A_ci = fbspontaneous(i, T_l)
        B_ci = fbradiative(i, T_s, T_l)
        β_ic = β[i,i]
        ∂q_ic = (bfcollision(i, T_l2) - q_ic)/dT_l
        ∂A_ic = (fbspontaneous(i, T_l2) - A_ci)/dT_l
        ∂B_ic = (fbradiative(i, T_s, T_l2) - B_ci)/dT_l
        X_i = χ1/T_l/i^2
        C_i = i^2*h^3/(2π*mₑ*kB*T_l)^(3/2)*exp(X_i)
        δσ[i] += n_e^3*C_i*∂q_ic
        δσ[i] += -n_e^3*C_i*(3/2 + X_i)/T_l*q_ic
        δσ[i] += β_ic*n_e^2*(∂A_ic + β_ic*W*∂B_ic)
        δσ[i] += -n_e*n[i]*∂q_ic
        for j = 1:i-1
            X_j = χ1/T_l/j^2
            q_ij = bbulcollision(i, j, T_l)
            ∂q_ij = (bbulcollision(i, j, T_l2) - q_ij)/dT_l
            δσ[i] += n_e*n[j]*i^2/j^2*exp(X_i - X_j)*(∂q_ij - (X_i - X_j)*q_ij/T_l)
            δσ[i] += -n_e*n[i]*∂q_ij
        end
        for k = i+1:n_levels
            X_k = χ1/T_l/k^2
            q_ik = bblucollision(i, k, T_l)
            ∂q_ik = (bblucollision(i, k, T_l2) - q_ik)/dT_l
            δσ[i] += n_e*n[k]*i^2/k^2*exp(X_i - X_k)*(∂q_ik - (X_i - X_k)*q_ik/T_l)
            δσ[i] += -n_e*n[i]*∂q_ik
        end
    end
    return δσ
end

"""
    solvestationary(n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100)

Solve stationary equation `σ_i = 0` in  expanding (v, z) 1d medium (n_H, T_l) with one radiation source (W, T_s)/

"""
function solvestationary(n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100, plane_parallel = false, simple_ne = false)
    n = zeros(n_levels)
    b = fill(1.0, n_levels)
    # n_new
    n_e = ltene(n_H, T_l)
    for i = 1:n_levels
        n[i] = lteni(i, n_e, T_l)
    end
    abs_Δb = 1e10 
    it = 0
    while (abs_Δb > ε) & (it < max_it)
        it += 1
        β = if plane_parallel 
            calculate_plane_parallel_escape_propability(v, z, n)
        else
            calculateescapepropabilities(v, z, n)
        end
        M, R = menzelmatrix(β, n_e, T_l, T_s, W)
        b_new = M\(-R)
        abs_Δb = sqrt(sum(@. (b_new/b - 1)^2))
        b .= b_new
        n_e = if simple_ne
            calc_simple_ne(b[1], n_H, T_l) 
        else
            menzelne(b, n_H, T_l)
        end
        for i = 1:n_levels
            n[i] = lteni(i, n_e, T_l)*b[i]
        end
    end
    return n
end

function solve_stationary_upper_levels(n_1, n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100, plane_parallel = false, simple_ne = false)
    n_upper = zeros(n_levels-1)
    b_upper = fill(1.0, n_levels-1)
    n_e = if simple_ne
        n_H - n_1
    else
        calc_upper_levels_ne(n_1, b_upper, n_H, T_l)
    end
    # n_new
    # n_e = ltene(n_H, T_l)
    for i = 2:n_levels
        n_upper[i-1] = lteni(i, n_e, T_l)
    end
    b_1 = n_1 / lteni(1, n_e, T_l)
    abs_Δb = 1e10 
    it = 0

    n = zeros(n_levels)
    n[1] = n_1

    while (abs_Δb > ε) & (it < max_it)
        n[2:n_levels] .= n_upper
        it += 1
        b_1 = n_1 / lteni(1, n_e, T_l)
        β = if plane_parallel 
            calculate_plane_parallel_escape_propability(v, z, n)
        else
            calculateescapepropabilities(v, z, n)
        end
        M, R = menzelmatrix(β, n_e, T_l, T_s, W)
        b_upper_new = M[2:end,2:end]\(-R[2:end] - M[2:end,1]*b_1)
        abs_Δb = sqrt(sum(@. (b_upper_new/b_upper - 1)^2))
        b_upper .= b_upper_new
        n_e = if simple_ne
            n_H - n_1
        else
            calc_upper_levels_ne(n_1, b_upper, n_H, T_l)
        end
        for i = 2:n_levels
            n_upper[i-1] = lteni(i, n_e, T_l)*b_upper[i-1]
        end
    end
    return [n_1; n_upper]
end

function calc_first_level_derivative(n_1, n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100, plane_parallel = false, simple_ne = false)
    n = solve_stationary_upper_levels(n_1, n_H, T_l, W, T_s, z, v; n_levels = n_levels, ε = ε, max_it = max_it, 
                                                                   plane_parallel = plane_parallel, simple_ne = simple_ne)
    
    β = if plane_parallel 
        calculate_plane_parallel_escape_propability(v, z, n)
    else
        calculateescapepropabilities(v, z, n)
    end
    n_e = if simple_ne
        n_H - n_1
    else
        n_H - sum(n)
    end
    M, R = menzelmatrix(β, n_e, T_l, T_s, W)
    b = [n[i]/lteni(i, n_e, T_l) for i = 1:n_levels]
    σ_1 = (M[1,:] ⋅ b + R[1])
    return σ_1 * h^3/(2π*mₑ*kB*T_l)^(3/2)*n_e^2*exp(χ1/T_l)
end

function calc_logarithmic_derivative(n_1, n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100, plane_parallel = false, simple_ne = false)
    σ_1 = calc_first_level_derivative(n_1, n_H, T_l, W, T_s, z, v; n_levels = n_levels, ε = ε, max_it = max_it, 
                                                                   plane_parallel = plane_parallel, simple_ne = simple_ne)
    df_ds = σ_1/n_H/v
    f = n_1/n_H
    l_f = log(f/(1-f))
    return 1/f*1/(1-f)*df_ds
end
