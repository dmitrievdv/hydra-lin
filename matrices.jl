include("coef.jl")

"""
    ltene(n_H, T; levels = 5)

Compute electron concentration in hydrogen gas with temperature `T` and hydroen concentration `n_H`. Assuming `levels` level atom.
"""
function ltene(n_H, T; levels = 5)
    sum_neutral = 0.0
    # println(lte_deviations)
    for i ∈ 1:levels
        sum_neutral = sum_neutral + i^2*menzel_cst/T^1.5*exp(χ1/i^2/T)
        # println("$i ", lte_deviations[i])
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
        println(n[i])
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
            M[i,j] = i^2/j^2*(A_ij*W*β_ji*J_ij + exp(χ1/T_l*(1/i^2 - 1/j^2))*q_ij*n_e)
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
            M[i,j] = exp(χ1/T_l*(1/j^2 - 1/i^2))*A_ij*W*β_ji*J_ij + q_ij*n_e
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


"""
    solvestationary(n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100)

Solve stationary equation `σ_i = 0` in  expanding (v, z) 1d medium (n_H, T_l) with one radiation source (W, T_s)/

"""
function solvestationary(n_H, T_l, W, T_s, z, v; n_levels = 15, ε = 1e-3, max_it = 100)
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
        β = calculateescapepropabilities(v, z, n)
        M, R = menzelmatrix(β, n_e, T_l, T_s, W)
        b_new = M\(-R)
        abs_Δb = sqrt(sum(@. (b_new/b - 1)^2))
        b .= b_new
        n_e = menzelne(b, n_H, T_l)
        for i = 1:n_levels
            n[i] = lteni(i, n_e, T_l)*b[i]
        end
        # println("$n_e $abs_Δb")
        # println(b[1:5])
        for i = 1:5
            println(β[i, 1:5])
        end
    end
    return n
end