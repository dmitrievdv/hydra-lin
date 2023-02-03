function g0(n :: Int)
    if n == 1
        return 1.1330
    elseif n == 2
        return 1.0785
    else
        return 0.9935 + 0.2328/n - 0.1296/n^2
    end
end

function g1(n :: Int)
    if n == 1
        return -0.4059
    elseif n == 2
        return -0.2319
    else
        return -(0.6282 - 0.5598/n + 0.5299/n^2)/n
    end
end

function g2(n :: Int)
    if n == 1
        return 0.07014
    elseif n == 2
        return 0.02947
    else
        return (0.3887 - 1.181/n + 1.470/n^2)/n^2
    end
end

function gaunt(n :: Int, x)
    return g0(n) + g1(n)/x + g2(n)/x^2
end

function oscillatorstrength(l :: Int, u :: Int, )
    x = 1 - l^2/u^2
    return 32/(3*√3*π)*l/u^3/x^3*gaunt(l, x)
end

function spontaneousrecombinationrate(n :: Int, T)
    D = 5.197e-14
    χn_T = χ1/n^2/T
    return D*(χn_T)^(3/2)*exp(χn_T)*(g0(n)*expint(1, χn_T) + g1(n)*expint(2, χn_T) + g2(n)*expint(3, χn_T))
end

function Alu(l :: Int, u :: Int)
    x = 1 - l^2/u^2
    return 2*l^2/x*oscillatorstrength(l, u)
end

function An(n :: Int)
    return 32/(3*√3*π)*n*(g0(n)/3 + g1(n)/4 + g2(n)/5)
end

function Blu(l :: Int, u :: Int)
    bn = if l == 1
        -0.603
    else
        (4 - 18.63/l + 36.24/l^2 - 28.09/l^3)/l
    end
    x = 1 - l^2/u^2
    return 4*l^4/u^3/x^2*(1 + 4/3/x + bn/x^2)
end

function Bn(n :: Int)
    bn = if n == 1
        -0.603
    else
        (4 - 18.63/n + 36.24/n^2 - 28.09/n^3)/n
    end
    return 2/3*n^2*(5 + bn)
end

function collisiontransitionrate(l :: Int, u :: Int, T)
    x = 1 - l^2/u^2
    constant = √(8*boltzmann_constant/π/electron_mass)*2π*bohr_radius^2
    y = (χ1/l^2 - χ1/u^2)/T
    rn = if l == 1
        0.45
    else
        1.94/l^(1.57)
    end
    z = rn*x + y
    return constant*√T*l^2*y^2/x*(Alu(l, u)*((1/y + 0.5)*expint(1,y) - (1/z + 0.5)*expint(1,z)) + 
                        (Blu(l, u) - Alu(l,u)*log(2l^2/x))*(expint(2,y)/y - expint(2,z)/z))
end

function expintξ(t)
    return expint(0, t) - 2expint(1, t) + expint(2, t)
end

function collisionionizationrate(n :: Int, T)
    constant = √(8*boltzmann_constant/π/electron_mass)*2π*bohr_radius^2
    y = χ1/n^2/T
    rn = if n == 1
        0.45
    else
        1.94/n^(1.57)
    end
    z = rn + y
    # println(z)
    # println(y)
    return constant*√T*n^2*y^2*(An(n)*(expint(1,y)/y - expint(1,z)/z) + 
                        (Bn(n) - An(n)*log(2n^2))*(expintξ(y) - expintξ(z)))
end

function johnsoncollisions(levels :: Int, T_local)
    q = zeros(levels, levels)
    qc = zeros(levels)
    qs = zeros(levels)
    RR = zeros(levels)
    exp_χ = initexponents(T_local, levels)
    for l = 1:levels
        qc[l] = collisionionizationrate(l, T_local)
        RR[l] = spontaneousrecombinationrate(l, T_local)/l^2/menzel_cst*T_local^(3/2)/exp(χ1/T_local/l^2)
        for u = l+1:levels
            q[l,u] = collisiontransitionrate(l, u, T_local)
            q[u,l] = q[l,u]*l^2/u^2*exp_χ[l]/exp_χ[u]
        end
    end
    return q, qc, qs, RR
end