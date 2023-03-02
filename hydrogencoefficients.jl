function Δik(i :: Number, k :: Number)
    z = - 4i*k/(i-k)^2
    _₂F₁(1-i, -k, 1, z)^2 - _₂F₁(1-k, -i, 1, z)^2
end

function gII(n, ν)
    E = planck_constant*ν
    Eₙ = rydberg_energy/n^2
    k = √(rydberg_energy/abs(E - Eₙ))
    if k > 1e10
        k = 1e10
    end
    Δ = Δik(n, 1im*k)
    return √3*π*k*n*exp(-4k*atan(n/k))/√(k^2+n^2)/(1-exp(-2π*k))*abs(Δ)
end

function limitfrequency(n :: Int)
    R = rydberg_energy/planck_constant
    return R/n^2
end


function boundfreeabsorbtioncoefficient(n :: Int, f) # f = ν/νₙ; νₙ = R/n²
    R = rydberg_energy/planck_constant
    ν = f*R/n^2
    # println(f)
    return 2.815e29*gII(n, ν)/n^5/ν^3
end

function gIIIkl(k, l)
    Δ = Δik(1im*k, 1im*l)
    return √3*π*k*l*exp(-2π*k)/(l-k)/(1-exp(-2π*k))/(1-exp(-2π*l))*abs(Δ)
end

function gIII(v, ν)
    Eₖ = electron_mass*v^2/2
    Eₗ = Eₖ + planck_constant*ν
    k = √(rydberg_energy/Eₖ)
    l = √(rydberg_energy/Eₗ)
    gIIIkl(k, l)
end

# super cool great approximation
function meangIIIaproximation(ν, T)
    E = planck_constant*ν
    return 1 + 0.1728*(E/rydberg_energy)^(1/3)*(1 + 2*boltzmann_constant*T/E)
end

function meangIII(ν, T)
    E = planck_constant*ν
    k = boltzmann_constant
    R = rydberg_energy
    f(x) = gIIIkl(√(R/(-log(x)*k*T)), √(R/(E-log(x)*k*T)))
    quadgk(f, 0, 1)
end

function freefreeansorbtioncoefficient(ν, T)
    E = planck_constant*ν
    k = boltzmann_constant
    return 3.69e8*meangIIIaproximation(ν, T)/ν^3/√(T)*(1 - exp(-E/(k*T)))
end

function stimulatedcorrection(ν, T)
    E = planck_constant*ν
    k = boltzmann_constant
    return (1 - exp(-E/(k*T)))
end

function hydrogenabsorbtion(populations, ne, ν, T; ε = 1e-3)
    R = rydberg_energy/planck_constant
    n_min = ceil(Int, √(R/ν))
    n_max = min(ceil(Int, n_min/ε^(1/5)), length(populations))
    bf_abs = 0.0
    for n = n_min:n_max
        # println(n)
        bf_abs += boundfreeabsorbtioncoefficient(n, ν/(R/n^2))*populations[n]
    end
    ff_abs = freefreeansorbtioncoefficient(ν, T)*ne^2
    return bf_abs + ff_abs
end

# function hydrogenemission(populations, ne, ν, T)