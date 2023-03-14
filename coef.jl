using HypergeometricFunctions
using SpecialFunctions
using QuadGK

const qₑ = electron_charge = 4.803e-10
const h = planck_constant = 6.626e-27
const mₑ = electron_mass = 9.109e-28
const c = light_speed = 2.998e10
const mₚ = proton_mass = 1.673e-24
const bohr_radius = 5.292e-9
const kB = boltzmann_constant = 1.3806e-16
const Rₑ = rydberg_energy = 2.180e-11
const χ1 = rydberg_energy/boltzmann_constant
const menzel_cst = h^3/(2π*mₑ*kB)^(3/2)



include("johnson.jl")
include("hydrogencoefficients.jl")


function spontaneouseinstein(i :: Int, j :: Int)
    u,l = if i < j
        j,i
    elseif i > j
        i,j
    else
        return 0.0
    end
    νᵤₗ = linefrequency(u, l)
    fₗᵤ = oscillatorstrength(l, u)
    return 8π^2*qₑ^2*νᵤₗ^2/(mₑ*c^3)*fₗᵤ*l/u*j/i
end 

function einsteincoef(u :: Int, l :: Int)
    r = l/u
    d = (1-r^2)^(2/3)
    A_ul = 1 - 0.1728*(1+r^2)/(l^(2/3)*d) - 0.0496*(1-4/3*r^2+r^4)/(l^(4/3)*d^2)
    A_ul = 1.571e10*A_ul/((u^2 - l^2)*u^3*l)
end

"""
Transition function notation:

{b,f}{b,f}{lu,ul}{process}(...)

Examples:
- `bbcollision`     -- collisional bound-bound transition radiative
- `bbulspontaneous` -- spontaneous deexitation (ul means from upper (u) level to lower level (l)), 
the arguments (levels) must be in the specified order. Here: upper level first, lower level second. 
If there is no {ul,lu} specification the order is arbitrary.
- `bbluradiative`   -- radiative excitation
- `fbcollision`     -- collisional free-bound transition rate (three-particle recombination e + e + i -> e + n)
- `bfradiative`     -- radiative ionization
"""
transitions_notation = "For function naming notation see [`transitions_notation`](@ref)"


"""
    bbulspontaneous(u :: Int, l :: Int)

Compute spontaneous u -> l transion rate (u > l)  

$transitions_notation
"""
function bbulspontaneous(u :: Int, l :: Int)
    νᵤₗ = linefrequency(u, l)
    fₗᵤ = oscillatorstrength(l, u)
    return 8π^2*qₑ^2*νᵤₗ^2/(mₑ*c^3)*fₗᵤ*l^2/u^2
end

"""
    bbluspontaneous(l :: Int, u :: Int)

Compute spontaneous u -> l transion rate (u > l) multiplied by u^2/l^2 

NOT A REAL PROCESS! Simplifies equations  

$transitions_notation
"""
function bbluspontaneous(l :: Int, u :: Int)
    νᵤₗ = linefrequency(u, l)
    fₗᵤ = oscillatorstrength(l, u)
    return 8π^2*qₑ^2*νᵤₗ^2/(mₑ*c^3)*fₗᵤ
end

"""
    bbcollision(i :: Int, j :: Int, T)

Compute collisional transitions rate at temerature T for i -> j transition

$transitions_notation
"""
function bbcollision(i :: Int, j :: Int, T)
    if i < j
        bblucollision(i, j, T)
    elseif i > j
        bbulcollision(i, j, T)
    else
        0.0
    end
end

"""
    bblucollision(l :: Int, u :: Int, T)

Compute collisional transition rate at temperature T for l -> u transition (l < u)

$transitions_notation
"""
function bblucollision(l :: Int, u :: Int, T)
    collisiontransitionrate(l, u, T)
end


"""
    bblucollision(l :: Int, u :: Int, T)

Compute collisional transition rate at temperature T for l -> u transition (l < u)

$transitions_notation
"""
function bbulcollision(u :: Int, l :: Int, T)
    l^2/u^2*exp((χ1/l^2 - χ1/u^2)/T)*collisiontransitionrate(l, u, T)
end

"""
    bfcollision(i :: Int, T)

Compute collisional ionization rate at temperature T

$transitions_notation
"""
function bfcollision(i :: Int, T)
    collisionionizationrate(i :: Int, T)
end

"""
    fbcollision(i :: Int, T)

Compute collisional recombination rate (NOT spontaneous, three-particle) at temperature T

$transitions_notation
"""
function fbcollision(i :: Int, T)
    i^2*h/(2π*mₑ*kB*T)^(3/2)*exp(χ1/i^2/T)*collisionionizationrate(i, T)
end

"""
    fbspontaneous(i :: Int, T)

Compute spontaneous recombination rate at temperature T

$transitions_notation
"""
function fbspontaneous(i :: Int, T)
    spontaneousrecombinationrate(i, T)
end

"""
    linefrequencyul(u :: Int, l :: Int)

Compute transition u -> l frequency (u > l).
"""
function linefrequencyul(u :: Int, l :: Int)
    Rₑ/h*(1/l^2 - 1/u^2)    
end

"""
    linefrequencyul(u :: Int, l :: Int)

Compute transition i -> j frequency (i ≂̸ j).
"""
function linefrequency(i :: Int, j :: Int)
    Rₑ/h*abs(1/i^2 - 1/j^2)   
end

"""
    blackbodybracket(i :: Int, j:: Int, T_s)

Compute black-body radiation bracket at i -> j transition frequency. Defined from the Plank's law B_ν = 2hν³/c² * blackbody(...).\
T_s is the temperature.
"""
function blackbodybracket(i :: Int, j :: Int, T_s)
    1/(exp(abs(1/i^2 - 1/j^2)*χ1/T_s) - 1)
end

"""
    blackbody(i :: Int, j:: Int, T_s)

Compute black-body radiation intensity at i -> j transition frequency. T_s is the temperature.
"""
function blackbody(i :: Int, j :: Int, T_s)
    ν = linefrequency(i, j)
    blackbody(ν, T_s)
end

"""
    blackbody(ν :: Real, T_s :: Real)

Compute black-body radiation intensity at frequency ν. T_s is the temperature.
"""
function blackbody(ν :: Real, T_s :: Real)
    2*h*ν^3/c^2/(exp(h*ν/(kB*T_s)) - 1)
end

"""
    bflimitabsorbtioncoefficient(n :: Int)

Absorbtion coefficient for n -> c transition (ionization) at threshold frequency.
"""
function bflimitabsorbtioncoefficient(n)
    boundfreeabsorbtioncoefficient(n, 1.0)
end

"""
    bfradiative(i :: Int, T_s)

Compute radiative ionization rate. Radiation field is assumed to be black-body with temperature T_s.

$transitions_notation
"""
function bfradiative(i :: Int, T_s)
    ν_ci = limitfrequency(i)
    int_fun(z) = 1/z/(exp(χ1/i^2/T_s/z) - 1)
    int_res = quadgk(int_fun, 0, 1)[1]
    return 8π*bflimitabsorbtioncoefficient(i)*ν_ci^3/c^2*int_res
end


"""
    fbradiative(i :: Int, T_s, T_l)

Compute recombinations induced by radiation field with temperature T_s (T_source) in hydrogen gas with temperarure T_l (T_local)

$transitions_notation
"""
function fbradiative(i :: Int, T_s, T_l)
    ν_ci = limitfrequency(i)
    int_fun(z) = 1/z/(exp(χ1/i^2/T_s/z) - 1)*exp(-χ1/i^2/T_l/z)
    int_res = quadgk(int_fun, 0, 1)[1]
    return 8π*bflimitabsorbtioncoefficient(i)*ν_ci^3/c^2*int_res*i^2*h^3/(2π*mₑ*kB*T_l)^(3/2)*exp(χ1/i^2/T_l)
end