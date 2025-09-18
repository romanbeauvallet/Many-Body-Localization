#!/usr/bin/env julia

##################################### Libraries #############

using ITensors
using ITensorMPS
using MBL
using LinearAlgebra
using ProgressMeter
using Statistics
using ITensors: Index, QN, Out, In, dag
using QuadGK
using Random
using Distributions

##################################### Functions #############
"""
N -- number of sites

return an initialized ancilla of length N with conserved quantum numbers at infinite temperature.
"""
function AncillaMPO(N)
    s = ITensors.siteinds("S=1/2", N; conserve_qns=true)
    rho = MPO(s, "Id") ./ √2
    return rho, s
end

"""
ancilla -- MPS on which gates will be applied
h -- disorder
s -- Index of ancilla
op -- String, to choose which model you want (Heisenberg or XY)
δτ -- step of Trotter-Suzuki

return the vector of gates (ITensors.ITensor type)
"""
function gatesTEBDancilla(ancilla, h, δτ, s, op::String)
    N = length(ancilla)
    if op == "XX"
        gates = ops([("exp-τXY", (n, n + 1), (τ=δτ / 2, h=h)) for n in 1:1:(N - 1)], s)
        append!(gates, reverse(gates))
        return gates
    elseif op == "SS"
        gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h)) for n in 1:1:(N - 1)], s)
        append!(gates, reverse(gates))
        return gates
    end
end

#==
it would be better to use the Julia strengh and define two funciton type depending on the operator type but this implies to change every function in all files 

function gatesTEBDancilla(ancilla, h, δτ, s,::OpName"XX")
    N = length(ancilla)
    gates = ops([("exp-τXY", (n, n + 1), (τ=δτ / 2, h=h)) for n in 1:1:(N - 1)], s)
    append!(gates, reverse(gates))
    return gates
end

function gatesTEBDancilla(ancilla, h, δτ, s,::OpName"SS")
    N = length(ancilla)
    gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h)) for n in 1:1:(N - 1)], s)
    append!(gates, reverse(gates))
    return gates
end
==#

"""
ancilla -- MPS on which gates will be applied
h -- disorder
s -- Index of ancilla
op -- String, to choose which model you want (Heisenberg or XY)
δτ -- step of Trotter-Suzuki
cutoff -- cutoff in the SVD when gates are applied
beta -- temperature goal

return the updated MPS at the temperature beta
"""
function TEBDancilla!(ancilla, gates, beta, cutoff, δτ, Dmax)
    if δτ <= 0
        return "δτ must be non negative"
    end
    k = floor(Int, beta / δτ)
    begin_time = time()
    for _ in 1:k
        ancilla = apply(gates, ancilla; cutoff, maxdim=Dmax)
        ancilla = ancilla / tr(ancilla)
    end
    end_time = time()
    println("time to evolve to $beta: ", end_time - begin_time, " s")
    return ancilla
end

"""
ancilla -- MPS
H -- MPO operator of the energy (Hamiltonian)

return Tr(ancilla*H)
"""
function energyMPO(ancilla, H)
    update = ancilla
    return inner(update, H)
end

"""
β -- inverse temperature
h -- disorder
y -- proportion of XX and YY in the model

exact energy at temperature beta for XY model at temperature β with disorder h 
"""
function exactenergyXX(β, h; y=0.0)
    function ε(k, h, y)
        return sqrt((cos(k) - h)^2 + (y * sin(k))^2)
    end
    integrand(k) = ε(k, h, y) * tanh(0.5 * β * ε(k, h, y)) / (2 * pi)
    val, _ = quadgk(integrand, -pi, pi; rtol=1e-9)
    return -val / 2
end

"""
spectre -- vap de l'Hamiltonien obtained by exact diagonalization
β -- inverse temperature
L -- number of sites in the Heinsenberg chain

return the exact energy from the exact diagonalization, could be for any Hamiltonian
"""
function energyexact(spectre, β, L)
    weights = exp.(-β .* spectre)
    Z = sum(weights)
    return sum(spectre .* weights) / (Z * L)
end

"""
return typeof::Vector{ITensor} =! ITensorMPS.MPO
"""
function ITensors.op(::OpName"exp-τSSdisorder", ::SiteType"S=1/2", s1::Index, s2::Index; h)
    H =
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2) +
        op("Sz", s1) * op("Sz", s2) +
        h * (op("Sz", s1) * op("Id", s2) + op("Id", s1) * op("Sz", s2))
    return H
end

"""
rng -- kernel to generate the random numbers 
ancilla -- MPS
s -- index sites
h -- disorder magnitude

return gates to measure and evolve your MPS with one fixed disorder magnitude 
"""
function evolutionwithrandomdisordergates(rng, ancilla, s, h::Float64, δτ)
    N = length(ancilla)
    if h < 0
        return "erreur, magnitude disorder has to be > 0"
    elseif h == 0
        disorder = [0 for i in 1:(N - 1)]
    elseif h > 0
        disorder = rand(rng, Uniform(-h, h), N - 1)  # utilise ce générateur local fixé
    end
    gatesmeasure = ops([("exp-τSSdisorder", (n, n + 1), (h=disorder[n],)) for n in 1:1:(N - 1)], s)
    gatesevolve = exp.(-δτ / 2 .* gatesmeasure)
    append!(gatesevolve, reverse(gatesevolve))
    return gatesmeasure, gatesevolve, disorder
end

"""
ancilla -- MPS
s -- index sites
h -- Vector of disorders (for each site)

return evolution and measure gates when you already the random disorder vector, to reproduce the results with a fixed and consistant vector
"""
function evolutionwithrandomdisordergates(ancilla, s, h::Vector, δτ)
    N = length(ancilla)
    gatesmeasure = ops([("exp-τSSdisorder", (n, n + 1), (h=h[n],)) for n in 1:1:(N - 1)], s)
    gatesevolve = exp.(-δτ / 2 .* gatesmeasure)
    append!(gatesevolve, reverse(gatesevolve))
    return gatesmeasure, gatesevolve
end

# ======================== DMRG ========================
"""
psi0 -- initial mps 
H -- operator on which you optimise : <psi/H/psi>/<psi/psi>
n_sweep -- number of sweep you'll optimize the function
dmax -- maximal bon dimension in the mps
cutoff -- cutoff in the svd process of the optimization
noise -- acceptance noise

works for mps only and return the ground state for the operator H and the ground value of the operator H
"""
function groundstateDMRG(psi0, H, n_sweep, dmax, cutoff, noise)
    n = length(psi0)
    sweeps = Sweeps(n_sweep)
    maxdim!(sweeps, dmax)
    cutoff!(sweeps, cutoff)
    noise!(sweeps, noise)
    energy, psi = dmrg(H, psi0, sweeps)
    return psi, energy / n
end

"""
x -- point
h -- upper boundary of the disorder magnitude
y -- pic center index
o -- pic width
H -- amplitude of the pic

define the disorder distribution: normal distribution on hard interval [0, h], h is the max disorder magnitude
"""
function samplingdisorder(x, h, y, o, H)
    if x < 0 || x > h
        return 0.0
    end
    return 1 + H * exp(-((x - y)^2) / (2o^2))
end

"""
N -- Int length of the sample
X -- Float64 maximum disorder magnitude
y -- Float64 pic center
init -- seed

return the disorder sample
"""
function rejection_sample(N::Int, X, y, init; o=0.005X, A=5.0)
    samples = Float64[]
    rng = MersenneTwister(init)

    # Trouver une borne supérieure sur la densité (pour normalisation du rejet)
    max_density = 1 + A  # p(x) <= 1 + A partout

    while length(samples) < N
        x = rand(rng) * X  # tirage uniforme sur [0, X]
        u = rand(rng) * max_density
        if u < samplingdisorder(x, X, y, o, A)
            push!(samples, x)
        end
    end
    return round.(samples, digits=4)
end

"""
energylist -- list of the energy
betalist -- list of beta

return the gradient of a list (energylist) with respect to another list (betalist) by finite differences method
"""
function specificheat(energylist, betalist)
    n, m = length(energylist), length(betalist)
    @assert n == m "Gradient not possible because two different lengths"
    Grandientlist = diff(energylist) ./ diff(betalist)
    return betalist[begin:(end - 1)], Grandientlist
end