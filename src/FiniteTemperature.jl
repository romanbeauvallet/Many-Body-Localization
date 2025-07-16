#!usr/bin/env julia

################## Librairies ####################
using ITensors
using ITensorMPS
using LinearAlgebra
using ProgressMeter
using Statistics
using ITensors: Index, QN, Out, In, dag
using Printf
using QuadGK
using Random
using Distributions
################## Functions #####################

"""
N -- number of sites

return an initialized ancilla of length N 
"""
function AncillaMPO(N)
    s = ITensors.siteinds("S=1/2", N; conserve_qns=true)
    rho = MPO(s, "Id") ./ √2
    return rho, s
end

"""
N -- number of sites

return the ancilla mps with quantum numbers and Bell state between auxiliary and physical indexes
"""
function Ancilla(N) #à corriger car sinteinds de Vector{ITensor ne marche pas}
    phys_sites = siteinds("S=1/2", N; conserve_qns=true)
    aux_sites = siteinds("S=1/2", N; conserve_qns=true)

    # Définir les liens avec directions : le premier en Out, le reste en In/Out alterné
    links = [Index([QN("Sz", +1) => 1, QN("Sz", -1) => 1], "Link, n=$i"; dir=Out) for i in 1:N+1]

    mps_tensors = Vector{ITensor}(undef, N)

    for i in 1:N
        l = dag(links[i])       # dir = In
        r = links[i+1]          # dir = Out
        s = phys_sites[i]       # dir = Out (par défaut)
        a = dag(aux_sites[i])   # dir = In (on suppose que c’est un "bra")

        A = ITensor(l, r, s, a)  # flux total = qn(r) + qn(s) - qn(l) - qn(a)

        # On remplit uniquement les blocs conservant le QN total
        A[l=>1, r=>1, s=>"Up", a=>"Up"] = 1.0   # Sz = +1
        A[l=>2, r=>2, s=>"Dn", a=>"Dn"] = 1.0   # Sz = -1

        mps_tensors[i] = A
    end

    return mps_tensors
end

"""
ancilla -- MPS on which gates will be applied
h -- disorder
s -- Index of ancilla
op -- String, to choose which model you want (Heisenberg or XY)
δτ -- step of Trotter-Suzuki

return the vector of gates (ITensor)
"""
function gatesTEBDancilla(ancilla, h, δτ, s, op::String)
    N = length(ancilla)
    if op == "XY"
        gates = ops([("exp-τXY", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:1:(N-1)], s)
        append!(gates, reverse(gates))
        return gates
    elseif op == "SS"
        gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:1:(N-1)], s)
        append!(gates, reverse(gates))
        return gates
    end
end

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
function TEBDancilla!(ancilla, gates, beta, cutoff, δτ)
    for β in 0:δτ:beta
        ancilla = apply(gates, ancilla; cutoff)
        #@printf("β = %.2f energy = %.8f\n", β, energyancilla)
        ancilla = ancilla / tr(ancilla)
    end
    return ancilla
end

"""
ancilla -- MPS
H -- MPO operator of the energy (Hamiltonian)

return the Tr(ancilla*H)
"""
function energyMPO(ancilla, H)
    return inner(ancilla, H)
end

"""
β -- inverse temperature
h -- disorder
γ -- proportion of XX and YY in the model
exact energy at temperature beta for XY model at temperature β with disorder h 
"""
function exactenergyXY(β, h, γ=0.0)
    function ε(k, h, γ)
        return sqrt((cos(k) - h)^2 + (γ * sin(k))^2)
    end
    integrand(k) = ε(k, h, γ) * tanh(0.5 * β * ε(k, h, γ)) / (2 * pi)
    val, _ = quadgk(integrand, -pi, pi, rtol=1e-9)
    return -val / 2
end


function exactenergyXY3(β)
    integrand(k) = cos(k) * tanh(0.5 * β * cos(k))
    val, _ = quadgk(integrand, 0, pi, rtol=1e-9)
    return -val / (2 * pi)
end

function exactenergyXY2(β)
    integrand(k) = cos(k) * tanh(β * cos(k))
    val, _ = quadgk(integrand, 0, pi / 2, rtol=1e-9)
    return -2 * val / (pi)
end

"""
spectre -- vap de l'Hamiltonien obtained by exact diagonalization
β -- inverse temperature
L -- number of sites in the Heinsenberg chain

return the exact energy 
"""
function energyexact(spectre, β, L)
    weights = exp.(-β .* spectre)
    Z = sum(weights)
    return sum(spectre .* weights) / (Z * L)
end

function ITensors.op(::OpName"exp-τSSdisorder", ::SiteType"S=1/2", s1::Index, s2::Index; h)
    H =
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2) +
        op("Sz", s1) * op("Sz", s2) +
        h * (op("Sz", s1) * op("Id", s2) + op("Id", s1) * op("Sz", s2))
    return H
end

"""
init -- integer to init the seed
ancilla -- MPS
"""
function evolutionwithrandomdisordergates(init::Int64, ancilla, s, h, δτ)
    rng = MersenneTwister(init)
    N = length(ancilla)
    if h < 0
        return "erreur, magnitude disorder has to be > 0"
    elseif h == 0
        disorder = [0 for i in 1:N-1]
    elseif h > 0
        disorder = rand(rng, Uniform(-h, h), N - 1)  # utilise ce générateur local fixé
    end
    gatesmeasure = ops([("exp-τSSdisorder", (n, n + 1), (h=disorder[n],)) for n in 1:1:(N-1)], s)
    gatesevolve = exp.(-δτ .* gatesmeasure)
    append!(gatesevolve, reverse(gatesevolve))
    return gatesmeasure, gatesevolve
end
