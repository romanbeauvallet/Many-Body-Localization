#!usr/bin/env julia

################## Librairies ####################
using ITensors
using ITensorMPS
using LinearAlgebra

################## Functions #####################

"""
N -- number of sites 
J -- coupling constant
h -- disorder constant

Define the operator expτ(-JSS - hSz)
"""
function ITensors.op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ, h)
    H =
        -1 / 2 * op("S+", s1) * op("S-", s2) +
        -1 / 2 * op("S-", s1) * op("S+", s2)
    -h * op("Sz", s1) * op("Sz", s2)
    return exp(τ * H)
end

"""
N -- number of sites 
J -- coupling constant
h -- disorder constant

return the Trotter Suzuki gates and the Hamiltonian to compute TEBD and energy
"""
function gateTrotterSuzukiandhamiltonian(N, h, δτ)
    s = ITensors.siteinds("S=1/2", N)
    gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2, h=h,)) for n in 1:(N-1)], s)
    @show typeof(gates)
    append!(gates, reverse(gates))
    ampo = AutoMPO()
    for j in 1:N-1
        add!(ampo, -1 / 2, "S+", j, "S-", j + 1)
        add!(ampo, -1 / 2, "S-", j, "S+", j + 1)
        add!(ampo, -1 / 2, "Sz", j, "Sz", j + 1)
        add!(ampo, -h, "Sz", j)
    end
    add!(ampo, -h, "Sz", N)
    H = MPO(ampo, s)
    return gates, H
end
"""
exact energy of the 1D Heisenberg Hamiltonian ground state
"""
function exactgroundenergy(J=1)
    return J * (1 / 4 - log(2))
end
