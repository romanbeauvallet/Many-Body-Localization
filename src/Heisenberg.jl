#!usr/bin/env julia

################## Librairies ####################
using ITensors
using ITensorMPS
using LinearAlgebra

################## Functions #####################

function heisenberghamiltonian(J, h ,N)
    sites = ITensors.siteinds("SpinHalf", N)
    ampo = AutoMPO()
    for j in 1:N-1
        add!(ampo, -J/2, "S+", j, "S-", j+1)
        add!(ampo, -J/2, "S-", j, "S+", j+1)
        add!(ampo, -J, "Sz", j, "Sz", j+1)
        add!(ampo, -h, "Sz", j)
    end
    H = MPO(ampo, sites)
    return H
end