#!usr/bin/env julia

################## Librairies ####################
using ITensors
using ITensorMPS
using LinearAlgebra
using ProgressMeter
using Statistics
################## Functions #####################

# ==================================================== Init
"""
s -- sites index
D -- bond dimension 

return an randmoly initialized MPS with the sites indexes in input 
"""
function random_initialized_MPS(N, D)
    s = ITensors.siteinds("S=1/2", N)
    psi = random_mps(s, linkdims=D)
    return psi, s
end

"""
N -- number of sites

return a Neel state of the form |↑↓↑↓...> with N sites
"""
function neelstate(N)
    s = siteinds("S=1/2", N, conserve_qns=true)
    mps = MPS(s, n -> isodd(n) ? "Up" : "Dn")
    return mps, s
end

# ==================================================== TEBD
"""
N -- number of sites 
J -- coupling constant
h -- disorder constant

Define the operator expτ(-JSS - hSz)
"""
function ITensors.op(::OpName"exp-τSS", ::SiteType"S=1/2", s1::Index, s2::Index; τ, h)
    H =
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2) +
        op("Sz", s1) * op("Sz", s2) +
        h * (op("Sz", s1) * op("Id", s2) + op("Id", s1) * op("Sz", s2))
    return exp(-τ * H)
end

"""
N -- number of sites 
J -- coupling constant
h -- disorder constant

Define the operator expτ(-JSS - hSz)
"""
function ITensors.op(::OpName"exp-τXY", ::SiteType"S=1/2", s1::Index, s2::Index; τ, h)
    H =
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2) +
        h * (op("Sz", s1) * op("Id", s2) + op("Id", s1) * op("Sz", s2))
    return exp(-τ * H)
end

"""
return the vector of Trotter Suzuki gates for Heinsenberg model in a row that means gates are in the order: (1,2) ; (2,3) ; ...
"""
function gateTrotterSuzukirow(mps, h, δτ)
    N = length(mps)
    s = siteinds(mps)
    gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:1:(N-1)], s)
    append!(gates, reverse(gates))
    return gates
end

"""
return the vector of Trotter Suzuki gates for XY model in a row that means gates are in the order: (1,2) ; (2,3) ; ...
"""
function gateTrotterSuzukirowXY(mps, h, δτ)
    N = length(mps)
    s = siteinds(mps)
    gates = ops([("exp-τXY", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:1:(N-1)], s)
    append!(gates, reverse(gates))
    return gates
end

"""
return the converged mps with the row application of gates
"""
function tebdstepHeisenbergRow!(nsweep, mps, h, δτ, cutoff, Dmax, op::String)
    if op == "SS"
        gate = gateTrotterSuzukirow(mps, h, δτ)
    elseif op == "XY"
        gate = gateTrotterSuzukirowXY(mps, h, δτ)
    end
    @showprogress desc = "tebd steps" for i in 1:nsweep
        mps = apply(gate, mps; cutoff, maxdim=Dmax)
        normalize!(mps)
    end
    return mps
end


# ==================================================== Operators
"""
mps -- mps on which you compute the energy
h -- disorder

return the Heisenberg Hamiltonian with disorder with the ITensorMPS.MPO type 
"""
function hamiltonianHeisenberg(mps, h::Float64, s)
    N = length(mps)
    ampo = AutoMPO()
    for j in 1:N-1
        add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
        add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
        add!(ampo, 1, "Sz", j, "Sz", j + 1)
        add!(ampo, h, "Sz", j)
    end
    add!(ampo, h, "Sz", N)
    H = MPO(ampo, s)
    return H
end

"""
mps -- mps on which you compute the energy
h -- disorder

return the Heisenberg Hamiltonian with disorder with the ITensorMPS.MPO type 
"""
function hamiltonianHeisenberg(mps, h::Vector, s)
    N = length(mps)
    ampo = AutoMPO()
    for j in 1:N-1
        add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
        add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
        add!(ampo, 1, "Sz", j, "Sz", j + 1)
        add!(ampo, h[j], "Sz", j)
    end
    add!(ampo, h[N], "Sz", N)
    H = MPO(ampo, s)
    return H
end

"""
mps -- mps on which you compute the energy
h -- disorder

return the Heisenberg Hamiltonian with disorder with the ITensorMPS.MPO type 
"""
function hamiltonianXY(mps, h, s)
    N = length(mps)
    ampo = AutoMPO()
    for j in 1:N-1
        add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
        add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
        add!(ampo, h, "Sz", j)
    end
    add!(ampo, h, "Sz", N)
    H = MPO(ampo, s)
    return H
end


"""
exact energy of the 1D Heisenberg Hamiltonian ground state
"""
function exactgroundenergy(J=1)
    return J * (1 / 4 - log(2))
end

# ============================================ Measure 
"""
psi -- MPS converged on which you make the measurement 
n -- site measure

return the Sz value on the site n 
"""
function measure_S(psi::MPS, n, j::String)
    string = "S" * j
    psicopy = orthogonalize(psi, n)
    sn = siteind(psicopy, n)
    Sz = scalar(dag(prime(psicopy[n], "Site")) * op(string, sn) * psicopy[n])
    return real(Sz)
end

"""
psi -- MPS converged on which you make the measurement 
n -- site measure

return the Sz value on the site n 
"""
function measure_S(psi::MPO, n, j::String)
    string = "S" * j
    psicopy = orthogonalize(psi, n)
    sn = siteind(psicopy, n)
    gate = replaceprime(op(string, sn), 1 => 2)
    inter = replaceprime(gate * psicopy[n], 2 => 0)
    Sz = scalar(dag(psicopy[n]) * inter)
    return real(Sz)
end

"""
psi -- MPS converged on which you make the measurement 
n -- site measure
H -- hamiltonian

WARNING: psi and H have to be compatible that means they have to be generated by the functions in this file

return the expectation value of H on site n
"""
function measure_H(psi, n, H)
    orthogonalize!(psi, n)
    e = inner(psi, H, psi)
    return e
end

"""
mps -- converged mps with tebd you want to access to the energy at the link sitemeasure
sitemeasure -- index of the site

return the energy on the site sitemeasure
"""
function energysite(mps, sitemeasure, h)
    copy = orthogonalize(mps, sitemeasure)
    sn = siteind(copy, sitemeasure)
    snn = siteind(copy, sitemeasure + 1)
    gate =
        1 / 2 * op("S+", sn) * op("S-", snn) +
        1 / 2 * op("S-", sn) * op("S+", snn) +
        op("Sz", sn) * op("Sz", snn) +
        h * (op("Sz", sn) * op("Id", snn) + op("Id", sn) * op("Sz", snn))
    inter = copy[sitemeasure] * copy[sitemeasure+1]
    normalize!(inter)
    e = scalar(dag(prime(inter, "Site")) * gate * inter)
    return real(e)
end

"""
mps -- MPS
sitemeasure -- index in the mps of the site you want to compute the energy
h -- disorder

return the energy of mps at the site sitemeasure 
"""
function energysiteMPO(mps, sitemeasure, h, operateur::String)
    copy = orthogonalize(mps, sitemeasure)
    sn = siteind(copy, sitemeasure)
    snn = siteind(copy, sitemeasure + 1)
    if operateur == "XY"
        gate =
            1 / 2 * op("S+", sn) * op("S-", snn) +
            1 / 2 * op("S-", sn) * op("S+", snn) +
            h * (op("Sz", sn) * op("Id", snn) + op("Id", sn) * op("Sz", snn))
    elseif operateur == "SS"
        gate =
            1 / 2 * op("S+", sn) * op("S-", snn) +
            1 / 2 * op("S-", sn) * op("S+", snn) +
            op("Sz", sn) * op("Sz", snn) +
            h * (op("Sz", sn) * op("Id", snn) + op("Id", sn) * op("Sz", snn))
    end
    inter = copy[sitemeasure] * copy[sitemeasure+1]
    normalize!(inter)
    adjust = replaceprime(gate, 0 => 2)
    double = replaceprime(adjust * inter, 2 => 1)
    return real(scalar(double * dag(inter)))
end

"""
mps -- MPS
sitemeasure -- index in the mps of the site you want to compute the energy
gate -- gates with random disorder

return the energy of mps at the site sitemeasure 
"""
function energysiteMPOdisorder(mps, sitemeasure, gate)
    newgate = replaceprime(gate, 0 => 2)
    copy = orthogonalize(mps, sitemeasure)
    inter = copy[sitemeasure] * copy[sitemeasure+1]
    normalize!(inter)
    adjust = replaceprime(inter, 1 => 2)
    double = newgate * adjust
    return real(scalar(double * dag(inter)))
end

"""
j -- axis of the spin operator 
n -- starting site
p -- the other site, the distance between the two sites is p-n+1

return the correlation function
"""
function correlationSpinoperator(mps, n, p, j)
    copy = orthogonalize(mps, n)
    sn = siteind(copy, n)
    sp = siteind(copy, p)
    if j == "z"
        S_n = op("Sz", sn)
        S_p = op("Sz", sp)
        gate = S_n * S_p
        psi_np = apply(gate, copy)
        return real(inner(copy, psi_np))
    elseif j == "x"
        S_n = 1 / 2 * (op("S+", sn) + op("S-", sn))
        S_p = 1 / 2 * (op("S+", sp) + op("S-", sp))
        gate = S_n * S_p
        psi_np = apply(gate, copy)
        return real(inner(copy, psi_np))
    elseif j == "y"
        S_n = -1im / 2 * (op("S+", sn) - op("S-", sn))
        S_p = -1im / 2 * (op("S+", sp) - op("S-", sp))
        gate = S_n * S_p
        psi_np = apply(gate, copy)
        return real(inner(copy, psi_np))
    end
end

"""
k -- length of the segment 
j -- axis against sprin is measured

the two boundary of the chain are excluded

return the correlation function over
"""
function correlationonlength(mps, k, j)
    N = length(mps)
    if k > N - 2
        return "mps is not long enough"
    end
    Listintercorrel = Vector{Float64}(undef, max(1, N - k - 2)) #car quand k = N-2 on veut quand même une valeur de correlation
    for i in 2:1:max(N - k - 1, 2) #car quand k = N-2 on a quand meme un point pour calculer la correlation
        #@show correlationSpinoperator(mps, i, i + k, j)
        Listintercorrel[i-1] = correlationSpinoperator(mps, i, i + k, j)
    end
    return mean(Listintercorrel)
end

"""
mps -- MPS
h -- disorder 
scale -- 0<scale<1, pourcentage of the chain you want to measure (from the middle chain)

return the MPS average energy for the model op, measured with gates 
"""
function energyagainstsiteMPO(mps, h, scale, op::String)
    N = length(mps)
    start, stop = MBL.section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    #@show sites
    Energypersite = Vector(undef, length(sites))
    @showprogress desc = "calcul energy over sites" for i in eachindex(sites)
        #@show i 
        Energypersite[i] = energysiteMPO(mps, sites[i], h, op)
    end
    return sites, mean(Energypersite)
end

"""
mps -- MPS
h -- disorder 
scale -- 0<scale<1, pourcentage of the chain you want to measure (from the middle chain)

return the MPS average energy for the model op, measured with gates 
"""
function energyagainstsiteMPOdisorder(mps, gates, scale)
    N = length(mps)
    start, stop = MBL.section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    #@show sites
    Energypersite = Vector(undef, length(sites))
    update = mps
    @showprogress desc = "calcul energy over sites" for i in eachindex(sites)
        #@show i 
        Energypersite[i] = energysiteMPOdisorder(update, sites[i], gates[sites[i]])
    end
    return sites, mean(Energypersite)
end

"""
mps -- MPS
h -- disorder 
scale -- 0<scale<1, pourcentage of the chain you want to measure (from the middle chain)

return the MPS average energy for the model op, measured with gates 
"""
function magnetagainstsiteMPOdisorder(mps, gates, scale)
    N = length(mps)
    start, stop = MBL.section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    #@show sites
    Magnetpersite = Vector(undef, length(sites))
    update = mps
    @showprogress desc = "calcul energy over sites" for i in eachindex(sites)
        #@show i 
        Magnetpersite[i] = energysiteMPOdisorder(update, sites[i], gates[sites[i]])
    end
    return sites, Magnetpersite
end

# =========================================== Else
"""
return the max bond dimension in the mps

"""
function maxbonddim(mps)
    maxdim = 0
    #@show typeof(mps)
    for i in 1:(length(mps)-1)
        #@show i
        s = commonind(mps[i], mps[i+1])
        #@show s
        if s === nothing
            continue
        end
        maxdim = max(maxdim, ITensors.dim(s))
    end
    return maxdim
end

"""
N -- length of a vector
scale -- between 0 and 1, overlap between the sliced list and the original list

return the indexes that slices a list of length N with the overlap scale
"""
function section_trunc(N, scale)
    q = div(N, 2)
    be, st = max(floor(Int, (1 + (1 - scale) * q)), 1), min(floor(Int, ((scale + 1) * q)), N)
    return be, st
end

"""
return a full random operator 
"""
function randomoperator(T, s1, s2)
    return ITensor(T, s1, s2, dag(s1), dag(s2))
end