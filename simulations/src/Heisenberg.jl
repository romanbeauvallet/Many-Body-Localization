#!usr/bin/env julia

################## Librairies ####################

using ITensors
using ITensorMPS
using LinearAlgebra
using ProgressMeter
using Statistics
# tune; 1–4 often best
################## Functions #####################

# ============================================= Initialization
"""
s -- sites index
D -- bond dimension 

return an randmoly initialized MPS with the sites indexes in input 
"""
function random_initialized_MPS(N, D)
    s = ITensors.siteinds("S=1/2", N)
    psi = random_mps(s; linkdims=D)
    return psi, s
end

"""
N -- number of sites

return a Neel state of the form |↑↓↑↓...> with N sites
"""
function neelstate(N)
    s = siteinds("S=1/2", N; conserve_qns=true)
    mps = MPS(s, n -> isodd(n) ? "Up" : "Dn")
    return mps, s
end

# ============================================= TEBD
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
function gateTrotterSuzukirow(mps, h, δτ, op::String)
    N = length(mps)
    s = siteinds(mps)
    if op == "SS"
        gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h)) for n in 1:1:(N - 1)], s)
        append!(gates, reverse(gates))
        return gates
    elseif op == "XX"
        gates = ops([("exp-τXY", (n, n + 1), (τ=δτ / 2, h=h)) for n in 1:1:(N - 1)], s)
        append!(gates, reverse(gates))
        return gates
    end
end

"""
return the converged mps with the row application of gates
"""
function tebdevolutionrow!(nsweep, mps, gates, cutoff, Dmax)
    @showprogress desc = "tebd steps" for i in 1:nsweep
        mps = apply(gates, mps; cutoff, maxdim=Dmax)
        normalize!(mps)
    end
    return mps
end

# ============================================= Hamiltonians
"""
mps -- mps on which you compute the energy
h -- disorder
s -- index of the mps
op -- choose between the Heisenberg and the XY model

return the Heisenberg Hamiltonian with disorder with the ITensorMPS.MPO type 
"""
function operator(mps, h::Float64, s, op::String)
    N = length(mps)
    ampo = AutoMPO()
    if op == "SS"
        for j in 1:(N - 1)
            add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
            add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
            add!(ampo, 1, "Sz", j, "Sz", j + 1)
            add!(ampo, h, "Sz", j)
        end
        add!(ampo, h, "Sz", N)
        H = MPO(ampo, s)
        return H
    elseif op == "XX"
        for j in 1:(N - 1)
            add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
            add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
            add!(ampo, h, "Sz", j)
        end
        add!(ampo, h, "Sz", N)
        H = MPO(ampo, s)
        return H
    end
end

"""
mps -- mps on which you compute the energy
h -- disorder
s -- index of the mps
op -- choose between the Heisenberg and the XY model

return the Heisenberg Hamiltonian with disorder with the ITensorMPS.MPO type 
"""
function operator(mps, h::Vector, s, op::String)
    N = length(mps)
    ampo = AutoMPO()
    if op == "SS"
        for j in 1:(N - 1)
            add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
            add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
            add!(ampo, 1, "Sz", j, "Sz", j + 1)
            add!(ampo, h[j], "Sz", j)
        end
        add!(ampo, h[N], "Sz", N)
        H = MPO(ampo, s)
        return H
    elseif op == "XX"
        for j in 1:(N - 1)
            add!(ampo, 1 / 2, "S+", j, "S-", j + 1)
            add!(ampo, 1 / 2, "S-", j, "S+", j + 1)
            add!(ampo, h[j], "Sz", j)
        end
        add!(ampo, h[N], "Sz", N)
        H = MPO(ampo, s)
        return H
    end
end

"""
exact energy of the 1D Heisenberg Hamiltonian ground state
"""
function exactgroundenergy(J=1)
    return J * (1 / 4 - log(2))
end

# ============================================= Measure
"""
psi -- MPS converged on which you make the measurement 
n -- site measure
j -- choose the axis (WARNING x and y doesnt work with quantum numbers)

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
j -- choose the axis (WARNING x and y doesnt work with quantum numbers)

return the Sz expectation value on the site n 
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
function energysite(mps::MPS, sitemeasure, h, operator::String)
    copy = orthogonalize(mps, sitemeasure)
    sn = siteind(copy, sitemeasure)
    snn = siteind(copy, sitemeasure + 1)
    if operator == "SS"
        gate =
            1 / 2 * op("S+", sn) * op("S-", snn) +
            1 / 2 * op("S-", sn) * op("S+", snn) +
            op("Sz", sn) * op("Sz", snn) +
            h * (op("Sz", sn) * op("Id", snn) + op("Id", sn) * op("Sz", snn))
    elseif operator == "XX"
        gate =
            1 / 2 * op("S+", sn) * op("S-", snn) +
            1 / 2 * op("S-", sn) * op("S+", snn) +
            h * (op("Sz", sn) * op("Id", snn) + op("Id", sn) * op("Sz", snn))
    end
    inter = copy[sitemeasure] * copy[sitemeasure + 1]
    normalize!(inter)
    e = scalar(dag(prime(inter, "Site")) * gate * inter)
    return real(e)
end

"""
mps -- MPO
sitemeasure -- index in the mps of the site you want to compute the energy
h -- disorder

return the energy of mps at the site sitemeasure 
"""
function energysite(mps::MPO, sitemeasure, h, operateur::String)
    copy = orthogonalize(mps, sitemeasure)
    sn = siteind(copy, sitemeasure)
    snn = siteind(copy, sitemeasure + 1)
    if operateur == "XX"
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
    inter = copy[sitemeasure] * copy[sitemeasure + 1]
    normalize!(inter)
    adjust = replaceprime(gate, 0 => 2)
    double = replaceprime(adjust * inter, 2 => 1)
    return real(scalar(double * dag(inter)))
end

"""
mps -- MPS or MPO
h -- disorder
scale -- percentage of the chain on which you measure
op -- choose between Heisenberg Hamiltonian and XY model

return the site list and the energy per site 
"""
function energyagainstsite(mps, h, scale, op::String)
    N = length(mps)
    start, stop = section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    Energypersite = Array{Float64}(undef, stop - start + 1)
    for i in eachindex(sites)
        Energypersite[i] = energysite(mps, sites[i], h, op)
    end
    return sites, Energypersite
end

"""
mps -- MPO
sitemeasure -- index in the mps of the site you want to compute the energy
gate -- gates with random disorder

return the energy of mps at the site sitemeasure 
"""
function energysitedisorder(mps::MPO, sitemeasure, gate)
    newgate = replaceprime(gate, 0 => 2)
    copy = orthogonalize(mps, sitemeasure)
    inter = copy[sitemeasure] * copy[sitemeasure + 1]
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
        Listintercorrel[i - 1] = correlationSpinoperator(mps, i, i + k, j)
    end
    return mean(Listintercorrel)
end

# ============================================= Else
"""
mps -- mps you want to have the max bon dimension 

return the max bond dimension in the mps
"""
function maxbonddim(input)
    maxdim = 0
    mps = input
    for i in 1:(length(mps) - 1)
        s = commonind(mps[i], mps[i + 1])
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
T -- a random array 
s1 -- index 
s2 -- index

return a full random operator 
"""
function randomoperator(T, s1, s2)
    return ITensor(T, s1, s2, dag(s1), dag(s2))
end