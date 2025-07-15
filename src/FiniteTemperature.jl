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
################## Functions #####################

"""
return an initialized ancilla of length N 
"""
function AncillaMPO(N)
    s = ITensors.siteinds("S=1/2", N; conserve_qns=true)
    rho = MPO(s, "Id") ./ √2
    return rho, s
end

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

"""
function gatesTEBDancilla(ancilla, h, δτ, s, op::String)
    N = length(ancilla)
    if op=="XY"
        gates = ops([("exp-τXY", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:1:(N-1)], s)
        append!(gates, reverse(gates))
        return gates
    elseif op=="SS"
        gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:1:(N-1)], s)
        append!(gates, reverse(gates))
        return gates
    end
end


"""

"""
function TEBDancilla!(ancilla, δτ, h, beta, s, cutoff, op)
    gates = gatesTEBDancilla(ancilla, h, δτ, s, op)
    for β in 0: δτ:beta
        ancilla = apply(gates, ancilla; cutoff)
        #@printf("β = %.2f energy = %.8f\n", β, energyancilla)
        ancilla = ancilla / tr(ancilla)
  end
  return ancilla
end

function energyMPO(ancilla, H)
    return inner(ancilla, H)
end

"""
exact energy at temperature beta for XY model

"""
function exactenergyXY1(β, h, γ)
    function ε(k, h, γ)
        return sqrt((cos(k) - h)^2 + (γ * sin(k))^2)
    end
    integrand(k) = ε(k, h, γ) * tanh(0.5 * β * ε(k, h, γ)) / (2*pi)
    val, _ = quadgk(integrand, -pi, pi, rtol=1e-9)
    return -val
end 

function exactenergyXY3(β)
    integrand(k) = cos(k) * tanh(0.5*β * cos(k))
    val, _ = quadgk(integrand, 0, pi, rtol=1e-9)
    return -val/(2*pi)
end 

function exactenergyXY2(β)
    integrand(k) = cos(k) * tanh(β * cos(k))
    val, _ = quadgk(integrand, 0, pi/2, rtol=1e-9)
    return -2*val/(pi)
end 
