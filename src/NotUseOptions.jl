#!usr/bin/env julia

using ITensors
using ITensorMPS
using ProgressMeter
using QuadGK

"""
N -- number of sites 
J -- coupling constant
h -- disorder constant

return the Trotter Suzuki gates (order 2) and the Hamiltonian to compute TEBD and energy
"""
function gateTrotterSuzukiandhamiltonian(mps, h, δτ, parity::String)#essayer l'autre contraction où on fait pair impair à la suite au lieu de tout pair puis impair
    N = length(mps)
    s = siteinds(mps)
    if mod(N, 2) == 0
        if parity == "even"
            gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 1:2:(N-1)], s)
        elseif parity == "odd"
            gates = ops([("exp-τSS", (n, n + 1), (τ=δτ, h=h,)) for n in 2:2:(N-2)], s)
        end
    elseif mod(N, 2) == 1
        if parity == "even"
            gates = ops([("exp-τSS", (n, n + 1), (τ=δτ / 2, h=h,)) for n in 2:2:(N-2)], s)
        elseif parity == "odd"
            gates = ops([("exp-τSS", (n, n + 1), (τ=δτ, h=h,)) for n in 1:2:(N-1)], s)
        end
    end
    #@show typeof(gates)
    #append!(gates, reverse(gates))
    return gates
end


"""
nsweep -- number of sweeps
mps -- initial mps
gates -- even and odd you want to apply
cutoff -- cutoff in the truncation part in the applying gate process

return the converged mps with n sweeps along the mps
"""
function tebdstepHeisenberg!(nsweep, mps, h, δτ, cutoff, Dmax)
    @showprogress for i in 1:nsweep
        gateseven = gateTrotterSuzukiandhamiltonian(mps, h, δτ, "even")
        mps = apply(gateseven, mps; cutoff, maxdim=Dmax)
        normalize!(mps)
        gatesodd = gateTrotterSuzukiandhamiltonian(mps, h, δτ, "odd")
        mps = apply(gatesodd, mps; cutoff, maxdim=Dmax)
        normalize!(mps)
        gateseven = gateTrotterSuzukiandhamiltonian(mps, h, δτ, "even")
        mps = apply(gateseven, mps; cutoff, maxdim=Dmax)
        normalize!(mps)
    end
    return mps
end

"""
j -- axis of the spin operator 
n -- starting site
p -- the other site, the distance between the two sites is p-n+1

return the correlation function
"""
function correlationSpinoperatorwrong(mps, n, p, j)
    copy = orthogonalize(mps, n)
    sn = siteind(copy, n)
    sp = siteind(copy, p)
    if j == "z"
        S_n = op("Sz", sn)
        S_p = op("Sz", sp)
        psi_n = apply(S_n, copy; site=n)
        psi_np = apply(S_p, psi_n; site=p)
        return inner(copy, psi_np)
    elseif j == "x"
        S_n = 1 / 2 * (op("S+", sn) + op("S-", sn))
        S_p = 1 / 2 * (op("S+", sp) + op("S-", sp))
        psi_n = apply(S_n, copy; site=n)
        psi_np = apply(S_p, psi_n; site=p)
        return inner(copy, psi_np)
    elseif j == "y"
        S_n = -1im / 2 * (op("S+", sn) - op("S-", sn))
        S_p = -1im / 2 * (op("S+", sp) - op("S-", sp))
        psi_n = apply(S_n, copy; site=n)
        psi_np = apply(S_p, psi_n; site=p)
        return inner(copy, psi_np)
    end
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