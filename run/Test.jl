#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots
using ITensorMPS
using Statistics
using Base.Threads
using LaTeXStrings
println("Nombre de threads disponibles : ", nthreads())
using JSON
################ Parameters ###############

N = 18
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
Dmax = 300
betamax = 5
step = 0.1
Beta = n_sweep * δτ
gammescale = 0.8
j = "z"
γ = 0.0
betalist = collect(0:step:betamax)
filename = joinpath("analyse_simulations_julia", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

# ========================= FUNCTIONS

function energyexact(spectre, β, L)
    weights = exp.(-β .* spectre)
    Z = sum(weights)
    return sum(spectre .* weights) / (Z * L)
end

function energysitetest(mps, sitemeasure, h)
    copy = orthogonalize(mps, sitemeasure)
    sn = siteind(copy, sitemeasure)
    snn = siteind(copy, sitemeasure + 1)
    gate =
        1 / 2 * op("S+", sn) * op("S-", snn) +
        1 / 2 * op("S-", sn) * op("S+", snn) +
        #op("Sz", sn) * op("Sz", snn) +
        h * (op("Sz", sn) * op("Id", snn) + op("Id", sn) * op("Sz", snn))
    inter = copy[sitemeasure] * copy[sitemeasure+1]
    normalize!(inter)
    #@show gate
    prout = replaceprime(gate, 0 => 2)
    #@show replaceprime(gate, 0=>2)
    #@show dag(inter)
    #@show prout*inter
    double = replaceprime(prout * inter, 2 => 1)
    #@show scalar(double*dag(inter))
    return real(scalar(double * dag(inter)))
end

function energyagainstsitetest(mps, h, scale)
    N = length(mps)
    start, stop = MBL.section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    #@show sites
    Energypersite = Vector(undef, length(sites))
    @showprogress desc = "calcul energy over sites" for i in eachindex(sites)
        #@show i 
        Energypersite[i] = energysitetest(mps, sites[i], h)
    end
    return sites, mean(Energypersite)
end


function energyforbetalisttest(betamax, step, ancilla, δτ, h, s, cutoff, op)
    betalist = collect(0:step:betamax)
    realbetalist = reverse(push!(diff(betalist), 0))
    Energylist = Vector{}(undef, length(realbetalist))
    update = ancilla
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, δτ, h, realbetalist[i], s, cutoff, op)
        _, Energylist[i] = energyagainstsitetest(update, h, gammescale)
    end
    return betalist, Energylist
end

# ============================== DATA
test, s = MBL.AncillaMPO(N)

xdatasites, ydatadites = energyforbetalisttest(betamax, step, test, δτ, h, s, cutoff, "XY")
xdataMPO, ydataMPO = MBL.energyforbetalist(betamax, step, test, δτ, h, s, cutoff, "XY")

exactenergy = [MBL.exactenergyXY(β, h, γ) for β in xdataMPO]

#exactdz = [energyexact(input["spectrum"], beta, N) for beta in xdataMPO]

gr()
scatter(xdatasites, ydatadites/2, label="mesure par site/2")
scatter!(xdatasites, ydatadites, label="mesure par site")
scatter!(xdataMPO, ydataMPO / N, label="mesure MPO")
#plot!(xdataMPO, exactdz, label="exact dz", xlabel="β", ylabel="<H>/N", title="N=$N, δτ=$δτ, cutoff=$cutoff, model XY")
plot!(xdataMPO, exactenergy, label=latexstring("-1/4π\\int_{-π}^{π} cos(k)tanh(βcos(k)/2)dk"), xlabel="β", ylabel="<H>/N", title="N=$N, δτ=$δτ, cutoff=$cutoff, model XY")
#hline!([1/4-log(2)], label="exact energy at zero temperature")