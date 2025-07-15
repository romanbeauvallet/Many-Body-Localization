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
println("Nombre de threads disponibles : ", nthreads())

################ Parameters ###############

N = 10
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

function energysitetest(mps, sitemeasure, h)
    copy = orthogonalize(mps, sitemeasure)
    sn = siteind(copy, sitemeasure)
    snn = siteind(copy, sitemeasure + 1)
    gate =
        1 / 2 * op("S+", sn) * op("S-", snn) +
        1 / 2 * op("S-", sn) * op("S+", snn) +
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
################# try to add quantum number conservation
test, s = MBL.AncillaMPO(N)
#@show siteinds(test)
#@show typeof(test)
#@show test[3]
#@show siteind(test, 3)
xdatasites, ydatadites = energyforbetalisttest(betamax, step, test, δτ, h, s, cutoff, "XY")
xdataMPO, ydataMPO = MBL.energyforbetalist(betamax, step, test, δτ, h, s, cutoff, "XY")
exacterngy1 = [MBL.exactenergyXY1(β, h, γ) for β in xdataMPO]
exacterngy2 = [MBL.exactenergyXY2(beta) for beta in xdataMPO]

gr()
scatter(xdatasites, ydatadites, label="mesure par site", xlabel="β", ylabel="energy", title="model XY, N=$N, cutoff=$cutoff, δτ=$δτ")
scatter!(xdataMPO, ydataMPO, label="mesure MPO")
plot!(xdataMPO, exacterngy1, label="exact energy ChatGPT")
plot!(xdataMPO, exacterngy2, label="S.Katsura exact energy")