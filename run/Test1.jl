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

N = 50
J = 1
h0 = 0
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 5000
cutoff = 1e-15
dmax = 300
betamaxstep = 10
betamax = 10
stepbeta1 = 0.1
stepbeta2 = 0.5
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 10
beta = 40
j = "z"
γ = 0.0
init = 4562
betalist1 = collect(0:stepbeta1:betamaxstep)
betalist2 = collect(betamaxstep:stepbeta2:betamax)
betalist = collect(Iterators.flatten((betalist1, betalist2)))
# ============================== DATA
ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)
Hamiltonian = MBL.hamiltonianXY(mps, h, smps)
###Energy

#energybetaMPO = MBL.energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, "XY", gammescale)
gates = MBL.gatesTEBDancilla(ancilla, h, δτ, s, "XY")
update = MBL.TEBDancilla!(ancilla, gates, 40, cutoff, δτ)
xdata, energysiteMPO = energyagainstsite(update, h, gammescale, "XY")
@show length(energybetaMPO)
###DMRG
st, dp = MBL.section_trunc(N, gammescale)
L = collect(st:1:dp)
psi, H = MBL.groundstateDMRG(mps, Hamiltonian, n_sweepDMRG, dmax, cutoff, noise)
xdata2, exactpersite = energyagainstsite(psi, h, gammescale, "XY")

###Exact energy

exact = MBL.exactenergyXY(beta, h, γ)
exactDMRG = mean(exactpersite)
exactpersite = mean(energysiteMPO)
###Graphes

gr()
p1 = plot()
scatter!(p1, betalist, energybetaMPO, label="TEBD+purification")
scatter!(p1, xdata2, exactpersite, label="DMRG energy on each site", xlabel="site index", ylabel="Energy of each site", title="N=$N, h=$h, δτ=$δτ, β=$β, Model XY", titlefont=font("Computer Modern", 16), legendfont=font("Computer Modern", 13))
scatter!(p1, [10], [exact],marker=(:star5, 5), label="Exact energy per site")
scatter!(p1, [10], [exactDMRG],marker=(:star5, 5), label="DMRG energy per site")
scatter!(p1, [10], [exactpersite],marker=(:star5, 5), label="TEBD energy per site")
#hline!(p1, [H], label="DMRG MPO/N")
display(p1)
savefig("run/Plots/verifenergypersiteXY.pdf")