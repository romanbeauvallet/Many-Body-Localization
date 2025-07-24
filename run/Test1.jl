#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using Distributed
using Base.Threads
using LaTeXStrings
println("Nombre de threads disponibles : ", nthreads())
using JSON
#@everywhere begin
using ITensors
using MBL
using ProgressMeter
using Plots
using ITensorMPS
using Statistics
#end

# Ajouter des processus seulement si tu es en local
#if nprocs() == 1
#addprocs(4)  # ne sera pas utilisé sur SLURM
#end
################ Parameters ###############
N = 50
J = 1
h0 = 0
h = 0.0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 5000
cutoff = 1e-15
dmax = 300
betamaxstep = 5
betamax = 10
stepbeta1 = 0.1
stepbeta2 = 0.5
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 5
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
Hamiltonian = operator(mps, h, smps, "SS")

###Energy

energybetaMPO = MBL.energyforbetalist(
    betalist, ancilla, δτ, h, s, cutoff, "SS", gammescale, dmax
)
@show size(energybetaMPO)
energy = vec(mean(energybetaMPO; dims=1))
@show typeof(energy)
println("done energy")
specificheat = MBL.specificheat(energy, betalist)
###Graphes

gr()
default(; fontfamily="Computer Modern")
p1 = plot()
scatter!(
    p1,
    betalist,
    energy;
    markershape=:star5,
    label="TEBD+purification",
    xlabel="β",
    ylabel="Average energy of one site",
    title="N=$N, h=$h, δτ=$δτ, Model Heisenberg",
)
scatter!(
    p1,
    betalist,
    specificheat;
    markershape=:star5,
    label="TEBD+purification+finite difference",
)
hline!(
    p1,
    [1 / 4 - log(2)];
    label="Analytical result for periodic boundary in thermodynamic limit",
)
display(p1)
savefig("run/Plots/specificheatSS.pdf")