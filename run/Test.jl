#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots
using ITensorMPS

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
betamax = 10
step = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"
γ=0.0

################# try to add quantum number conservation
MBL.exactenergyXY(1, 0, 0)
test, s = MBL.AncillaMPO(N)
beta, energy = MBL.energyforbetalist(betamax, step, test, δτ, h, s, cutoff, "XY")
exactenergy = collect([MBL.exactenergyXY(beta[i], h, γ) for i in eachindex(beta)])
gr()
scatter(beta, energy/N, ylabel="<H>", title="N=$N, cutoff=$cutoff, δτ=$δτ, step=$step", xlabel="β", label="tebd")
scatter!(beta, exactenergy, label="exact energy")
