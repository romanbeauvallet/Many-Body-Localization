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
h0 = 0
h = 3.5
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 1000
cutoff = 1e-15
dmax = 300
betamax = 10
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 10
j = "z"
γ = 0.0
init = 1234
betalist = collect(0:stepbeta:betamax)

# ============================== DATA
ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)
H = hamiltonianHeisenberg(mps, h, smps)

xdata1, ydata1 = MBL.energyforbestalistdisorder(betalist, ancilla, δτ, h, s, cutoff, gammescale, init)
#xdata2, ydata2 = MBL.energyforbetalist(betalist, ancilla, δτ, h0, s, cutoff, "SS", gammescale)
#exactenergy = [MBL.exactenergyXY(β, h, γ) for β in xdataMPOstep]

gr()
p = plot()
scatter!(p, xdata1, ydata1, xlabel="β", ylabel="energy moyenne par site", label ="TEBD step=0.5 (h)", title="N=$N, h=$h, cutoff=$cutoff, δτ=$δτ, model SS")
#scatter!(p, xdata2, ydata2, label ="TEBD step=0.5")
#hline!(p, [1/4-log(2)], label="exact energy at 0K without disorder")

