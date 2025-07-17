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
dmax = 300
betamax = 20
stepbeta = 0.5
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 10
j = "z"
γ = 0.0
init = 1234
betalist = collect(0:stepbeta:betamax)

# ============================== DATA
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)
H = MBL.hamiltonianXY(mps, h, smps)

_, H = MBL.groundstateDMRG(mps, H, n_sweepDMRG, dmax, cutoff, noise)
xdata2, ydata2 = MBL.energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, "XY", gammescale)
exactenergy = [MBL.exactenergyXY(β, h, γ) for β in xdata2]


#exactdz = [energyexact(input["spectrum"], beta, N) for beta in xdataMPO]
gr()
p = plot()
scatter!(p, xdata2, ydata2, xlabel="β", ylabel="energy moyenne par site", label ="TEBD", title="N=$N, h=$h, cutoff=$cutoff, δτ=$δτ, model XY")
hline!(p, [H], label="exact energy at 0K DMRG MPO")
plot!(p, xdata2, exactenergy, label="exact energy")
