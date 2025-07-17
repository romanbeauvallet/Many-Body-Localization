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
betamax = 10
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"
γ = 0.0
init = 1234
betalist1 = collect(0:0.1:betamax)
betalist2 = collect(0:0.5:betamax)
betalist3 = collect(0:1:betamax)

filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

# ============================== DATA
ancilla, s = MBL.AncillaMPO(N)

xdata1, ydata1 = MBL.energyforbetalist(betalist1, ancilla, δτ, h, s, cutoff, "SS", gammescale)
xdata2, ydata2 = MBL.energyforbetalist(betalist2, ancilla, δτ, h, s, cutoff, "SS", gammescale)
xdata3, ydata3 = MBL.energyforbetalist(betalist3, ancilla, δτ, h, s, cutoff, "SS", gammescale)

#exactenergy = [MBL.exactenergyXY(β, h, γ) for β in xdataMPOstep]


#exactdz = [energyexact(input["spectrum"], beta, N) for beta in xdataMPO]
gr()
p = plot()
scatter!(p, xdata1, ydata1, xlabel="β", ylabel="energy moyenne par site", label ="TEBD step=0.1", title="N=$N, h=$h, cutoff=$cutoff, δτ=$δτ, model SS")
scatter!(p, xdata2, ydata2, label ="TEBD step=0.5")
scatter!(p, xdata3, ydata3, label ="TEBD step=1")
hline!(p, [1/4-log(2)], label="exact energy at 0K without disorder")