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
h = 10
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
Dmax = 300
betamax = 5
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"
γ = 0.0
init = 1234
betalist = collect(0:stepbeta:betamax)
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

# ============================== DATA
test, s = MBL.AncillaMPO(N)

xdata, ydata = MBL.energyforbestalistdisorder(betamax, stepbeta, test, δτ, h, s, cutoff, gammescale, init)
#exactenergy = [MBL.exactenergyXY(β, h, γ) for β in xdataMPOstep]


#exactdz = [energyexact(input["spectrum"], beta, N) for beta in xdataMPO]
gr()
p = plot()
scatter!(p, xdata, ydata, xlabel="β", ylabel="energy moyenne par site", label ="TEBD with random disorder", title="N=$N, h=$h, cutoff=$cutoff, model SS")
hline!(p, [1/4-log(2)], label="exact energy at 0K without disorder")