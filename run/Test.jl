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
step = 0.5
Beta = n_sweep * δτ
gammescale = 0.8
j = "z"
γ = 0.0
betalist = collect(0:step:betamax)
filename = joinpath("analyse_simulations_julia", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

# ============================== DATA
test, s = MBL.AncillaMPO(N)

xdatasites, ydatadites = MBL.energyforbetalist(betamax, step, test, δτ, h, s, cutoff, "SS", gammescale)
xdataMPO, ydataMPO = MBL.energyforbetalistMPO(betamax, step, test, δτ, h, s, cutoff, "SS")

#exactenergy = [MBL.exactenergyXY(β, h, γ) for β in xdataMPO]

#exactdz = [energyexact(input["spectrum"], beta, N) for beta in xdataMPO]

gr()
p = plot()
scatter!(p, xdatasites, ydatadites, label="mesure par site")
scatter!(p, xdataMPO, ydataMPO / N, label="mesure MPO")
#plot!(xdataMPO, exactdz, label="exact dz", xlabel="β", ylabel="<H>/N", title="N=$N, δτ=$δτ, cutoff=$cutoff, model XY")
#plot!(p, xdataMPO, exactenergy, label=latexstring("-\\frac{1}{4π}\\int_{-π}^{π} cos(k)tanh(βcos(k)/2)dk"), xlabel="β", ylabel="<H>/N", title="N=$N, δτ=$δτ, cutoff=$cutoff, model XY")
#hline!([1/4-log(2)], label="exact energy at zero temperature")
display(p)