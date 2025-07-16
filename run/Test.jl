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
step0 = 1
step1 = 0.1
step2 = 0.01
step3 = 0.001
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"
γ = 0.0
betalist = collect(0:step:betamax)
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

# ============================== DATA
test, s = MBL.AncillaMPO(N)

xdataMPOstep1, ydataMPOstep1 = MBL.energyforbetalistMPO(betamax, step1, test, δτ, h, s, cutoff, "XY")
xdataMPOstep2, ydataMPOstep2 = MBL.energyforbetalistMPO(betamax, step2, test, δτ, h, s, cutoff, "XY")
xdataMPOstep3, ydataMPOstep3 = MBL.energyforbetalistMPO(betamax, step3, test, δτ, h, s, cutoff, "XY")
xdataMPOstep4, ydataMPOstep4 = MBL.energyforbetalistMPO(betamax, step0, test, δτ, h, s, cutoff, "XY")


#exactenergy1 = [MBL.exactenergyXY(β, h, γ) for β in xdataMPOstep1]
#exactenergy2 = [MBL.exactenergyXY(β, h, γ) for β in xdataMPOstep2]


#exactdz = [energyexact(input["spectrum"], beta, N) for beta in xdataMPO]

xdatasites, ydatasites = MBL.energyforbetalist(betamax, step1, test, δτ, h, s, cutoff, "XY", gammescale)
xdatasites2, ydatasites2 = MBL.energyforbetalist(betamax, step2, test, δτ, h, s, cutoff, "XY", gammescale)
xdatasites3, ydatasites3 = MBL.energyforbetalist(betamax, step3, test, δτ, h, s, cutoff, "XY", gammescale)
xdatasites4, ydatasites4 = MBL.energyforbetalist(betamax, step0, test, δτ, h, s, cutoff, "XY", gammescale)

gr()
p = plot()
scatter!(p, xdatasites, ydatasites, label="mesure par site step=0.1", xlabel="β", ylabel="<H>/N", title="N=$N, δτ=$δτ, cutoff=$cutoff, model XY")
scatter!(p, xdatasites2, ydatasites2, label="mesure par site step=0.01")
scatter!(p, xdatasites3, ydatasites3, label="mesure par site step=0.001")
scatter!(p, xdatasites4, ydatasites4, label="mesure par site step=1")
plot!(p, xdataMPOstep1, ydataMPOstep1, label="mpo step=0.1")
plot!(p, xdataMPOstep2, ydataMPOstep2, label="mpo step=0.01")
plot!(p, xdataMPOstep3, ydataMPOstep3, label="mpo step=0.001")
plot!(p, xdataMPOstep4, ydataMPOstep4, label="mpo step=1")
display(p)