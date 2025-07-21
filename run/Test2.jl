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

N = 20
J = 1
h = 10
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
dmax = 300
betamax = 1
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 20
j = "z"
γ = 0.0
betalist = collect(0:stepbeta:betamax)
seed1 = 1234432
seed2 = 457823451
# ============================== DATA
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

ancilla, s = MBL.AncillaMPO(N)
println("init done")
copy = deepcopy(ancilla)
gatemeasure1, gateevolve1 = MBL.evolutionwithrandomdisordergates(seed1, ancilla, s, h, δτ)
gatemeasure2, gateevolve2 = MBL.evolutionwithrandomdisordergates(seed2, ancilla, s, h, δτ)
println("gates done")
update1 = MBL.trackonesite(betalist, ancilla, δτ, h, s, cutoff, site_measure, seed1)
update2 = MBL.trackonesite(betalist, copy, δτ, h, s, cutoff, site_measure, seed2)
println("update done")

gr()
p = plot()
scatter!(p, betalist, update1, label="seed $seed1", xlabel="β", ylabel="energy of the site $site_measure")
scatter!(p, betalist, update2, label = "seed $seed2")
display(p)