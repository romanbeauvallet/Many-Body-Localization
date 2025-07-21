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
h = 100
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
dmax = 300
betamax = 5
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 20
j = "z"
γ = 0.0
betalist = collect(0:stepbeta:betamax)
seed1 = 1234432
seed2 = 4578
# ============================== DATA
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)
gatesmeasure, gatesevolve, disorder = MBL.evolutionwithrandomdisordergates(seed1, ancilla, s, h, δτ)
Hamiltonian = MBL.hamiltonianXY(mps, h, smps)
println("init done")
gates = MBL.gatesTEBDancilla(ancilla, h, δτ, s, "XY")
update = MBL.TEBDancilla!(ancilla, gates, 20, cutoff, δτ)
xdata, energysiteMPO = magnetagainstsite(ancilla, "z", gammescale)
st, dp = MBL.section_trunc(N, gammescale)
L = collect(st:1:dp)

psi, e = MBL.groundstateDMRG(mps, Hamiltonian, n_sweepDMRG, dmax, cutoff, noise)
xdata, exactenergysite = magnetagainstsite(psi, "z", gammescale)

#@show update1[5] == update2[5]

gr()
p = plot()
scatter!(p, [10], [mean(energysiteMPO)], label="TEBD", xlabel="chain site", ylabel="magnetization per site")
scatter!(p, [10], [mean(exactenergysite)], label = "DMRG")
display(p)