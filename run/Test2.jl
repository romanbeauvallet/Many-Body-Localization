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
betamax = 15
stepbeta = 1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 50
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
copy = deepcopy(ancilla)

H, gates = MBL.evolutionwithrandomdisordergates(init, ancilla, s, h, δτ)

update = MBL.TEBDancilla(ancilla, gates, 10, cutoff, δτ)
psi, H = MBL.groundstateDMRG(copy, H, n_sweepDMRG, dmax, cutoff, noise)

L = MBL.section_trunc(N, gammescale)
xdataDMRG = Vector()
xdata = Vector()

for i in 1:length(L)-1
    T = randn(2,2,2,2)
    s1, s2 = siteind(update, L[i]), siteind(update, L[i+1])
    p1, p2 = siteind(psi, L[i]), siteind(psi, L[i+1])
    S = MBL.randomoperator(T,s1,s2)
    P = MBL.randomoperator(T, p1,p2)
    push!(xdata, energysiteMPOdisorder(update, L[i], S))
    push!(xdataDMRG, energysiteMPOdisorder(psi, L[i], P))
end
