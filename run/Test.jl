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

N = 20
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

################# try to add quantum number conservation
test, s = MBL.AncillaMPO(N)
beta, energy = MBL.energyforbetalist(betamax, step, test, δτ, h, s, cutoff, "SS")

gr()
scatter(beta, energy, label="<H>", title="N=$N, cutoff=$cutoff, δτ=$δτ, step=$step")