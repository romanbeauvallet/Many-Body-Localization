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

N = 100
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"

################# try to add quantum number conservation
mpsrandom, _ = random_initialized_MPS(N, D0)
mps, _ = neelstate(N)

update = tebdstepHeisenbergRow!(n_sweep, mpsrandom, h, δτ, cutoff, Dmax)
xdata, ydata = correlationagainstsite(update, j)

gr()
scatter(xdata, abs.(ydata))
