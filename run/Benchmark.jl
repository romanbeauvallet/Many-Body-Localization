#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using ProgressMeter
using Plots

################# Parameters ###############
N = 100
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 400
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ

################# Scaling ################
mps_random_debut, _ = random_initialized_MPS(N, D0)
gammelength = (div(N, 10), N)
gammescale = 0.5
j = "z"
gammesweep = (100, 500, 50) #(start, stop, step)

xdata, ydata = MBL.magnetaverageagainstsweep(j, mps_random_debut, gammesweep, gammescale, h, δτ, cutoff, Dmax)

gr()

scatter(xdata, ydata, xlabel="tebd sweep number", ylabel="average magnet (Sz) on $gammescale of $N", title="init with random mps")
