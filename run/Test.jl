#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL

################ Parameters ###############

N = 20
J = 1
h = 0.1
δτ = 1e-3
D = 10
site_measure = div(N, 2)
n_sweep = 400
cutoff = 1e-15
Dmax = 100
Beta = n_sweep*δτ

################ Run ################
mps, s = random_initialized_MPS(N, D)
@show typeof(mps)
update = tebdstepHeisenberg!(n_sweep, mps, h, δτ, cutoff, Dmax)
@show update
e = measure_H(update, site_measure, H)
