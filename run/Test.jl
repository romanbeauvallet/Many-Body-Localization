#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL

################ Parameters ###############

N = 20
J = 1
h = 0.1
δτ = 1e-2
D = 10
site_measure = div(N, 2)
n_sweep = 300
cutoff = 1e-20
Dmax = 300
Beta = n_sweep*δτ

################ Run basic ################
mps, s = random_initialized_MPS(N, D)
@show typeof(mps)
update = tebdstepHeisenberg!(n_sweep, mps, h, δτ, cutoff, Dmax)
@show update
H = hamiltonianHeisenberg(update, h)
@show typeof(H)
e = measure_H(update, site_measure, H)
@show e 
############## Run convergence #################

time_list = reverse(collect(1e-4:1e-3:1.1e-2))
@show time_list
