#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots

using Base.Threads
println("Nombre de threads disponibles : ", nthreads())

################ Parameters ###############

N = 10
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 1000
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ
gammescale=1

################# try to add quantum number conservation

mps = neelstate(N)

@show mps
@show typeof(mps)
