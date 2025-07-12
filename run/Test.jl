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

N = 50
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 1000
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ
gammescale = 1

################# try to add quantum number conservation
mpsrandom, _ = random_initialized_MPS(N, D0)
mps, _ = neelstate(N)

update =tebdstepHeisenbergRow!(3000, mps, 0, 1e-3, 1e-15, 200)
@show energyagainstsite(update, 0, 0.5)

s = ITensors.siteinds("S=1/2", N; conserve_qns=true)
rho = MPO(s, "Id") ./ √2
@show rho