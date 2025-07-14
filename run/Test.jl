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

N = 10
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
Dmax = 300
beta = 1
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"

################# try to add quantum number conservation
test, s = MBL.AncillaMPO(N)
@show typeof(test)
@show length(test)
gates= MBL.gatesTEBDancilla(test, 0, 1e-2, s, "SS")
@show typeof(gates)
_, Energytry = MBL.EnergyAncilla(test, δτ, h, beta, s, cutoff, "SS")
@show Energytry