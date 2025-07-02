#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############

using MBL

################ Parameters ###############

N = 20
J = 1
h = 0.1
δτ = 1e-3

################ Run ################

H = gateTrotterSuzukiandhamiltonian(N, h, δτ, "even")
@show H[1]
