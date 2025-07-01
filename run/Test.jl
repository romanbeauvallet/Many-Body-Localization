#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############


using MBL

################ Parameters ###############

N = 20
J = 1
h = 0.1

################ Run ################

H = heisenberghamiltonian(J, h, N)