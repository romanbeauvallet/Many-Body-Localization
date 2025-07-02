#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############

using MBL

################ Parameters ###############

N = 20
J = 1
h = 0.1
δτ = 1e-3
D = 10
site_measure = div(N, 2)

################ Run ################

gateseven, H, s = gateTrotterSuzukiandhamiltonian(N, h, δτ, "even")
gatesodd, _, _ = gateTrotterSuzukiandhamiltonian(N, h, δτ, "odd")
@show typeof(H)
mps = random_initialized_MPS(s, D)
e = measure_H(mps, site_measure, H)
