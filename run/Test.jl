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

N = 100
J = 1
h = 0
δτ = 1e-3
D = 10
site_measure = div(N, 2)
n_sweep = 400
cutoff = 1e-30
Dmax = 300
Beta = n_sweep * δτ

################ Run basic ################
function test()
    mps, s = random_initialized_MPS(N, D)
    @show typeof(mps)
    update_test = tebdstepHeisenberg!(n_sweep, mps, h, δτ, cutoff, Dmax)
    @show update_test
    H_test = hamiltonianHeisenberg(update_test, h)
    @show typeof(H_test)
    e_test = measure_H(update_test, site_measure, H_test)
    @show e_test
end
############## Run convergence #################

mps_start, s = random_initialized_MPS(N, D)
copy = deepcopy(mps_start)
converged_mps_2step = tebdstepHeisenberg!(n_sweep, mps_start, h, δτ, cutoff, Dmax)
converged_mps = tebdstepHeisenbergRow!(n_sweep, copy, h, δτ, cutoff, Dmax)

sites_2steps, magnetlist_2step = magnetagainstsite(converged_mps_2step)
sites_1step, magnetlist_1step = magnetagainstsite(converged_mps)

############## Graphs ##############

gr()

plot1 = scatter(sites_1step, magnetlist_1step, label="tebd row", xlabel="site", ylabel="Sz", title="N = $N, nsweep = $n_sweep, cutoff = $cutoff")
scatter!(sites_2steps, magnetlist_2step, label="ordered tebd")
display(plot1)



