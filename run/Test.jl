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
D0 = 10
site_measure = div(N, 2)
n_sweep = 100
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ
gammescale=1
################ Run basic ################
function test()
    mps, s = random_initialized_MPS(N, D0)
    @show typeof(mps)
    update_test = tebdstepHeisenberg!(n_sweep, mps, h, δτ, cutoff, Dmax)
    @show update_test
    H_test = hamiltonianHeisenberg(update_test, h)
    @show typeof(H_test)
    e_test = measure_H(update_test, site_measure, H_test)
    @show e_test
end
############## Run convergence #################
mps_neel_start, s_neel = neelstate(N)
mps_random_start, s_random = random_initialized_MPS(N, D0)

converged_mps_random = tebdstepHeisenbergRow!(n_sweep, mps_random_start, h, δτ, cutoff, Dmax)
converged_mps_neel = tebdstepHeisenbergRow!(n_sweep, mps_neel_start, h, δτ, cutoff, Dmax)

sites_neel_magnet, magnetlist_neel = magnetagainstsite(converged_mps_neel, "z", gammescale)
sites_random_magnet, magnetlist_random = magnetagainstsite(converged_mps_random, "z", gammescale)

site_neel_energy, energy_neel = energyagainstsite(converged_mps_neel, h)
site_random_energy, energy_random = energyagainstsite(converged_mps_random, h)

#@show mps_start
#@show converged_mps
#@show converged_mps_2step
############## Graphs ##############

#les deux tebd sont bien les mêmes, la version 2 steps est légèrement plus lente

gr()

plotneel1 = scatter(sites_neel_magnet, magnetlist_neel, label="tebd", xlabel="site", ylabel="Sz", title="N = $N, nsweep = $n_sweep, cutoff = $cutoff, neel")
plotneel2 = scatter(site_neel_energy, energy_neel, label="tebd", xlabel="site", ylabel="ϵ", title="N = $N, nsweep = $n_sweep, cutoff = $cutoff, neel")

plotrandom2 = scatter(site_neel_energy, energy_neel, label="tebd", xlabel="site", ylabel="ϵ", title="random")
plotrandom1 = scatter(site_random_energy, energy_random, label="tebd", xlabel="site", ylabel="Sz", title="random")

scatter(plotneel1, plotneel2, plotrandom1, plotrandom2, layout = (2, 2), legend = :topleft)