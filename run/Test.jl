#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
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
EnergyList = Vector{}()
mps_init, _ = random_initialized_MPS(N,D)
time_list = reverse(collect(1e-4:1e-3:1.1e-2))
function void()
    @showprogress for k in eachindex(time_list)
        update = tebdstepHeisenberg!(n_sweep, mps_init, h, time_list[k], cutoff, Dmax)
        H = hamiltonianHeisenberg(update, h)
        e = measure_H(update, site_measure, H)
        push!(EnergyList, e)
    end
end
void()
############ Plots #################
gr()

plot(time_list, EnergyList)