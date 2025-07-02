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

N = 20
J = 1
h = 0
δτ = 1e-2
D = 10
site_measure = div(N, 2)
n_sweep = 300
cutoff = 1e-20
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
EnergyList = zeros(Float64, nthreads())
mps_init, _ = random_initialized_MPS(N, D)
time_list = reverse(collect(1e-5:1e-3:2e-2))
@show length(time_list)
function void()
    update = tebdstepHeisenberg!(n_sweep, mps_init, h, 1e-2, cutoff, Dmax)
    used = falses(nthreads())
    @threads for k in eachindex(time_list)
        tid = threadid()
        println("Thread ", tid, " travaille sur i = ", k)
        used[tid] = true
        update = tebdstepHeisenberg!(n_sweep, update, h, time_list[k], cutoff, Dmax)
        H = hamiltonianHeisenberg(update, h)
        e = measure_H(update, site_measure, H)
        EnergyList[tid] = e
    end
end
void()
############ Plots #################
gr()

scatter(time_list, EnergyList, label="TEBD", xlabel="\$δτ\$", ylabel="energy", title="n_sweep = 300, N=20, cutoff=1e-20, h=0")