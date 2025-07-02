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

N = 40
J = 1
h = 0
δτ = 5e-3
D = 10
site_measure = div(N, 2)
n_sweep = 200
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
EnergytotList = Vector{}()
EnergysiteList = Vector{}()
Sz_list = Vector{}()
mps_init, _ = random_initialized_MPS(N, D)
time_list = reverse(collect(1e-5:1e-3:2e-2))
site_list = collect(1:1:N)
function voidenergy()
    update = tebdstepHeisenbergRow!(n_sweep, mps_init, h, 2e-2, cutoff, Dmax)
    @showprogress for k in eachindex(time_list)
        update = tebdstepHeisenbergRow!(n_sweep, update, h, time_list[k], cutoff, Dmax)
        #H = hamiltonianHeisenberg(update, h)
        e = energysite(update, site_measure, h)
        push!(EnergyList, e)
    end
end

voidenergy()

function voidmagnet()
    update = tebdstepHeisenbergRow!(n_sweep, mps_init, h, δτ, cutoff, Dmax)
    @showprogress for k in eachindex(site_list)
        m = measure_Sz(update, site_list[k])
        push!(Sz_list, m)
    end
end
#voidmagnet()

function voidenergysite()
    update = tebdstepHeisenbergRow!(n_sweep, mps_init, h, δτ, cutoff, Dmax)
    @showprogress for k in 1:2:N-2
        e = energysite(update, k, h)
        push!(EnergysiteList, e)
    end
end
#voidenergysite()
############ Plots #################
gr()

scatter(time_list, EnergyList, label="TEBD", xlabel="\$δτ\$", ylabel="energy in the middle", title="n_sweep = 100, N=30, cutoff=1e-20, h=0")
#scatter(site_list, Sz_list, label="TEBD", xlabel="site", ylabel="\$S_z\$", title="n_sweep = 100, N=100, cutoff = 1e-20, h=0")
#scatter(site_list, EnergyList, label="TEBD", xlabel="site", ylabel="\$ϵ\$", title="n_sweep = 100, N=100, cutoff = 1e-20, h=0")