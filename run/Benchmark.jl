#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using MKL
using ProgressMeter
using JSON
using Statistics
using Dates
using LinearAlgebra
using Pkg
# ===================== log
println(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
println("Julia $VERSION")
@show Base.julia_cmd()
@show Threads.nthreads()
@show LinearAlgebra.BLAS.get_num_threads()
Pkg.status()
# ===================== parameters
N = 20
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 100
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ

################# Scaling ################
gammelength = (div(N, 10), N)
gammescale = 0.5
j = "z"
gammesweep = (1000, 3000, 500) #(start, stop, step)

################# Run ###############
sweep_list = collect(gammesweep[1]:gammesweep[3]:gammesweep[2])
realsweeplist = [gammesweep[3] for k in 1:floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3]))+1]
realsweeplist[1] = gammesweep[1]
mps_random_debut, _ = random_initialized_MPS(N, D0)

# ===================== data 

metadata = Dict{String,Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "initial bond dimension" => D0,
    "Dmax" => Dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "sweep range" => sweep_list,
    "effective sweep list" => realsweeplist,
    "maximum bond dimension per tebd step" => nothing
)
println("\nmetadata:")
display(metadata)

results = Dict{String,Any}(
    "energy sweep list" => nothing,
    "magnetization sweep list" => nothing,
)

Energytebd = Vector()
Magnettebd = Vector()
Maxbonddim = Vector()
#####evolv
function void()
    update_tebd = deepcopy(mps_random_debut)
    for i in eachindex(realsweeplist)
        println("Time evolution with tebd")
        update_tebd = tebdstepHeisenbergRow!(i, update_tebd, h, δτ, cutoff, Dmax)
        push!(Maxbonddim,maxbonddim(update_tebd))
        metadata["maximum bond dimension per tebd step"] = Maxbonddim

        #####measure
        println("Measure average energy")
        _, e = energyagainstsite(update_tebd, h, gammescale)
        push!(Energytebd, mean(e))
        results["energy sweep list"] = Energytebd

        println("Mesure average magnet")
        _, magnet = magnetagainstsite(update_tebd, j, gammescale)
        push!(Magnettebd, mean(magnet))
        results["magnetization sweep list"] = Magnettebd

        #####data saving
        output_data = merge(metadata, results)
        #savefile = get_savefile(output_data)
        open(savefile, "w") do io
            JSON.print(io, output_data, 4)
        end
        println("\nResults saved in $savefile")
        flush(stdout)
    end
end

void()
println("simulation finie")
