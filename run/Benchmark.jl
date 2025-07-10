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
if length(ARGS) < 1
    println("Missing input file: use default")
    throw(ErrorException)
else
    json_input = String(ARGS[1])
end

println("\nLoad input parameters from file $json_input")

input_data = JSON.parsefile(json_input)
map(k -> println(k, ": ", input_data[k]), sort(collect(keys(input_data))))

N = input_data["N"]
J = input_data["J"]
D0 = input_data["D0"]
h = input_data["disorder"]
δτ = input_data["Trotter-Suzuki step"]
Dmax = input_data["max bond dimension"]
gammesweep = input_data["nsweep range"]
gammescale = input_data["gammescale"]
cutoff = input_data["cutoff"]
n_sweep = input_data["fixed number of sweep"]
j = input_data["axis"]
savefile = String(input_data["savefile"])

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
        metadata["maximum bond dimension"] = Maxbonddim

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
