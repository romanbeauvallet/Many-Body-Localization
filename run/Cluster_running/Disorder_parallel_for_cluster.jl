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

# ===================================== parameters

N = input_data["N"]
J = input_data["J"]
h = input_data["fixed disorder"]
δτ = input_data["Trotter-Suzuki step"]
dmax = input_data["max bond dimension"]
gammescale = input_data["gammescale"]
cutoff = input_data["cutoff"]
j = input_data["axis"]
betamax = input_data["maximum value for beta"]
step = input_data["step beta"]
savefile = input_data["savefile"]
seed = input_data["job_id"]

# ====================================== Dict

betalist = collect(0:step:betamax)

metadata = Dict{String,Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "Given maximal bond dimension" => Dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "disorder" => h,
    "number of spins measured" => gammescale,
    "maximum bond dimension per tebd step" => nothing,
    "seed" => seed,
    "beta list values" => betalist
)
println("\nmetadata:")
display(metadata)

results = Dict{String,Any}(
    "energy [site, beta]" => nothing,
    "magnet [site, beta]" => nothing,
    "disorder list" => nothing,
    "maximum bond dim at each beta" => nothing
)

# ====================================== begin

ancilla, s = MBL.AncillaMPO(N)

Energyarray, Magnetarray, disorder, Bonddimlist = MBL.magnetandenergyforbetalistdisorder(betalist, ancilla, δτ, h, s, cutoff, gammescale, dmax, j, seed)

# ====================================== data

results["energy [site, beta]"] = Energyarray
results["magnet [site, beta]"] = Magnetarray
results["disorder list"] = disorder
results["maximum bond dim at each beta"] = Bonddimlist

output_data = merge(metadata, results)

open(savefile, "w") do io
        JSON.print(io, output_data, 2)
    end
println("\nResults saved in $savefile")
return flush(stdout)

println("Simulation terminée")

