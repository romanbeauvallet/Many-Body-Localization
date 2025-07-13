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
init = input_data["initialization"]
savefile = String(input_data["savefile"])

################# Run ###############
mpsdebutrandom = random_initialized_MPS(N, D0)
mpsneel = neelstate(N)
# ===================== data 

metadata = Dict{String,Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "initial bond dimension" => D0,
    "Dmax" => Dmax,
    "J" => J,
    "cutoff" => cutoff,
    "disorder" => h,
    "maximum bond dimension per tebd step" => nothing,
)
println("\nmetadata:")
display(metadata)

results = Dict{String,Any}(
    "correlationN x" => nothing,
    "correlationN y" => nothing,
    "correlationN z" => nothing,
    "correlationR x" => nothing,
    "correlationR y" => nothing,
    "correlationR z" => nothing
)

updateN = tebdstepHeisenbergRow!(n_sweep, mpsneel, h, δτ, cutoff, Dmax)
updateR = tebdstepHeisenbergRow!(n_sweep, mpsdebutrandom, h, δτ, cutoff, Dmax)

results["correlationN x"] = correlationagainstsite(updateN, "x")
results["correlationN y"] = correlationagainstsite(updateN, "y")
results["correlationN z"] = correlationagainstsite(updateN, "z")
results["correlationR x"] = correlationagainstsite(updateR, "x")
results["correlationR y"] = correlationagainstsite(updateR, "y")
results["correlationR z"] = correlationagainstsite(updateR, "z")


output_data = merge(metadata, results)
open(savefile, "w") do io
    JSON.print(io, output_data, 4)
end
println("\nResults saved in $savefile")
flush(stdout)

