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
using Random
# ===================== log
println(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
println("Julia $VERSION")
@show Base.julia_cmd()
@show Threads.nthreads()
@show LinearAlgebra.BLAS.get_num_threads()
Pkg.status()
# ===================================== parameters

#=
json_input = joinpath(
    "analyse_simulations_julia", "DATA_Cluster", "Disorder_parallel", "output_h_4.973_beta10.json"
)
input_data = JSON.parsefile(json_input)
#map(k -> println(k, ": ", input_data[k]), sort(collect(keys(input_data))))

N = input_data["N"]
J = input_data["J"]
h = input_data["disorder list"]
δτ = input_data["Trotter-Suzuki time step"]
dmax = input_data["Given maximal bond dimension"]
gammescale = 0.8 #missing in the output data (found in the input data on the cluster)
cutoff = input_data["cutoff"]
j = "z"
seedlist = input_data["seed"]
random_draw = input_data["nombre de tirage"]
=#

N = 100
J = 1
h = 0.309
δτ = 1e-3
dmax = 300
gammescale = 0.8 #missing in the output data (found in the input data on the cluster)
cutoff = 1e-15
j = "z"
seedlist = [123578574309]
random_draw = 3

# ====================================== Dict

#betalist = input_data["beta list values"]
betalist = collect(0:0.1:3)

metadata = Dict{String, Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "Given maximal bond dimension" => dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "disorder" => h,
    "number of spins measured" => gammescale,
    "maximum bond dimension per tebd step" => nothing,
    "seed" => seedlist,
    "beta list values" => betalist,
    "random draw number" => random_draw,
)
println("\nmetadata:")
display(metadata)

results = Dict{String, Any}(
    "energy [site, beta]" => nothing,
    "magnet [site, beta]" => nothing,
    "maximum bond dim at each beta" => nothing,
)

# ====================================== init
savefile = joinpath("analyse_simulations_julia", "DATA_Local", "verifcodeh0309.json")
st, dp = MBL.section_trunc(N, gammescale)

Energyarray = fill(1.0, dp - st + 1, length(betalist), random_draw)
Magnetarray = fill(1.0, dp - st + 1, length(betalist), random_draw)
Bonddimlist = fill(1.0, length(betalist), random_draw)
Disorderlist = fill(1.0, N-1, random_draw)
rng = MersenneTwister(seed)
ancilla, s = MBL.AncillaMPO(N)
# ===================================== run

for i in 1:random_draw
    println("\n# " * "="^90)
    println("random draw number = ", i)
    flush(stdout)
    e, m, d, b = MBL.magnetandenergyforbetalistdisorder(
    betalist, ancilla, δτ, h, s, cutoff, gammescale, rng, dmax, j
)
    Energyarray[:, :, i] = e
    Magnetarray[:, :, i] = m
    Bonddimlist[:, i] = b
    Disorderlist[:, i] = d
    results["energy [site, beta]"] = Energyarray
    results["magnet [site, beta]"] = Magnetarray
    results["maximum bond dim at each beta"] = Bonddimlist
    results["disorderlist"] = disorderlist
    output_data = merge(metadata, results)
    open(savefile, "w") do io
        JSON.print(io, output_data, 2)
    end
    println("\nResults saved in $savefile")
    flush(stdout)
end
