#!/usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Libraries #################
using MBL
import LinearAlgebra: BLAS
@show BLAS.vendor()
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

N = 20
J = 1
h = 0.0
δτ = 1e-3
dmax = 300
gammescale = 0.8 #missing in the output data (found in the input data on the cluster)
cutoff = 1e-15
j = "z"
init = 123578574309
random_draw = 5

# ====================================== Dict

#betalist = input_data["beta list values"]
betalist = collect(0:0.1:5)

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
    "seed" => init,
    "beta list values" => betalist,
    "random draw number" => random_draw
)
println("\nmetadata:")
display(metadata)

results = Dict{String, Any}(
    "energy [site, beta]" => nothing,
    "magnet [site, beta]" => nothing,
    "maximum bond dim at each beta" => nothing
)

# ====================================== init
savefile = joinpath("analyse_simulations_julia", "DATA_Local", "verifcodeh0.json")
st, dp = MBL.section_trunc(N, gammescale)

Energyarray = fill(1.0, dp - st + 1, length(betalist), random_draw)
Magnetarray = fill(1.0, dp - st + 1, length(betalist), random_draw)
Bonddimlist = fill(1.0, length(betalist), random_draw)
Disorderlist = fill(1.0, N-1, random_draw)

results["energy [site, beta]"] = Energyarray
results["magnet [site, beta]"] = Magnetarray
results["maximum bond dim at each beta"] = Bonddimlist
results["disorderlist"] = Disorderlist

rng = MersenneTwister(init)
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
    results["disorderlist"] = Disorderlist
    output_data = merge(metadata, results)
    open(savefile, "w") do io
        JSON.print(io, output_data, 2)
    end
    println("\nResults saved in $savefile")
    flush(stdout)
end
