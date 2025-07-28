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

N = 50
J = 1
h = 2.5
δτ = 1e-3
dmax = 400
gammescale = 0.7
cutoff = 1e-15
j = "z"
betamax = 5
step = 0.1
seedlist = [34632244, 8789, 876688, 6679092, 12234]
random_draw = 5

# ====================================== Dict

betalist = collect(0:step:betamax)

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
    "random draw number" => random_draw
)
println("\nmetadata:")
display(metadata)

results = Dict{String, Any}(
    "energy [site, beta]" => nothing,
    "magnet [site, beta]" => nothing,
    "disorder list" => nothing,
    "maximum bond dim at each beta" => nothing,
)

# ====================================== begin

st, dp = MBL.section_trunc(N, gammescale)

ancilla, s = MBL.AncillaMPO(N)
Energyarray = Array{Float64}(undef, dp - st + 1, length(betalist), 5, length(seedlist))
Magnetarray = Array{Float64}(undef, dp - st + 1, length(betalist), 5, length(seedlist))
Bonddimlist = Array{Float64}(undef, length(betalist), 5, length(seedlist))
disorderarray = Array{Float64}(undef, N - 1, 5, length(seedlist))
@showprogress desc = "run over seed" for p in eachindex(seedlist)
    println("seed =", seedlist[p])
    rgn = MersenneTwister(seedlist[p])
    for i in 1:random_draw
        println("tirage numéro : ", p)
        e, m, disorder, b = MBL.magnetandenergyforbetalistdisorder(
            betalist, ancilla, δτ, h, s, cutoff, gammescale, dmax, j, rgn
        )
        Energyarray[:, :, i, p] = e
        Magnetarray[:, :, i, p] = m
        Bonddimlist[:, i, p] = b
        disorderarray[:, i, p] = disorder
        results["energy [site, beta]"] = Energyarray
        results["magnet [site, beta]"] = Magnetarray
        results["disorder list"] = disorderarray
        results["maximum bond dim at each beta"] = Bonddimlist
        output_data = merge(metadata, results)
        savefile = joinpath("analyse_simulations_julia", "DATA_Local", "tryrandomdisorder.json")
        open(savefile, "w") do io
                JSON.print(io, output_data, 2)
            end
    end
end

println("Simulation terminée")
