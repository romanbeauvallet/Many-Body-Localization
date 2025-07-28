#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using MKL
using ProgressMeter
using JSON
using Random
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
seedlist = input_data["seeds"]
random_draw =  input_data["nombre de tirage"]

# ====================================== Dict

betalist = collect(0:step:betamax)

metadata = Dict{String,Any}(
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
    "nombre de tirage" => random_draw
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
st, dp = MBL.section_trunc(N, gammescale)

###########
Energyarray = Array{Float64}(undef, dp - st + 1, length(betalist), 5, length(seedlist))
Magnetarray = Array{Float64}(undef, dp - st + 1, length(betalist), 5, length(seedlist))
Bonddimlist = Array{Float64}(undef, length(betalist), 5, length(seedlist))
disorderarray = Array{Float64}(undef, N - 1, 5, length(seedlist))
##########

output_data = merge(metadata, results)
open(savefile, "w") do io
    JSON.print(io, output_data, 2)
end
println("\nResults saved in $savefile")
return flush(stdout)

# ====================================== run

@showprogress desc="run over seed" for p in eachindex(seedlist)
    println("seed = ", seedlist[p])
    rgn = MersenneTwister(seedlist[p])
    @showprogress desc="run over tirage" for i in 1:random_draw
        println("tirage numéro : ", i)
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

        # ======================================= save data

        output_data = merge(metadata, results)
        open(savefile, "w") do io
                JSON.print(io, output_data, 2)
            end
        println("\nResults saved in $savefile")
        return flush(stdout)
    end
end

println("Simulation terminée")

