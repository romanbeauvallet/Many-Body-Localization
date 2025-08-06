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
using HDF5
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
init = input_data["seeds"]
random_draw = input_data["nombre de tirage"]
savemps = input_data["h5 file"]

# ====================================== Dict

betalist = collect(0:step:betamax)
realbetalist = pushfirst!(diff(betalist), 0)

metadata = Dict{String, Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "Given maximal bond dimension" => dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "gammescale" => gammescale,
    "disorder" => h,
    "number of spins measured" => gammescale,
    "seed" => init,
    "beta list values" => betalist,
    "nombre de tirage" => random_draw,
    "fichier de sauvegarde MPS" => savemps,
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

st, dp = MBL.section_trunc(N, gammescale)

###########
Energyarray = fill(1.0, dp - st + 1, length(betalist), random_draw)
Magnetarray = fill(1.0, dp - st + 1, length(betalist), random_draw)
Bonddimlist = fill(1.0, length(betalist), random_draw)
Disorderlist = fill(1.0, N-1, random_draw)
##########

results["energy [site, beta]"] = Energyarray
results["magnet [site, beta]"] = Magnetarray
results["maximum bond dim at each beta"] = Bonddimlist
results["disorder list"] = Disorderlist

# ====================================== init
rng = MersenneTwister(init)
# ====================================== run
f_h5 = h5open(savemps, "cw")

for l in 1:random_draw
    ancilla, s = MBL.AncillaMPO(N)
    update = ancilla
    println("\n# " * "="^90)
    println("random draw number = ", l)
    flush(stdout)
    gatesmeasure, gatesevolve, Disorder = MBL.evolutionwithrandomdisordergates(
        rng, update, s, h, δτ
    )
    Disorderlist[:, l] = Disorder
    results["disorder list"] = Disorderlist
    for i in eachindex(realbetalist)
        println("\n# " * "="^30)
        beta = betalist[i]
        @info "β[$i]" beta
        update = TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)

        # Write to the already opened HDF5 file (f_h5)
        # The group name will be unique for each beta and draw combination
        write(f_h5, "MPO beta = $(betalist[i]), draw = $l", update)
        _, Energyarray[:, i, l] = MBL.energyagainstsiteMPOdisorder(update, gatesmeasure, gammescale)
        results["energy [site, beta, draw]"] = Energyarray
        println("Average energy at β=$beta : ", mean(Energyarray; dims=1)[i])
        #println("\n# " * "="^30)
        _, Magnetarray[:, i, l] = magnetagainstsite(update, j, gammescale)
        results["magnet [site, beta, draw]"] = Magnetarray
        println("Average Sz at β=$beta : ", mean(Magnetarray; dims=1)[i])
        Bonddimlist[i, l] = maxbonddim(update)
        println("Maximum bond dimension at β=$beta : ", Bonddimlist[i, l])
        results["maximum bond dim at each beta"] = Bonddimlist
        flush(stdout)
    end
    output_data = merge(metadata, results)
    open(savefile, "w") do io
        JSON.print(io, output_data, 2)
    end
    println("\nResults saved in $savefile")
    flush(stdout)
end

# Close the HDF5 file AFTER the entire simulation is finished
close(f_h5)

println("Simulation terminée")
