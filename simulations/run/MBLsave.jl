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
using ITensorMPS
# ===================== log
println(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
println("Julia $VERSION")
@show Base.julia_cmd()
@show Threads.nthreads()
@show LinearAlgebra.BLAS.get_num_threads()
Pkg.status()
# ===================== parameters
#=
if length(ARGS) < 1
    println("Missing input file: use default")
    throw(ErrorException)
else
    json_input = String(ARGS[1])
end

println("\nLoad input parameters from file $json_input")

input_data = JSON.parsefile(json_input)
map(k -> println(k, ": ", input_data[k]), sort(collect(keys(input_data))))
=#
# ===================================== parameters

N = 20
J = 1
h = 2.78
δτ = 1e-3
dmax = 300
gammescale = 0.8 #missing in the output data (found in the input data on the cluster)
cutoff = 1e-15
j = "z"
init = 123578574309
random_draw = 5
betamax = 5
step = 0.1
random_draw = 5

# ====================================== Dict

betalist = collect(0:step:betamax)
realbetalist = pushfirst!(diff(betalist), 0)

savefile = joinpath("..", "analyse_simulations_julia", "DATA_Local", "testsavemps.json")
savemps = joinpath("..", "analyse_simulations_julia", "DATA_Cluster", "Disorder_parallel", "test", "mps_h_3.247_beta10.h5")

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

results = Dict{String, Any}("disorder list" => nothing, "maximum bond dim at each beta" => nothing)

# ====================================== begin

###########
Disorderlist = fill(1.0, N - 1, random_draw)
##########

results["disorder list"] = Disorderlist

# ====================================== init
function void()
    rng = MersenneTwister(init)
    # ====================================== run
    # Open the HDF5 file ONCE in append mode before the main loop
    # This ensures all MPOs are written to the same file without overwriting.
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
            update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)

            # Write to the already opened HDF5 file (f_h5)
            # The group name will be unique for each beta and draw combination
            write(f_h5, "MPO beta = $(betalist[i]), draw = $l", update)

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
end

# Now, you can read from the file
f = h5open(savemps, "r")
psi = read(f, "MPO beta = 1.4, draw = 2", MPO)
close(f)

println("Reading test worked")
