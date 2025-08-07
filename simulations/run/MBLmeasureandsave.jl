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
using PProf
using Profile
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
h = 4.5
δτ = 1e-3
dmax = 300
gammescale = 0.8 #missing in the output data (found in the input data on the cluster)
cutoff = 1e-15
j = "z"
init = 123578574309
betamax = 1
step = 0.1
random_draw = 1
savefile = joinpath("..","analyse_simulations_julia","DATA_Local","debugsimu","outputh0_3.json")
savemps = joinpath("..","analyse_simulations_julia","DATA_Local","debugsimu","mpsh4_5.h5")

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
function void()
    f_h5 = h5open(savemps, "cw")

    for l in 1:random_draw
        update, s = MBL.AncillaMPO(N)
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
end

@pprof void()


println("Simulation terminée")
