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
gammelength = input_data["length range"]
δτ = input_data["Trotter-Suzuki step"]
Dmax = input_data["max bond dimension"]
gammesweep = input_data["nsweep range"]
gammescale = input_data["gammescale"]
cutoff = input_data["cutoff"]
n_sweep = input_data["fixed number of sweep"]
j = input_data["axis"]
init = input_data["initialization"]
initseed = input_data["graine"]
stepbeta = input_data["pas de beta"]
betamax = input_data["limite beta"]
savefile = String(input_data["savefile"])

lengthlist = collect(gammelength[1]:gammelength[3]:gammelength[2])

# ===================== data 

metadata = Dict{String,Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "initial bond dimension" => D0,
    "Dmax" => Dmax,
    "J" => J,
    "axis spin" => j,
    "length list" => lengthlist,
    "cutoff" => cutoff,
    "disorder" => h,
    "proportion spin average" => gammescale,
    "maximum bond dimension per tebd step" => nothing,
    "type d'initialisation" => init,
    "seed for random" => initseed,
    "pas pour beta" => stepbeta,
    "limite max pour beta" => betamax
)
println("\nmetadata:")
display(metadata)

results = Dict{String,Any}(
    "energy sweep list no disorder" => nothing,
    "magnetization sweep list no disorder" => nothing,
    "energy sweep list random disorder" => nothing,
    "magnetization sweep list random disorder" => nothing,
)

ancilla, s = MBL.AncillaMPO(N)
betalist = collect(0:stepbeta:betamax)
realbetalist = pushfirst!(diff(betalist), 0)

Energylist = Vector()
AverageMagnetlist = Vector()
Magnetlist = Vector()

function voiddisorder()
    update = ancilla
    gatesmeasure, gatesevolve = evolutionwithrandomdisordergates(initseed, update, s, h, δτ)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ)
        _, E = energyagainstsiteMPOdisorder(update, gatesmeasure, gammescale)
        push!(Energylist, E)
        results["energy sweep list random disorder"] = Energylist
        _, M =
            push!(Magnetlist, M)

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

function void()

end

voiddisorder()
void()
println("simulation finie")
