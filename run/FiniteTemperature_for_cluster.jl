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
γ = input_data["gamma"]
D0 = input_data["D0"]
h = input_data["disorder"]
δτ = input_data["Trotter-Suzuki step"]
Dmax = input_data["max bond dimension"]
gammescale = input_data["gammescale"]
cutoff = input_data["cutoff"]
j = input_data["axis"]
noise = input_data["noise"]
sweep_DMRG = input_data["number sweep dmrg"]
beta = input_data["beta fixe"]
savefile = String(input_data["savefile"])

# ===================== data 

metadata = Dict{String,Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "initial bond dimension" => D0,
    "Dmax" => Dmax,
    "J" => J,
    "beta fixe" => beta,
    "axis spin" => j,
    "cutoff" => cutoff,
    "disorder" => h,
    "n_sweepDMRG" => sweep_DMRG,
    "proportion spin average" => gammescale,
    "maximum bond dimension per tebd step" => nothing
)
println("\nmetadata:")
display(metadata)

results = Dict{String,Any}(
    "energy sites tebd" => nothing,
    "energy sites dmrg" => nothing,
    "exact energy" => nothing,
    "exact energy exact" => nothing,
    "exact energy dmrg" => nothing, 
    "liste sites" => nothing
)

Energylist = Vector()
AverageMagnetlist = Vector()
Magnetlist = Vector()

function void()
    ancilla, s = MBL.AncillaMPO(N)
    mps, smps = neelstate(N)
    Hamiltonian = MBL.hamiltonianXY(mps, h, smps)
    ###Energy

    #energybetaMPO = MBL.energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, "XY", gammescale)
    gates = MBL.gatesTEBDancilla(ancilla, h, δτ, s, "XY")
    println("gate generated")
    update = MBL.TEBDancilla!(ancilla, gates, beta, cutoff, δτ)
    println("update tebd done")
    xdata, energysiteMPO = energyagainstsite(update, h, gammescale, "XY")
    println("measure energy per site done")
    results["energy sites tebd"] = energysiteMPO
    ###DMRG
    results["liste sites"] = xdata
    psi, H = MBL.groundstateDMRG(mps, Hamiltonian, sweep_DMRG, dmax, cutoff, noise)
    xdata2, exactpersite = energyagainstsite(psi, h, gammescale, "XY")
    results["energy sites dmrg"] = exactpersite
    ###Exact energy

    exact = MBL.exactenergyXY(beta, h, γ)
    results["exact energy exact"] = exact  
    exactDMRG = mean(exactpersite)
    results["exact energy dmrg"] = exactDMRG
    exactpersite = mean(energysiteMPO)
    results["exact energy"] = exactpersite
     #####data saving
    output_data = merge(metadata, results)
    #savefile = get_savefile(output_data)
    open(savefile, "w") do io
        JSON.print(io, output_data, 4)
    end
    println("\nResults saved in $savefile")
    flush(stdout)
end

void()

println("simulation finie")
