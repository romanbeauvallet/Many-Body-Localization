#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Libraries #################
using MBL
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
betamax1 = input_data["beta max 1"]
betamax2 = input_data["beta max 2"]
step1 = input_data["beta step 1"]
step2 = input_data["beta step 2"]
savefilexxh0 = String(input_data["savefilexxh0"])
savefilexxzh0 = String(input_data["savefilexxzh0"])

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
    "maximum bond dimension per tebd step" => nothing,
)
println("\nmetadata:")
display(metadata)

resultsdmrg = Dict{String,Any}(
    "energy sites tebd" => nothing,
    "energy sites dmrg" => nothing,
    "exact energy" => nothing,
    "exact energy exact" => nothing,
    "exact energy dmrg" => nothing,
    "liste sites" => nothing,
    "magnet sites tebd" => nothing,
    "magnet sites dmrg" => nothing,
)

function voiddmrg()
    ancilla, s = MBL.AncillaMPO(N)
    mps, smps = neelstate(N)
    Hamiltonian = MBL.operator(mps, h, smps, "XX")
    ###Energy

    #energybetaMPO = MBL.energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, "XY", gammescale)
    gates = MBL.gatesTEBDancilla(ancilla, h, δτ, s, "XY")
    println("gate generated")
    update = TEBDancilla!(ancilla, gates, beta, cutoff, δτ, Dmax)
    println("update tebd done")
    xdata, energysiteMPO = energyagainstsite(update, h, gammescale, "XX")
    results["energy sites tebd"] = energysiteMPO
    println("measure energy per site done")
    _, magnetsiteMPO = magnetagainstsite(update, "z", gammescale)
    results["magnet sites tebd"] = magnetsiteMPO
    println("measure magnet per site done")
    ###DMRG
    results["liste sites"] = xdata
    psi, H = groundstateDMRG(mps, Hamiltonian, sweep_DMRG, Dmax, cutoff, noise)
    xdata2, exactpersite = energyagainstsite(psi, h, gammescale, "XX")
    _, magnetDMRG = magnetagainstsite(psi, "z", gammescale)
    results["energy sites dmrg"] = exactpersite
    results["magnet sites tebd"] = magnetDMRG
    println("DMRG done")
    ###Exact energy

    exact = MBL.exactenergyXX(beta, h; γ)
    results["exact energy exact"] = exact
    exactDMRG = mean(exactpersite)
    results["exact energy dmrg"] = exactDMRG
    exactpersite = mean(energysiteMPO)
    results["exact energy"] = exactpersite
    #####data saving
    output_data = merge(metadata, results)
    #savefile = get_savefile(output_data)
    open(savefile, "w") do io
        JSON.print(io, output_data, 2)
    end
    println("\nResults saved in $savefile")
    return flush(stdout)
end

gammelength = [10, 100, 10]
sites = collect(gammelength[1]:gammelength[3]:gammelength[2])

function voidmpo()
    averageenergysites = Array{Float64}(undef, length(sites))
    averageenergympoXX = Array{Float64}(undef, length(sites))
    averageenergympoXXZ = Array{Float64}(undef, length(sites))
    @showprogress for i in eachindex(sites)
        mpstransit, s = random_initialized_MPS(sites[i], D0)
        HXX = operator(mpstransit, h, s, "XX")
        HXXZ = operator(mpstransit, h, s, "SS")
        gates = gateTrotterSuzukirow(mpstransit, h, δτ, op)
        converged = tebdevolutionrow!(numbersweep, mpstransit, gates, cutoff, dmax)
        _, averageenergytemp = energyagainstsite(converged, h, gammescale, op)
        energyMPOxx = energyMPO(converged, HXX)/sites[i]
        energyMPOxxz = energyMPO(converged, HXXZ)/sites[i]
        averageenergysites[i] = mean(averageenergytemp)
        averageenergympoXX[i] = energyMPOxx
        averageenergympoXXZ[i] = energyMPOxxz
    end

    open(savefilexxh0, "w") do io
        for i in eachindex(sites)
            println(io, "$(sites[i]) $(averageenergympoXX[i])")
        end
    end

    open(savefilexxzh0, "w") do io
        for i in eachindex(sites)
            println(io, "$(sites[i]) $(averageenergympoXXZ[i])")
        end
    end
end

voidmpo()

println("simulation finie")
