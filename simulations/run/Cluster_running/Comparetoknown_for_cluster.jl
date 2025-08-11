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

function voidmpo()
    st, dp = MBL.section_trunc(N, gammescale)
    ancilla, s = MBL.AncillaMPO(N)
    HXX = operator(mps, h, smps, "XX")
    HXXZ = operator(mps, h, smps, "SS")

    EnergylistXX = fill(1.0,  dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    EnergylistXXZ = fill(1.0,  dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gates = gatesTEBDancilla(update, h, δτ, s, op)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = TEBDancilla!(update, gates, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, EnergylistXX[:, i] = energyMPO(update, HXX)
        _, EnergylistXXZ[:, i] = energyMPO(update, HXXZ)
        println("Average energy at β=$(betalist[i]): ", mean(Energylist; dims=1)[i])
    end

    open(savefilexxh0, "w") do io
        avgxx = vec(mean(EnergylistXX; dims=1))
        for i in eachindex(betalist)
            println(io, "$(betalist[i]) $(avgxx[i])")
        end
    end

    open(savefilexxzh0, "w") do io
        avgxxz = vec(mean(EnergylistXXZ; dims=1))
        for i in eachindex(betalist)
            println(io, "$(betalist[i]) $(avgxxz[i])")
        end
    end
end

voidmpo()

println("simulation finie")
