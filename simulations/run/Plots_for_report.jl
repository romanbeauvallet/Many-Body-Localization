#!/usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Libraries #################
using MBL
import LinearAlgebra: BLAS
@show BLAS.vendor()
using ProgressMeter
using JSON
using Random
using Statistics
using Dates
using LinearAlgebra
using Pkg
# Ensure we use the simulation project environment
try
    using PProf
    using Profile
catch
    @warn "PProf/Profile not available; profiling disabled"
end
# ===================== log
println(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
println("Julia $VERSION")
@show Base.julia_cmd()
@show Threads.nthreads()
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

N = input_data["N"]
J = input_data["J"]
γ = input_data["gamma"]
h = input_data["disorder"]
δτ = input_data["Trotter-Suzuki step"]
dmax = input_data["max bond dimension"]
gammescale = input_data["gammescale"]
cutoff = input_data["cutoff"]
j = input_data["axis"]
noise = input_data["noise"]
n_sweep = input_data["number sweep dmrg"]
betamax1 = input_data["beta fixe"]
betamax2 = input_data["beta fixe"]
=#

N = 100
J = 1
γ = 0.0
h = 0.0
δτ = 1e-3
dmax = 400
gammescale = 0.8
cutoff = 1e-15
j = "z"
noise = 1e-8
n_sweep = 10
betamax1 = 10
betamax2 = 30
step1 = 0.1
step2 = 1

savefilexxh0 = joinpath("..", "analyse_simulations_julia", "DATA_Local", "Plots_txt","XX", "energyxxh0.txt")
savefilexxzh0 = joinpath("..", "analyse_simulations_julia", "DATA_Local", "Plots_txt","XXZ", "energyxxzh0.txt")
savefileexactenergy = joinpath("..", "analyse_simulations_julia", "DATA_Local", "Plots_txt","XX", "energyexactxxh0.txt")
savefilexxzdmrgh0 = joinpath("..", "analyse_simulations_julia", "DATA_Local", "Plots_txt","XXZ", "energyxxzDMRGh0.txt")
savefilexxdmrgh0 = joinpath("..", "analyse_simulations_julia", "DATA_Local", "Plots_txt","XX", "energyexactxxDMRGh0.txt")

# ====================================== Dict

betalist1 = collect(0:step1:betamax1)
betalist2 = collect(betamax1:step2:betamax2)

betalist = vcat(betalist1, betalist2)

realbetalist = pushfirst!(diff(betalist), 0)

metadata = Dict{String, Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "Given maximal bond dimension" => dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "disorder" => h,
    "number of spins measured" => gammescale,
    "beta list values" => betalist,
    "noise DMRG" => noise, 
    "nombre de sweep DMRG" =>n_sweep
)
println("\nmetadata:")
display(metadata)

# ====================================== begin

st, dp = MBL.section_trunc(N, gammescale)
ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)
HXX = operator(mps, h, smps, "XX")
HXXZ = operator(mps, h, smps, "SS")

energyxx = MBL.energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, "XX", gammescale, dmax)
open(savefilexxh0, "w") do io
    avgxx = vec(mean(energyxx; dims=1))
    for i in eachindex(betalist)
        println(io, "$(betalist[i]) $(avgxx[i])")
    end
end
psi, energyXXDMRG = MBL.groundstateDMRG(mps, HXX, n_sweep, dmax, cutoff, noise)
open(savefilexxdmrgh0, "w") do io
    for i in eachindex(betalist)
        println(io, "$(betalist[i]) $(energyXXDMRG)")
    end
end
psi, energyXXZDMRG = MBL.groundstateDMRG(mps, HXXZ, n_sweep, dmax, cutoff, noise)
open(savefilexxzdmrgh0, "w") do io
    for i in eachindex(betalist)
        println(io, "$(betalist[i]) $(energyXXZDMRG)")
    end
end
energyxxz = MBL.energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, "SS", gammescale, dmax)
open(savefilexxzh0, "w") do io
    avgxxz = vec(mean(energyxxz; dims=1))
    for i in eachindex(betalist)
        println(io, "$(betalist[i]) $(avgxxz[i])")
    end
end
exactenergy = [MBL.exactenergyXX(β, h; γ=0.0) for β in betalist]
open(savefileexactenergy, "w") do io
    for i in eachindex(betalist)
        println(io, "$(betalist[i]) $(exactenergy[i])")
    end
end




