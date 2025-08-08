#!/usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
if Sys.isapple()
    ENV["JULIA_BLAS_VENDOR"] = "accelerate"
end
import LinearAlgebra: BLAS
@show BLAS.vendor()          # should print :accelerate
BLAS.set_num_threads(8)      # tune; 1–4 often best
using ProgressMeter
using JSON
using Random
using Statistics
using Dates
using LinearAlgebra
using Pkg
# Ensure we use the simulation project environment
#Pkg.activate(joinpath(@__DIR__, "..")) not worth it because julia is started in the project with julia --project=.
try
    using HDF5
catch
    @warn "HDF5 not available; continuing without it"
end
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
#@show BLAS.get_config()
#@show BLAS.get_num_threads()
@show Threads.nthreads()
#@show LinearAlgebra.BLAS.get_num_threads()
Pkg.status()
# ===================== parameters

N = 100
J = 1
h = 0.0
δτ = 1e-3
dmax = 300
gammescale = 0.8 #missing in the output data (found in the input data on the cluster)
cutoff = 1e-15
j = "z"
betamax1 = 10
step1 = 0.1
step2 = 1
betamax2 = 15
noise = 1e-8
n_sweep = 20

base_plots = joinpath(@__DIR__, "..", "..", "analyse_simulations_julia", "DATA_Local", "Plots_txt")
dir_xx = joinpath(base_plots, "XX")
dir_xxz = joinpath(base_plots, "XXZ")
mkpath(dir_xx)
mkpath(dir_xxz)

savefilexxh0 = joinpath(dir_xx, "energyxxh0.txt")
savefilexxzh0 = joinpath(dir_xxz, "energyxxzh0.txt")
savefileexactenergy = joinpath(dir_xx, "energyexactxxh0.txt")
savefilexxzdmrgh0 = joinpath(dir_xxz, "energyxxzDMRGh0.txt")
savefilexxdmrgh0 = joinpath(dir_xx, "energyexactxxDMRGh0.txt")

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
    "gammescale" => gammescale,
    "disorder" => h,
    "number of spins measured" => gammescale,
    "beta list values" => betalist,
    "noise" => noise, 
    "nombre de sweep" =>n_sweep
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




