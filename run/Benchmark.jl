#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using ProgressMeter
using Plots
using JSON
################# Parameters ###############
N = 100
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
Dmax = 500
Beta = n_sweep * δτ
gammelength = (div(N, 10), N)
gammescale = 0.5
j = "z"
gammesweep = (100, 4000, 500) #(start, stop, step)

params = Dict(
    "N" => N,
    "J" => 1,
    "D0" => D0,
    "site measure" => site_measure,
    "cutoff" => cutoff,
    "max bond dimension" => Dmax,
    "Trotter-Suzuki step" => δτ,
    "disorder" => h,
    "nsweep range" => gammesweep,
    "length range" => gammelength,
    "fixed number of sweep" => n_sweep
)

################# Scaling ################
mps_random_debut, _ = random_initialized_MPS(N, D0)

x1data, y1data = MBL.energyaverageagainstsweep(mps_random_debut, gammesweep, gammescale, cutoff, Dmax, δτ, h)
x2data, y4data = MBL.magnetaverageagainstlength(j, gammelength, gammescale, n_sweep, cutoff, Dmax, D0, δτ, h)
_, y2data = MBL.magnetaverageagainstsweep(j, mps_random_debut, gammesweep, gammescale, cutoff, Dmax, δτ, h)
_, y3data = MBL.energyaverageagainstlength(gammelength, gammescale, n_sweep, cutoff, Dmax, D0, δτ, h)

#@show "it works"

results = Dict(
    "sweep list" => x1data,
    "energy sweep list" => y1data,
    "length list" => x2data,
    "magnetization sweep list" => y2data,
    "energy length list" => y3data,
    "magnetization length list" => y4data
)

simulation = Dict(
    "parameters" => params,
    "results" => results
)

json_path = joinpath(@__DIR__, "..", "analyse_simulations_julia", "data.json")

open(json_path, "w") do io
    write(io, JSON.json(simulation, 2))
end