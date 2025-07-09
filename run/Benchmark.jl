#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using ProgressMeter
using Plots
using JSON
################# Parameters ###############
# Lire les arguments
input_file = ARGS[1]
output_file = ARGS[2]
sim_id = parse(Int, ARGS[3])

params = JSON.parsefile(input_file)

# Extraire les paramètres
N = params["N"]
D0 = params["D0"]
h = params["disorder"]
δτ = params["Trotter-Suzuki step"]
Dmax = params["max bond dimension"]
gammesweep = params["nsweep range"]
gammescale = params["gammescale"]
cutoff = params["cutoff"]
n_sweep = params["fixed number of sweep"]
j = params["axis"]


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
    "input" => params,
    "results" => results
)

json_path = joinpath(@__DIR__,"..", "analyse_simulations_julia", "data.json")

open(json_path, "w") do io
    write(io, JSON.json(simulation, 2))
end