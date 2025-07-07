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
n_sweep = 100
cutoff = 1e-15
Dmax = 500
Beta = n_sweep * δτ
gammelength = (div(N, 10), N)
gammescale = 0.5
j = "z"
gammesweep = (1000, 3000, 500) #(start, stop, step)

params = Dict(
    "N" => N,
    "J" => 1,
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
final = tebdstepHeisenbergRow!(5000, mps_random_debut, h, δτ, cutoff, Dmax)
_, energyinf = mean(energyagainstsite(final, h, gammescale))


xdata, ydata = MBL.energyaverageagainstsweep(mps_random_debut, gammesweep, gammescale, cutoff, Dmax, δτ, h)

results = Dict(



)


simulation = Dict(
    "parameters" => params,
    "results" => results
)

open("simulation_output.json", "w") do io
    JSON.print(io, simulation_data)
end