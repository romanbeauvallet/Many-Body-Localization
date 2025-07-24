#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using Distributed
using Base.Threads
using LaTeXStrings
println("Nombre de threads disponibles : ", nthreads())
using JSON
@everywhere begin
    using ITensors
    using MBL
    using ProgressMeter
    using Plots
    using ITensorMPS
    using Statistics
end

# Ajouter des processus seulement si tu es en local
if nprocs() == 1
    addprocs(4)  # ne sera pas utilisé sur SLURM
end
################ Parameters ###############

N = 20
J = 1
h = 10
h0 = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
dmax = 300
betamax = 5
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 20
j = "z"
γ = 0.0
betalist = collect(0:stepbeta:betamax)
seed1 = 1234432
seed2 = 4578
# ============================== DATA
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)

println("init done")
#ydata, _ = MBL.magnetforbestalistdisorder(betalist, ancilla, δτ, h, s, cutoff, gammescale, seed1, "z")
result = pmap(
    t -> MBL.magnetforbestalist(t, mps, δτ, h, cutoff, gammescale, "SS", "z", dmax),
    betalist,
)
ydata2 = MBL.magnetforbestalist(betalist, mps, δτ, h, cutoff, gammescale, "SS", "z", dmax)
st, dp = MBL.section_trunc(N, gammescale)
L = collect(st:1:dp)
@show magnetagainstsite(ancilla, "z", gammescale)
@show magnetagainstsite(mps, "z", gammescale)

#psi, e = MBL.groundstateDMRG(mps, Hamiltonian, n_sweepDMRG, dmax, cutoff, noise)
#xdata, exactenergysite = magnetagainstsite(psi, "z", gammescale)

#@show update1[5] == update2[5]
gr()
heatmap(
    betalist,
    L,
    hcat(ydata2...);
    xlabel="β",
    ylabel="site list",
    title="N=$N, cutoff=$cutoff, δτ=$δτ, h=$h",
    colorbar=true,
    color=:viridis,  # Palette de couleurs
)