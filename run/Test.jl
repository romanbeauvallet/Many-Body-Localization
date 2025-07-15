#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots
using ITensorMPS

using Base.Threads
println("Nombre de threads disponibles : ", nthreads())

################ Parameters ###############

N = 10
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
Dmax = 300
betamax = 5
step = 0.5
Beta = n_sweep * δτ
gammescale = 0.6
j = "z"
γ=0.0
betalist = collect(0:step:betamax)
################# try to add quantum number conservation
test, s = MBL.AncillaMPO(N)
@show siteinds(test)
#@show energyagainstsite(test, h, gammescale)
beta, energy = MBL.energyforbetalist(betamax, step, test, δτ, h, s, cutoff, "SS")
exactenergy = [MBL.exactenergyXY(b, h, γ) for b in betalist]
gr()
scatter(beta, energy, ylabel="<H>/N", title="N=$N, cutoff=$cutoff, δτ=$δτ, step=$step", xlabel="β", label="tebd")
scatter!(betalist, exactenergy, label="exact energy")
plot!(betalist, -1 .* betalist, label="-betalist")
#hline!([1/4-log(2)], y="exact")
