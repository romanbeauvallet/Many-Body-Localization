#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using Base.Threads
using LaTeXStrings
println("Nombre de threads disponibles : ", nthreads())
using JSON
using ITensors
using MBL
using ProgressMeter
using Plots
using ITensorMPS
using Statistics
using Random
using Distributions
# =========================
N = 20
Ndisorder = 30
J = 1
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
dmax = 300
# =========================
betamax = 5
stepbeta = 1
betaeff = n_sweep * δτ
gammescale = 0.7
noise = 1e-8
n_sweepDMRG = 20
j = "z"
γ = 0.0
# =========================
h = 20
centerpic = 3.5
init = 4375432
# =========================
betalist = collect(0:stepbeta:betamax)
betafixe = 5
place = findfirst(==(betafixe), betalist)
@assert betafixe in betalist "La valeur $betafixe n'est pas dans la liste"

fullseed = 1234567890
rng = MersenneTwister(fullseed)
seeds = rand(rng, Uniform(-h, h), N - 1)
# =========================
ancilla, s = MBL.AncillaMPO(N)
disorder = sort(MBL.rejection_sample(Ndisorder, h, centerpic; σ=0.03h, A=1.0, init))
st, dp = MBL.section_trunc(N, gammescale)
L = collect(st:dp)
#pour chaque champ, pour chaque beta : mesurer l'énergie de chaque site, moyenne et ecart type 

Magnetlist = Array{Float64}(undef, length(L), length(betalist), length(disorder))

# ========================= SIMU 

@showprogress desc ="run over disorder" for i in eachindex(disorder)
    value, _ = MBL.magnetforbestalistdisorder(betalist, ancilla, δτ, disorder[i], s, cutoff, gammescale, init, j, dmax)
    Magnetlist[:, :, i] = value
end 
Meanvalue = mean(Magnetlist; dims=1)
Stdvalue = std(Magnetlist; dims=1)
@show disorder

output_data = Dict{String,Any}(
    "energy" => Energylist,
    "beta list" => betalist,
    "champ" => disorder, 
    "site list" => L
)

#savefile = joinpath("../analyse_simulations_julia/DATA_Local/", "disorderbeta.json")

#open(savefile, "w") do io
    #JSON.print(io, output_data, 4)
#end
#println("\nResults saved in $savefile")


gr()
default(fontfamily="Computer Modern")
p = plot()
scatter!(p, disorder, vec(Meanvalue[:, place, :]), label="mean value", xlabel="disorder magnitude", ylabel="Sz expectation value", title="N=$N, cutoff=$cutoff, δτ=$δτ, β=$betafixe")
scatter!(p, disorder, vec(Stdvalue[:, place, :]), label="standart deviation")
display(p)
savefig("run/Plots/meanandstdmagnetSS.pdf")





