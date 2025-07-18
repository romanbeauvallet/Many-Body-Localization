#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots
using ITensorMPS
using Statistics
using Base.Threads
using LaTeXStrings
println("Nombre de threads disponibles : ", nthreads())
using JSON
################ Parameters ###############

N = 50
J = 1
h0 = 0
h = 100
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 5000
cutoff = 1e-15
dmax = 300
betamax = 10
stepbeta = 0.1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 10
j = "z"
γ = 0.0
init = 1234
betalist = collect(0:stepbeta:betamax)

# ============================== DATA
ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N)
H = hamiltonianHeisenberg(mps, h, smps)
ydata=Vector{Vector{Float64}}()
function void()
    list= [100 for i in 1:n_sweep/100]
    update = mps
    for i in list
        update = tebdstepHeisenbergRow!(i, update, h, δτ, cutoff, dmax, "SS")
        _ , m = magnetagainstsite(update, j, gammescale)
        push!(ydata, m)
    end
    return ydata
end

ydata = MBL.magnetforbestalistdisorder(betalist, ancilla, δτ, h, s, cutoff, gammescale, init, j)
matrix = abs.(hcat(ydata...))
st, dp = MBL.section_trunc(N, gammescale)
L = collect(st:dp)
inter = matrix./maximum(matrix)
gr()

heatmap(betalist, L, inter,
        xlabel="β", ylabel="L",
        title="Magnétisation",
        colorbar_title="M",
        aspect_ratio=:auto,
        c=:viridis)  # palette de couleur
#scatter!(p, xdata2, ydata2, label ="TEBD step=0.5")
#hline!(p, [1/4-log(2)], label="exact energy at 0K without disorder")

