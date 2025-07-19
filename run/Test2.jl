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

N = 18
J = 1
h = 10
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 3000
cutoff = 1e-15
dmax = 300
betamax = 1
stepbeta = 1
Beta = n_sweep * δτ
gammescale = 0.6
noise = 1e-8
n_sweepDMRG = 20
j = "z"
γ = 0.0
init = 1234
betalist = collect(0:stepbeta:betamax)

# ============================== DATA
filename = joinpath("analyse_simulations_julia", "DATA", "spec_XX_N18.json")
json_string = read(filename, String)
input = JSON.parse(json_string)

ancilla, s = MBL.AncillaMPO(N)
mps, smps = neelstate(N-1)
println("init")
_, gates, champ = MBL.evolutionwithrandomdisordergates(init, ancilla, s, h, δτ)
println("gates generated")
update = MBL.TEBDancilla!(ancilla, gates, betamax, cutoff, δτ)
println("MPO done")
ampo = AutoMPO()
for j in 1:(N - 2)
    # ampo .+= (operator_coefficient, operator_name, site_index, ...)
    ampo .+= (1.0, "Sz", j, "Sz", j + 1)
    ampo .+= (0.5, "S+", j, "S-", j + 1)
    ampo .+= (0.5, "S-", j, "S+", j + 1)
    ampo .+= (champ[j], "Sz", j)
end
ampo .+= (champ[N-1], "Sz", N-1)
# Convert the AutoMPO to an MPO
H = MPO(ampo, smps)
psi, E = MBL.groundstateDMRG(mps, H, n_sweepDMRG, dmax, cutoff, noise)
println("update done")
st, dp= MBL.section_trunc(N-1, gammescale)
LDMRG = collect(st:dp)
st, dp= MBL.section_trunc(N, gammescale)
L = collect(st:dp)
_, xdataDMRG = magnetagainstsite(psi, j, gammescale)
_, xdata = magnetagainstsite(update, j, gammescale)

gr()
p = plot()
scatter!(p, LDMRG, xdataDMRG, label="DMRG")
scatter!(p, L, xdata)
display(p)


#for i in 1:length(L)-1  

    #T = randn(2,2,2,2)
    #s1, s2 = siteind(update, L[i]), siteind(update, L[i+1])
    #p1, p2 = siteind(psi, L[i]), siteind(psi, L[i+1])
    #S = MBL.randomoperator(T,s1,s2)    
    #P = MBL.randomoperator(T, p1,p2)
    #push!(xdata, energysiteMPOdisorder(update, L[i], S))
    #push!(xdataDMRG, energysiteMPOdisorder(psi, L[i], P))
#end
