#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using ProgressMeter
using Plots

################# Parameters ###############
N = 100
J = 1
h = 0
δτ = 1e-3
D0 = 10
site_measure = div(N, 2)
n_sweep = 400
cutoff = 1e-15
Dmax = 300
Beta = n_sweep * δτ

################# Scaling ################

gammelength = (div(N, 10), N)
@show gammelength
gammescale = 0.5
j = "z"

xdata, ydata = MBL.averagespinoverlength(j, gammelength, gammescale, n_sweep, cutoff, Dmax, D0, δτ, h)
gr()

scatter(xdata, ydata)
