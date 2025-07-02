#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots

############### Functions ################

"""
return the site list and the energy per site 
"""
function energyagainstsite()

end

"""
return the energy list of the site i with respect to gates time step
"""
function energyagainstdeltatime()

end

"""
return the site list and the magnet (Sz) per site
"""
function magnetagainstsite()

end
