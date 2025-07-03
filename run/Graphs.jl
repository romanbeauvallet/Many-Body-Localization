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
function energyagainstsite(mps, h)
    N = length(mps)
    sites = collect(1:2:N-2)
    Energypersite = Vector{}()
    @showprogress for i in sites
        push!(Energypersite, energysite!(mps, i, h))
    end
    return sites, Energypersite
end

"""
return the energy list of the site i with respect to gates time step
"""
function energyagainstdeltatime!(site_measure, gamme::Tuple, mpsinit, step, numbersweep, cutoff, Dmax)
    timesteplist = reverse(collect(gamme[1]:step:gamme[2]))
    EnergyList = Vector{}()
    @showprogress for j in eachindex(timesteplist)
        mpsinit = tebdstepHeisenbergRow!(numbersweep, mpsinit, h, timesteplist[j], cutoff, Dmax) 
        e= energysite(mpsinit, site_measure, h)
        push!(EnergyList, e)
    end
    return timesteplist, EnergyList
end 

"""
return the site list and the magnet (Sz) per site
"""
function magnetagainstsite(mps)
    N = length(mps)
    sites = collect(1:1:N)
    Magnetpersite = Vector{}()
    @showprogress for i in sites
        push!(Magnetpersite, measure_Sz(mps, i))
    end
    return sites, Magnetpersite
end

