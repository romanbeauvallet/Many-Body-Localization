#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Plots
using Statistics
############### Functions ################

"""
return the site list and the energy per site 
"""
function energyagainstsite(mps, h)
    N = length(mps)
    sites = collect(1:2:N-2)
    Energypersite = Vector(undef, length(sites))
    @showprogress desc = "calcul energy over sites" for i in sites
        Energypersite[i] =  energysite(mps, i, h)
    end
    return sites, Energypersite
end

"""
return the energy list of the site i with respect to gates time step
"""
function energyagainstdeltatime!(site_measure, gamme::Tuple, mpsinit, step, numbersweep, cutoff, Dmax)
    timesteplist = reverse(collect(gamme[1]:step:gamme[2]))
    EnergyList = Vector(undef, length(timesteplist))
    @showprogress for j in eachindex(timesteplist)
        mpsinit = tebdstepHeisenbergRow!(numbersweep, mpsinit, h, timesteplist[j], cutoff, Dmax)
        e = energysite(mpsinit, site_measure, h)
        EnergyList[j] = e
    end
    return timesteplist, EnergyList
end

"""
return the site list and the magnet (Sz) per site
"""
function magnetagainstsite(mps, j::String, scale)
    N = length(mps)
    start, stop = section_trunc(N, scale)
    sites = collect(start:1:stop)
    Magnetpersite = Vector{}(undef, length(sites))
    @showprogress desc = "calcul magnet over sites" for p in eachindex(sites)
        Magnetpersite[p] = measure_S(mps, p, j)
    end
    return sites, Magnetpersite
end
"""
N -- length of a vector
scale -- between 0 and 1, overlap between the sliced list and the original list

return the indexes that slices a list of length N with the overlap scale
"""
function section_trunc(N, scale)
    q = div(N, 2)
    be, st = max(floor(Int, (1 + (1 - scale) * q)), 1) , min(floor(Int, ((scale + 1) * q)), N)
    return be, st
end

"""
return the average spin against j axis for mps of length N
"""
function averagespinoverlength(j::String, gammelength::Tuple, gammescale, numbersweep, cutoff, Dmax, D0, δτ, h)
    sites = collect(gammelength[1]:1:gammelength[2])
    averagespin = Vector(undef, length(sites))
    for i in gammelength[1]:1:gammelength[2]
        @show i
        mpstransit, _ = random_initialized_MPS(i, D0)
        converged = tebdstepHeisenbergRow!(numbersweep, mpstransit, h, δτ, cutoff, Dmax)
        _, averagespintemp = magnetagainstsite(converged, j, gammescale)
        averagespin[i] = mean(averagespintemp)
    end
    return sites, averagespin
end

function magnetaverageagainstsweep(j::String, mps_init_sweep, gammesweep, gammescale, cutoff, Dmax, δτ, h)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2])
    realsweeplist = [gammesweep[3] for k in 1:((gammesweep[2]-gammesweep[1])/gammesweep[3])+1]
    realsweeplist[1] = gammesweep[1]
    #@show realsweeplist
    meanvalues = Vector(undef, length(sweeplist))
    update = mps_init_sweep
    @show update
    for p in eachindex(realsweeplist)
        @show p, realsweeplist[p]
        update = tebdstepHeisenbergRow!(realsweeplist[p], update, h, δτ, cutoff, Dmax)
        _, magnet = magnetagainstsite(update, j, gammescale)
        meanvalues[p] = mean(magnet)
    end
end