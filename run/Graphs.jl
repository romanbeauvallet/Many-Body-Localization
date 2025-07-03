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
    Energypersite = Vector{}()
    @showprogress desc = "calcul energy over sites" for i in sites
        push!(Energypersite, energysite(mps, i, h))
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
        e = energysite(mpsinit, site_measure, h)
        push!(EnergyList, e)
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
    Magnetpersite = Vector{}()
    @showprogress desc = "calcul magnet over sites" for i in sites
        push!(Magnetpersite, measure_S(mps, i, j))
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
    be, st = floor(Int,(1 + (1 - scale) * q)), floor(Int,((scale + 1) * q))
    return be, st
end

"""
return the average spin against j axis for mps of length N
"""
function averagespinoverlength(j, gammelength::Tuple, gammescale, numbersweep, cutoff, Dmax, D0, δτ, h)
    sites = collect(gammelength[1]:1:gammelength[2])
    averagespin = Vector{}()
    for i in gammelength[1]:1:gammelength[2]
        @show i 
        mpstransit, _ = random_initialized_MPS(i, D0)
        converged = tebdstepHeisenbergRow!(numbersweep, mpstransit, h, δτ, cutoff, Dmax)
        _, averagespintemp = magnetagainstsite(converged, j, gammescale)
        push!(averagespin, mean(averagespintemp[gammemeasure[1]:gammemeasure[2]]))
    end
    return sites, averagespin
end

function magnetaverageagainstsweep(j, mps_init, gammesweep, gammescale, cutoff, Dmax, δτ, h)
    sweeplist = collect(gammesweep[1]:1:gammesweep[2])
    meanvalues = Vector(undef, length(sweeplist))
    update = mps_init
    for p in eachindex(sweeplist)
        update = tebdstepHeisenbergRow!(sweeplist[p], update, h, δτ, cutoff, Dmax)
        _, magnet = magnetagainstsite(update, j, gammescale)
        meanvalues[p] = mean(magnet)
    end
end