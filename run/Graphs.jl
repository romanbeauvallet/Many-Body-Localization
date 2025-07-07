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
function energyagainstsite(mps, h, scale)
    N = length(mps)
    start, stop = section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:2:stop)
    #@show sites
    Energypersite = Vector(undef, length(sites))
    @showprogress desc = "calcul energy over sites" for i in eachindex(sites)
        #@show i 
        Energypersite[i] = energysite(mps, sites[i], h)
    end
    return sites, Energypersite
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
    be, st = max(floor(Int, (1 + (1 - scale) * q)), 1), min(floor(Int, ((scale + 1) * q)), N)
    return be, st
end

############################ Tracer ############################

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
return the average spin  (over gammescale*length spins) against j axis for mps of length in the specif range 
"""
function magnetaverageagainstlength(j::String, gammelength::Tuple, gammescale, numbersweep, cutoff, Dmax, D0, δτ, h)
    sites = collect(gammelength[1]:1:gammelength[2])
    averagespin = Vector(undef, length(sites))
    for i in eachindex(sites)
        @show i, gammelength[2]
        mpstransit, _ = random_initialized_MPS(sites[i], D0)
        converged = tebdstepHeisenbergRow!(numbersweep, mpstransit, h, δτ, cutoff, Dmax)
        _, averagespintemp = magnetagainstsite(converged, j, gammescale)
        averagespin[i] = mean(averagespintemp)
    end
    return sites, averagespin
end

"""
return the average spin (over gammescale*length spins) value gainst the axis j for a mps  with a fixed length but updated through tebd algorithm with a number of sweep in the gammesweep range
"""
function magnetaverageagainstsweep(j::String, mps_init_sweep, gammesweep, gammescale, cutoff, Dmax, δτ, h)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2])
    realsweeplist = [gammesweep[3] for k in 1:floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3]))+1]
    realsweeplist[1] = gammesweep[1]
    #@show realsweeplist
    meanvalues = Vector(undef, length(sweeplist))
    update = mps_init_sweep
    @show update
    for p in eachindex(realsweeplist)
        @show p, realsweeplist[end]
        update = tebdstepHeisenbergRow!(realsweeplist[p], update, h, δτ, cutoff, Dmax)
        _, magnet = magnetagainstsite(update, j, gammescale)
        meanvalues[p] = mean(magnet)
    end
end


"""
return the average energy on x% of the total number of spins with respect to number of spins
"""
function energyaverageagainstsweep(mps_init_sweep, gammesweep, gammescale, cutoff, Dmax, δτ, h)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2]) #on crée la liste des sweeps (nombre total de sweep de chaque pas)
    realsweeplist = [gammesweep[3] for k in 1:floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3]))+1] #on est plus efficace si on garde le même mps et qu'on ajoute des sweep
    realsweeplist[1] = gammesweep[1]
    #@show realsweeplist
    meanvalues = Vector(undef, length(sweeplist))
    update = mps_init_sweep
    @show update
    for p in eachindex(realsweeplist)
        @show p, sweeplist[p]
        update = tebdstepHeisenbergRow!(realsweeplist[p], update, h, δτ, cutoff, Dmax)
        _, magnet = energyagainstsite(update, h, gammescale)
        #@show magnet
        meanvalues[p] = mean(magnet)
    end
    return sweeplist, meanvalues
end

"""
return the average energy (over gammescale*length spins) for a fixed number of tebd steps but with different length of mps
"""
function energyaverageagainstlength(gammelength::Tuple, gammescale, numbersweep, cutoff, Dmax, D0, δτ, h)
    sites = collect(gammelength[1]:1:gammelength[2])
    averageenergy = Vector(undef, length(sites))
    for i in eachindex(sites)
        @show i, gammelength[i]
        mpstransit, _ = random_initialized_MPS(sites[i], D0)
        converged = tebdstepHeisenbergRow!(numbersweep, mpstransit, h, δτ, cutoff, Dmax)
        _, averageenergytemp = energyagainstsite(converged, j, gammescale)
        averageenergy[i] = mean(averageenergytemp)
    end
    return sites, averageenergy
end