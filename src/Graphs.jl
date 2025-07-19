#!usr/bin/env julia

################ Librairies #############
using ITensors
using MBL
using ProgressMeter
using Statistics
############### Functions ################

# ============================================= Energy
"""
mps -- MPS
h -- disorder
scale -- percentage of the chain on which you measure

return the site list and the energy per site 
"""
function energyagainstsite(mps, h, scale)
    N = length(mps)
    start, stop = section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    Energypersite = Vector(undef, length(sites))
    @showprogress desc = "calcul energy over sites" for i in eachindex(sites)
        Energypersite[i] = energysite(mps, sites[i], h)
    end
    return sites, Energypersite
end

"""
return the average energy on x% of the total number of spins with respect to number of spins
"""
function energyaverageagainstsweep(mps_init_sweep, gammesweep, gammescale, cutoff, Dmax, δτ, h, op::String)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2]) #on crée la liste des sweeps (nombre total de sweep de chaque pas)
    realsweeplist = [gammesweep[3] for k in 1:floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3]))+1] #on est plus efficace si on garde le même mps et qu'on ajoute des sweep
    realsweeplist[1] = gammesweep[1]
    meanvalues = Vector(undef, length(sweeplist))
    update = mps_init_sweep
    for p in eachindex(realsweeplist)
        update = tebdstepHeisenbergRow!(realsweeplist[p], update, h, δτ, cutoff, Dmax, op)
        _, magnet = energyagainstsite(update, h, gammescale)
        meanvalues[p] = mean(magnet)
    end
    return sweeplist, meanvalues
end

"""
return the average energy (over gammescale*length spins) for a fixed number of tebd steps but with different length of mps
"""
function energyaverageagainstlength(gammelength::Tuple, gammescale, numbersweep, cutoff, Dmax, D0, δτ, h, op::String)
    sites = collect(gammelength[1]:gammelength[3]:gammelength[2])
    averageenergy = Vector(undef, length(sites))
    @showprogress for i in eachindex(sites)
        mpstransit, _ = random_initialized_MPS(sites[i], D0)
        converged = tebdstepHeisenbergRow!(numbersweep, mpstransit, h, δτ, cutoff, Dmax, op)
        _, averageenergytemp = energyagainstsite(converged, h, gammescale)
        averageenergy[i] = mean(averageenergytemp)
    end
    return sites, averageenergy
end

##################### Tracer Finite Temperature #################
"""
betamax -- maximal beta you want to reach
step -- step in the beta list
ancilla -- MPS

return the average energy of the MPS at temperature β ∈ [0:step:betamax] for the operator op, computed with the op as an MPO
"""
function energyforbetalistMPO(betamax, step, ancilla, δτ, h, s, cutoff, op)
    N = length(ancilla)
    betalist = collect(0:step:betamax)
    realbetalist = reverse(push!(diff(betalist), 0))
    Energylist = Vector{}(undef, length(realbetalist))
    if op == "XY"
        H = hamiltonianXY(ancilla, h, s)
    elseif op == "SS"
        H = hamiltonianHeisenberg(ancilla, h, s)
    end
    update = ancilla
    gates = gatesTEBDancilla(update, h, δτ, s, op)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gates, realbetalist[i], cutoff, δτ)
        Energylist[i] = MBL.energyMPO(update, H) / N
    end
    return betalist, Energylist
end

"""
betamax -- maximal beta you want to reach
step -- step in the beta list
ancilla -- MPS

return the MPS average energy measured with gates on gammescale*length(MPS) number of sites taken from the MPS center
"""
function energyforbetalist(betalist, ancilla, δτ, h, s, cutoff, op::String, gammescale)
    realbetalist = pushfirst!(diff(betalist), 0)
    Energylist = Vector{}(undef, length(realbetalist))
    update = ancilla
    gates = gatesTEBDancilla(update, h, δτ, s, op)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gates, realbetalist[i] / 2, cutoff, δτ)
        _, Energylist[i] = energyagainstsiteMPO(update, h, gammescale, op)
    end
    return Energylist
end

"""
return the energy list of the site i with respect to gates time step
"""
function energyagainstdeltatime!(site_measure, gamme::Tuple, mpsinit, step, numbersweep, cutoff, Dmax, op::String)
    timesteplist = reverse(collect(gamme[1]:step:gamme[2]))
    EnergyList = Vector(undef, length(timesteplist))
    @showprogress for j in eachindex(timesteplist)
        mpsinit = tebdstepHeisenbergRow!(numbersweep, mpsinit, h, timesteplist[j], cutoff, Dmax, op)
        e = energysite(mpsinit, site_measure, h)
        EnergyList[j] = e
    end
    return timesteplist, EnergyList
end


####################### Random disorder #######################
"""
return the betalist and the energy list with a random uniform on each site
"""
function energyforbestalistdisorder(betalist, ancilla, δτ, h, s, cutoff, gammescale, init)
    realbetalist = pushfirst!(diff(betalist), 0)
    Energylist = Vector{}(undef, length(realbetalist))
    update = ancilla
    gatesmeasure, gatesevolve, _ = evolutionwithrandomdisordergates(init::Int64, update, s, h, δτ)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ)
        _, Energylist[i] = energyagainstsiteMPOdisorder(update, gatesmeasure, gammescale)
    end
    return Energylist
end

# =============================================== Magnetization
"""
return the site list and the magnet (Sz) per site
"""
function magnetagainstsite(mps, j::String, scale)
    N = length(mps)
    start, stop = section_trunc(N, scale)
    sites = collect(start:1:stop)
    Magnetpersite = Vector{Float64}(undef, length(sites))
    @showprogress desc = "calcul magnet over sites" for p in eachindex(sites)
        Magnetpersite[p] = measure_S(mps, p, j)
    end
    return sites, Magnetpersite
end

"""
return the average spin  (over gammescale*length spins) against j axis for mps of length in the specif range 
"""
function magnetaverageagainstlength(j::String, gammelength::Tuple, gammescale, numbersweep, cutoff, Dmax, D0, δτ, h, op::String)
    sites = collect(gammelength[1]:1:gammelength[2])
    averagespin = Vector(undef, length(sites))
    for i in eachindex(sites)
        @show i, gammelength[2], sites[i]
        mpstransit, _ = random_initialized_MPS(sites[i], D0)
        converged = tebdstepHeisenbergRow!(numbersweep, mpstransit, h, δτ, cutoff, Dmax, op)
        _, averagespintemp = magnetagainstsite(converged, j, gammescale)
        averagespin[i] = mean(averagespintemp)
    end
    return sites, averagespin
end

"""
return the average spin (over gammescale*length spins) value gainst the axis j for a mps  with a fixed length but updated through tebd algorithm with a number of sweep in the gammesweep range
"""
function magnetaverageagainstsweep(j::String, mps_init_sweep, gammesweep, gammescale, cutoff, Dmax, δτ, h, op::String)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2])
    realsweeplist = [gammesweep[3] for k in 1:floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3]))+1]
    realsweeplist[1] = gammesweep[1]
    meanvalues = Vector(undef, length(realsweeplist))
    update = mps_init_sweep
    for p in eachindex(realsweeplist)
        update = tebdstepHeisenbergRow!(realsweeplist[p], update, h, δτ, cutoff, Dmax, op)
        _, magnet = magnetagainstsite(update, j, gammescale)
        meanvalues[p] = mean(magnet)
    end
    return sweeplist, meanvalues
end

"""
return the betalist and the energy list with a random uniform on each site
"""
function magnetforbestalistdisorder(betalist, ancilla, δτ, h, s, cutoff, gammescale, init, j::String)
    realbetalist = pushfirst!(diff(betalist), 0)
    Magnetlist = Vector{Vector{Float64}}(undef, length(realbetalist))
    update = ancilla
    _, gatesevolve, _ = evolutionwithrandomdisordergates(init::Int64, update, s, h, δτ)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ)
        _, Magnetlist[i] = magnetagainstsite(ancilla, j, gammescale)
    end
    return Magnetlist
end

"""
return the betalist and the energy list with a random uniform on each site
"""
function magnetforbestalist(betalist, ancilla, δτ, h, s, cutoff, gammescale, op, j::String)
    realbetalist = pushfirst!(diff(betalist), 0)
    Magnetlist = Vector{Vector{Float64}}(undef, length(realbetalist))
    update = ancilla
    gatesevolve = gatesTEBDancilla(update, h, δτ, s, op)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ)
        _, Magnetlist[i] = magnetagainstsite(ancilla, j, gammescale)
    end
    return Magnetlist
end
# ======================================== Correlation
"""
return the list of correlation function on the whole chain with the two boundaries excluded
"""
function correlationagainstsite(mps, j)
    N = length(mps)
    lengthlist = collect(1:1:N-2)
    Correlation = Vector{Float64}(undef, length(lengthlist))
    @showprogress desc = "correlation" for p in eachindex(lengthlist)
        Correlation[p] = correlationonlength(mps, lengthlist[p], j)
    end
    return lengthlist, Correlation
end




