#!usr/bin/env julia

##################################### Librairies #############

using ITensors
using MBL
using ProgressMeter
using Statistics

##################################### Commentaires #############

#=
on pourrait refactoriser encore le code en evitant de regénérer les gates à chaque fois qu'on mesure
comme on fait quand on ajoute le champ aléatoire mais générer les gates est o devant l'évolution et 
flemme de réécrire pour l'instant
=#

##################################### Functions #############

# ============================================= Energy
"""
mps_init_sweep -- initial boundary mps 
gammesweep -- start and end of the number of TEBD sweep with the step [begin:end:step]
gammescale -- 0<x<1 in order to compute the operator on gammescale*length(mps_init_sweep) from the center
cutoff -- cutoff in the svd 
dmax -- maximal bond dimension
δτ -- Trotter Suzuki step 
h -- disorder 
op -- choose between Heisenberg Hamiltonian and XY model

return the average energy on x% of the total number of spins with respect to number of spins
"""
function energyaverageagainstsweep(
    mps_init_sweep, gammesweep, gammescale, cutoff, dmax, δτ, h, op::String
)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2]) #on crée la liste des sweeps (nombre total de sweep de chaque pas)
    realsweeplist = [
        gammesweep[3] for k in 1:(floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3])) + 1)
    ] #on est plus efficace si on garde le même mps et qu'on ajoute des sweep
    realsweeplist[1] = gammesweep[1]
    meanvalues = Vector(undef, length(sweeplist))
    update = mps_init_sweep
    gates = gateTrotterSuzukirow(mps_init_sweep, h, δτ, op)
    for p in eachindex(realsweeplist)
        update = tebdevolutionrow!(realsweeplist[p], update, gates, cutoff, dmax)
        _, magnet = energyagainstsite(update, h, gammescale, op)
        meanvalues[p] = mean(magnet)
    end
    return sweeplist, meanvalues
end

"""
gammelength -- range of length of the mps ([begin:end:step])
numbersweep -- fixed number of TEBD sweeps
gammescale -- 0<x<1 in order to compute the operator on gammescale*length(mps_init_sweep) from the center
cutoff -- cutoff in the svd 
dmax -- maximal bond dimension
δτ -- Trotter Suzuki step 
h -- disorder 
op -- choose between Heisenberg Hamiltonian and XY model

return the average energy (over gammescale*length spins) for a fixed number of tebd steps but with different length of mps
"""
function energyaverageagainstlength(
    gammelength::Tuple, gammescale, numbersweep, cutoff, dmax, D0, δτ, h, op::String
)
    sites = collect(gammelength[1]:gammelength[3]:gammelength[2])
    averageenergy = Vector(undef, length(sites))
    @showprogress for i in eachindex(sites)
        mpstransit, _ = random_initialized_MPS(sites[i], D0)
        gates = gateTrotterSuzukirow(mpstransit, h, δτ, op)
        converged = tebdevolutionrow!(numbersweep, mpstransit, gates, cutoff, dmax)
        _, averageenergytemp = energyagainstsite(converged, h, gammescale, op)
        averageenergy[i] = mean(averageenergytemp)
    end
    return sites, averageenergy
end

##################################### Tracer Finite Temperature #############
"""
betamax -- maximal beta you want to reach
step -- step in the beta list
ancilla -- MPS

return the MPS average energy measured with gates on gammescale*length(MPS) number of sites taken from the MPS center
"""
function energyforbetalist(betalist, ancilla::MPO, δτ, h, s, cutoff, op::String, gammescale, dmax)
    realbetalist = pushfirst!(diff(betalist), 0) #realbetalist and betalist same length
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Energylist = fill(1.0,  dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gates = gatesTEBDancilla(update, h, δτ, s, op)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = TEBDancilla!(update, gates, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Energylist[:, i] = energyagainstsite(update, h, gammescale, op)
        println("Average energy at β=$(betalist[i]): ", mean(Energylist; dims=1)[i])
    end
    return Energylist
end

"""
betamax -- maximal beta you want to reach
step -- step in the beta list
ancilla -- MPS

return the MPS average energy measured with gates on gammescale*length(MPS) number of sites taken from the MPS center
"""
function energyforbetalist(betalist, mps::MPS, δτ, h, cutoff, op::String, gammescale, dmax)
    realbetalist = pushfirst!(diff(betalist), 0) #realbetalist and betalist same length
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Energylist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gatesevolve = gateTrotterSuzukirow(mps, h, δτ, op)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        nsweep = floor(realbetalist[i] / 2 / δτ)
        update = tebdevolutionrow!(nsweep, update, gatesevolve, cutoff, dmax)
        _, Energylist[:, i] = energyagainstsite(update, h, gammescale, op)
        println("Average energy at β=$betalist[i] = ", mean(Energylist; dims=1))
    end
    return Energylist
end

"""
return the energy list of the site i with respect to gates time step
"""
function energyagainstdeltatime(
    site_measure, gamme::Tuple, mpsinit, step, numbersweep, cutoff, dmax, op::String, h
)
    timesteplist = reverse(collect(gamme[1]:step:gamme[2]))
    EnergyList = Vector(undef, length(timesteplist))
    update = mpsinit
    @showprogress for j in eachindex(timesteplist)
        gates = gateTrotterSuzukirow(update, h, timesteplist[i], op)
        update = tebdevolutionrow!(numbersweep, update, gates, cutoff, dmax)
        e = energysite(mpsinit, site_measure, h, op)
        EnergyList[j] = e
    end
    return timesteplist, EnergyList
end

##################################### Random disorder #############

"""
mps -- MPS
h -- disorder 
scale -- 0<scale<1, pourcentage of the chain you want to measure (from the middle chain)

return the MPS average energy for the model op, measured with gates 
"""
function energyagainstsiteMPOdisorder(mps, gates, scale)
    N = length(mps)
    start, stop = MBL.section_trunc(N, scale)
    stop = stop < N - 2 ? stop : N - 2
    sites = collect(start:1:stop)
    Energypersite = Array{Float64}(undef, length(sites))
    update = mps
    for i in eachindex(sites)
        Energypersite[i] = energysitedisorder(update, sites[i], gates[sites[i]])
    end
    return sites, Energypersite
end

"""
return the energy list with a random uniform on each site
"""
function energyforbestalistdisorder(
    betalist, ancilla, δτ, h, s, cutoff, gammescale, init::Int64, dmax
)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Energylist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gatesmeasure, gatesevolve, _ = MBL.evolutionwithrandomdisordergates(init, update, s, h, δτ)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Energylist[:, i] = energyagainstsiteMPOdisorder(update, gatesmeasure, gammescale)
        println("Average energy at β=$betalist[i] = ", mean(Energylist; dims=1))
    end
    return Energylist
end

function magnetandenergyforbetalistdisorder(
    betalist, ancilla, δτ, h::Float64, s, cutoff, gammescale, init, dmax, j::String
)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Energylist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    Magnetlist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    Dimensionlist = Array{Float64}(undef, length(realbetalist))
    update = ancilla
    gatesmeasure, gatesevolve, Disorder = MBL.evolutionwithrandomdisordergates(
        init, update, s, h, δτ
    )
    for i in eachindex(realbetalist)
        println("\n# " * "="^30)
        beta = betalist[i]
        @info "β[$i]" beta
        update = TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Energylist[:, i] = MBL.energyagainstsiteMPOdisorder(update, gatesmeasure, gammescale)
        println("Average energy at β=$beta : ", mean(Energylist; dims=1)[i])
        #println("\n# " * "="^30)
        _, Magnetlist[:, i] = magnetagainstsite(update, j, gammescale)
        println("Average Sz at β=$beta : ", mean(Magnetlist; dims=1)[i])
        #println("\n# " * "="^30)
        Dimensionlist[i] = maxbonddim(update)
        println("Maximum bond dimension at β=$beta : ", Dimensionlist[i])
        flush(stdout)
    end
    return Energylist, Magnetlist, Disorder, Dimensionlist
end

function magnetandenergyforbetalistdisorder(
    betalist, ancilla, δτ, h::Vector, s, cutoff, gammescale, dmax, j::String
)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Energylist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    Magnetlist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    Dimensionlist = Array{Float64}(undef, length(realbetalist))
    update = ancilla
    gatesmeasure, gatesevolve = MBL.evolutionwithrandomdisordergates(update, s, h, δτ)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDanxcilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Energylist[:, i] = MBL.energyagainstsiteMPOdisorder(update, gatesmeasure, gammescale)
        println("Average energy at β=$betalist[i] = ", mean(Energylist; dims=1))
        _, Magnetlist[:, i] = MBL.magnetagainstsite(update, j, gammescale)
        println("Average Sz at β=$betalist[i] = ", mean(Magnetlist; dims=1))
        Dimensionlist[i] = maxbonddim(update)
        println("Maximum bond dimension at β=$betalist[i] = ", Dimensionlist[i])
    end
    return Energylist, Magnetlist, Dimensionlist
end

function magnetandenergyforbetalist(
    betalist, ancilla, δτ, h, s, cutoff, gammescale, dmax, j::String, op::String
)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Energylist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    Dimensionlist = Array{Float64}(undef, length(realbetalist))
    Magnetlist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gates = gatesTEBDancilla(ancilla, h, δτ, s, op)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = TEBDancilla!(update, gates, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Energylist[:, i] = MBL.energyagainstsite(update, h, gammescale, op)
        println("Average energy at β=$betalist[i] = ", mean(Energylist; dims=1))
        _, Magnetlist[:, i] = MBL.magnetagainstsite(update, j, gammescale)
        println("Average Sz at β=$betalist[i] = ", mean(Magnetlist; dims=1))
        Dimensionlist[i] = maxbonddim(update)
        println("Maximum bond dimension at β=$betalist[i] = ", Dimensionlist[i])
    end
    return Energylist, Magnetlist, Dimensionlist
end

"""
betalist -- list des beta 
ancilla -- boudary mps
δτ -- step de Trotter Suzuki
h -- disorder
s -- site index
cutoff -- cutoff de la svd
sitemeasure -- site sur lequel on mesure
init -- valeur de la seed

return la liste des énergies sur les sites de mesures pour une liste de valeur de beta
"""
function trackonesite(betalist, ancilla, δτ, h, s, cutoff, sitemeasure, init::Int, dmax)
    realbetalist = pushfirst!(diff(betalist), 0)
    Energylist = Vector{}(undef, length(realbetalist))
    update = ancilla
    gatesmeasure, gatesevolve, _ = MBL.evolutionwithrandomdisordergates(init, update, s, h, δτ)
    @showprogress desc = "compute energy for β" for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)
        Energylist[i] = energysitedisorder(update, sitemeasure, gatesmeasure[sitemeasure])
    end
    return Energylist
end

# ============================================= Magnetization
"""
mps
j::String
scale

return the site list and the magnet (Sz) per site
"""
function magnetagainstsite(mps, j::String, scale)
    N = length(mps)
    update = mps
    start, stop = section_trunc(N, scale)
    sites = collect(start:1:stop)
    Magnetpersite = Array{Float64}(undef, length(sites))
    for p in eachindex(sites)
        Magnetpersite[p] = measure_S(update, p, j)
    end
    return sites, Magnetpersite
end

"""
return the average spin  (over gammescale*length spins) against j axis for mps of length in the specif range 
"""
function magnetaverageagainstlength(
    j::String, gammelength::Tuple, gammescale, numbersweep, cutoff, dmax, D0, δτ, h, op::String
)
    sites = collect(gammelength[1]:1:gammelength[2])
    averagespin = Vector(undef, length(sites))
    for i in eachindex(sites)
        @show i, gammelength[2], sites[i]
        mpstransit, _ = random_initialized_MPS(sites[i], D0)
        gates = gateTrotterSuzukirow(mpstransit, h, δτ, op)
        converged = tebdevolutionrow!(numbersweep, mpstransit, gates, cutoff, dmax)
        _, averagespintemp = magnetagainstsite(converged, j, gammescale)
        averagespin[i] = mean(averagespintemp)
    end
    return sites, averagespin
end

"""
return the average spin (over gammescale*length spins) value gainst the axis j for a mps  with a fixed length but updated through tebd algorithm with a number of sweep in the gammesweep range
"""
function magnetaverageagainstsweep(
    j::String, mps_init_sweep, gammesweep, gammescale, cutoff, dmax, δτ, h, op::String
)
    sweeplist = collect(gammesweep[1]:gammesweep[3]:gammesweep[2])
    realsweeplist = [
        gammesweep[3] for k in 1:(floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3])) + 1)
    ]
    realsweeplist[1] = gammesweep[1]
    meanvalues = Vector(undef, length(realsweeplist))
    update = mps_init_sweep
    for p in eachindex(realsweeplist)
        gates = gateTrotterSuzukirow(update, h, δτ, op)
        update = tebdevolutionrow!(realsweeplist[p], update, gates, cutoff, dmax)
        _, magnet = magnetagainstsite(update, j, gammescale)
        meanvalues[p] = mean(magnet)
    end
    return sweeplist, meanvalues
end

"""
return the betalist and the energy list with a random uniform on each site
"""
function magnetforbestalistdisorder(
    betalist, ancilla, δτ, h, s, cutoff, gammescale, init::Int64, j::String, dmax
)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Magnetlist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    _, gatesevolve, disorder = MBL.evolutionwithrandomdisordergates(init, update, s, h, δτ)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Magnetlist[:, i] = magnetagainstsite(update, j, gammescale)
        println("Average Sz at β=$betalist[i] = ", mean(Magnetlist; dims=1))
    end
    return Magnetlist, disorder
end

"""

"""
function magnetforbestalist(
    betalist, ancilla::MPO, δτ, h, s, cutoff, gammescale, op, j::String, dmax
)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Magnetlist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gatesevolve = gatesTEBDancilla(update, h, δτ, s, op)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        update = MBL.TEBDancilla!(update, gatesevolve, realbetalist[i] / 2, cutoff, δτ, dmax)
        _, Magnetlist[:, i] = magnetagainstsite(update, j, gammescale)
        println("Average Sz at β=$betalist[i] = ", mean(Magnetlist; dims=1))
    end
    return Magnetlist
end

"""

"""
function magnetforbestalist(betalist, ancilla::MPS, δτ, h, cutoff, gammescale, op, j::String, dmax)
    realbetalist = pushfirst!(diff(betalist), 0)
    N = length(ancilla)
    st, dp = MBL.section_trunc(N, gammescale)
    Magnetlist = Array{Float64}(undef, dp - st + 1, length(realbetalist)) #fist dimension for the spin chain, second dimension for the beta list
    update = ancilla
    gatesevolve = gateTrotterSuzukirow(update, h, δτ, op)
    for i in eachindex(realbetalist)
        @info "β[$i]" betalist[i]
        nsweep = floor(realbetalist[i] / 2 / δτ)
        update = tebdevolutionrow!(nsweep, update, gatesevolve, cutoff, dmax)
        _, Magnetlist[:, i] = magnetagainstsite(update, j, gammescale)
        println("Average Sz at β=$betalist[i] = ", mean(Magnetlist; dims=1))
    end
    return Magnetlist
end

# ======================================== Correlation
"""
return the list of correlation function on the whole chain with the two boundaries excluded
"""
function correlationagainstsite(mps, j)
    N = length(mps)
    lengthlist = collect(1:1:(N - 2))
    Correlation = Vector{Float64}(undef, length(lengthlist))
    @showprogress desc = "correlation" for p in eachindex(lengthlist)
        Correlation[p] = correlationonlength(mps, lengthlist[p], j)
    end
    return lengthlist, Correlation
end
