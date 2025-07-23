#!usr/bin/env julia
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
############### Librairies #################
using MBL
using MKL
using ProgressMeter
using JSON
using Statistics
using Dates
using LinearAlgebra
using Pkg
# ===================== log
println(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
println("Julia $VERSION")
@show Base.julia_cmd()
@show Threads.nthreads()
@show LinearAlgebra.BLAS.get_num_threads()
Pkg.status()
# ===================== parameters
if length(ARGS) < 1
    println("Missing input file: use default")
    throw(ErrorException)
else
    json_input = String(ARGS[1])
end

println("\nLoad input parameters from file $json_input")

input_data = JSON.parsefile(json_input)
map(k -> println(k, ": ", input_data[k]), sort(collect(keys(input_data))))

N = input_data["N"]
Ndisorder = input_data["Ndisorder"]
J = input_data["J"]
D0 = input_data["D0"]
betamax = input_data["maximum value for beta"]
stepbeta = input_data["step value for beta list"]
h0 = input_data["fixed disorder"]
h = input_data["maximum disorder magnitude"]
δτ = input_data["Trotter-Suzuki step"]
Dmax = input_data["max bond dimension"]
gammesweep = input_data["nsweep range"]
gammescale = input_data["gammescale"]
cutoff = input_data["cutoff"]
centerpic = input_data["phase transition point"]
initseed = input_data["fixed seed"]
n_sweep = input_data["fixed number of sweep"]
j = input_data["axis"]
init = input_data["initialization"]
savefile = String(input_data["savefile"])

################# Run ###############
sweep_list = collect(gammesweep[1]:gammesweep[3]:gammesweep[2])
realsweeplist = [gammesweep[3] for k in 1:floor(Int, ((gammesweep[2] - gammesweep[1]) / gammesweep[3]))+1]
realsweeplist[1] = gammesweep[1]
if init == "random"
    mps_random_debut, _ = random_initialized_MPS(N, D0)
elseif init == "neel"
    mps_random_debut, _ = neelstate(N)
end
# ===================== data 

metadata = Dict{String,Any}(
    "N" => N,
    "Trotter-Suzuki time step" => δτ,
    "initial bond dimension" => D0,
    "Dmax" => Dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "disorder fixed" => h0,
    "proportion spin average" => gammescale,
    "sweep range" => sweep_list,
    "effective sweep list" => realsweeplist,
    "maximum bond dimension per tebd step" => nothing,
    "type d'initialisation" => init
)
#println("\nmetadata:")
#display(metadata)

results = Dict{String,Any}(
    "energy sweep list" => nothing,
    "magnetization sweep list" => nothing,
)

Energytebd = Vector()
Magnettebd = Vector()
Maxbonddim = Vector()
#####evolv
function void()
    update_tebd = deepcopy(mps_random_debut)
    for i in eachindex(realsweeplist)
        println("Time evolution with tebd")
        update_tebd = tebdstepHeisenbergRow!(realsweeplist[i], update_tebd, h0, δτ, cutoff, Dmax)
        push!(Maxbonddim, maxbonddim(update_tebd))
        metadata["maximum bond dimension per tebd step"] = Maxbonddim

        #####measure
        println("Measure average energy")
        _, e = energyagainstsite(update_tebd, h0, gammescale)
        push!(Energytebd, mean(e))
        results["energy sweep list"] = Energytebd

        println("Mesure average magnet")
        _, magnet = magnetagainstsite(update_tebd, j, gammescale)
        push!(Magnettebd, mean(magnet))
        results["magnetization sweep list"] = Magnettebd

        #####data saving
        output_data_local = merge(metadata, results)
        #savefile = get_savefile(output_data)
        open(savefile, "w") do io
            JSON.print(io, output_data_local, 4)
        end
        println("\nResults saved in $savefile")
        flush(stdout)
    end
end

metadatadisorder = Dict{String,Any}(
    "N" => N,
    "N disorder" => Ndisorder,
    "Trotter-Suzuki time step" => δτ,
    "Dmax" => Dmax,
    "J" => J,
    "axis spin" => j,
    "cutoff" => cutoff,
    "betamax" => betamax,
    "seed" => initseed,
    "step beta" => stepbeta,
    "maximum disorder magnitude" => h,
    "seed" => initseed,
    "phase transition point" => centerpic,
    "proportion spin average" => gammescale,
    "maximum bond dimension per tebd step" => nothing,
)
println("\nmetadatadisorder:")
display(metadatadisorder)

ancilla, s = MBL.AncillaMPO(N)
disorder = sort(MBL.rejection_sample(Ndisorder, h, centerpic, initseed; σ=0.03h, A=1.0))
st, dp = MBL.section_trunc(N, gammescale)
L = collect(st:dp)
# =========================
betalist = collect(0:stepbeta:betamax)
# =========================

output_data = Dict{String,Any}(
    "energy" => nothing,
    "magnet" => nothing,
    "beta list" => betalist,
    "champ" => disorder,
    "site list" => L,
    "mean value list" => nothing,
    "std list" => nothing
)

function voidmeanandstd()
    #pour chaque champ, pour chaque beta : mesurer l'énergie de chaque site, moyenne et ecart type 

    Energylist = Array{Float64}(undef, length(L), length(betalist), length(disorder))
    Magnetlist = Array{Float64}(undef, length(L), length(betalist), length(disorder))

    # ========================= SIMU 

    @showprogress desc = "run over disorder" for i in eachindex(disorder)
        println("h = ", disorder[i])
        valuee, valuem = MBL.magnetandenergyforbetalistdisorder(betalist, ancilla, δτ, disorder[i], s, cutoff, gammescale, initseed, Dmax, j)
        Magnetlist[:, :, i] = valuem
        output_data["magnet"] = Magnetlist
        printl("Sz part done")
        Energylist[:, :, i] = valuee
        output_data["energy"] = Energylist
        printl("Energy part done")

        #####data saving
        final_data = merge(metadatadisorder, output_data)
        #savefile = get_savefile(output_data)
        open(savefile, "w") do io
            JSON.print(io, final_data, 4)
        end
        println("\nResults saved in $savefile")
        flush(stdout)
    end
    Meanvalue = mean(Magnetlist; dims=1)
    Stdvalue = std(Magnetlist; dims=1)
    output_data["mean value list"] = Meanvalue
    output_data["std value list"] = Stdvalue
    final_data = merge(metadatadisorder, output_data)
    #savefile = get_savefile(output_data)
    open(savefile, "w") do io
        JSON.print(io, final_data, 4)
    end
    println("\nResults saved in $savefile")
    flush(stdout)
end

voidmeanandstd()

println("simulation finie")
