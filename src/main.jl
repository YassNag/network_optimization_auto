########### Importations ###########

# Système
using Dates
using Random
using CSV

# Affichage graphique
using PyCall
using PyPlot
pygui(:qt5)
Path = PyPlot.matplotlib.path.Path

# Modélisation
using JuMP

# Solveurs libres
# using GLPK
# using Cbc
# using SCIP

# Solveurs commerciaux
#using Gurobi
using CPLEX
using DataFrames
# using MosekTools

#Graph library
using LightGraphs

# Divers
include("ascii_art.jl")

# Essentiels
include("structures.jl")
include("parser.jl")
include("utilitaires.jl")
include("visualisation.jl")
include("controles_saisies.jl")
include("parametres_solveurs.jl")

include("model.jl")


########### Programme principal ###########

function main()
    # Nettoyage de la console
    clear()
    ascii_art_pmcn()
    # Lecture des bibliothèques d'instances
    println("Saisissez le nom du répértoire des expérémentations souhaité: ")
    expefolder = readline()
    bibliotheques = lecture_bibliotheques(expefolder)
    connexity_module, ajout_R3, ajout_R4, callback_vide, ajout_R5, only_on_root, frequency, separation_exacte, solver, display_log, temps_limite =
    lecture_parametres_resolution(expefolder)
    # println("ajout_R3: $ajout_R3, ajout_R4: $ajout_R4, callback_vide: $callback_vide, ajout_R5: $ajout_R5, only_on_root: $only_on_root, frequency: $frequency, separation_exacte: $separation_exacte")
    df=DataFrame(Bib=String[], Inst=String[], Var=Int[], Const=Int[], DomConst=Int[], CoverI=Int[], Nodes=Int[], Gap=Float64[], CPU=Float64[])
    for bib in bibliotheques
        println("bib: $bib")
        liste_scenarios = readdir("bibliotheques_instances/"*bib*"/scenarios/")
        for scenario in liste_scenarios
            try
                ResultDict=Dict(:Bib=>"", :Inst=>"", :Var=>0, :Const=>0, :DomConst=>0, :CoverI=>0, :Nodes=>0, :Gap=>0, :CPU=>0)
                ResultDict[:Bib]=bib
                ResultDict[:Inst]=scenario
                println("scenario: $scenario")
                instance = lecture_instance(bib,scenario)
                solution = resolution_modele(instance, connexity_module, ajout_R3, ajout_R4, callback_vide, ajout_R5, only_on_root, frequency, separation_exacte, solver, display_log, temps_limite, ResultDict, expefolder)
                push!(df, ResultDict)
                CSV.write("Expes/"*expefolder*"/Results/IntResults.csv", df)
            catch
                break
            end
        end
    end
    CSV.write("Expes/"*expefolder*"/Results/results.csv", df)
    ascii_art_fin()
end

main()
