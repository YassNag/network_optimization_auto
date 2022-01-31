
########## Lecture d'une instance ##########

function lecture_instance(bib, scenario)
    # Lecture du contexte
    chemin_contexte = "bibliotheques_instances/"*bib*"/contexte.dat"
    contexte = lire_contexte(chemin_contexte)

    # Lecture du scénario
    chemin_scenario = "bibliotheques_instances/"*bib*"/scenarios/"*scenario
    scenario = lire_scenario(chemin_scenario, contexte)

    # Création de l'instance
    instance = Instance(contexte, scenario)

    return instance
end

function lecture_bibliotheques(expefolder)
    #lecture des nom de bibliothèques à résoudre
    lignes = Vector{String}
    open("Expes/"*expefolder*"/parametres_expes.txt") do fichier
        lignes = readlines(fichier)
    end
    bibliotheques = [bib for bib in split(lignes[2], ", ")] # nRc[c] : nombre de relais de type c
    return bibliotheques
end


function lecture_parametres_resolution(expefolder)
    lignes = Vector{String}
    open("Expes/"*expefolder*"/parametres_expes.txt") do fichier
        lignes = readlines(fichier)
    end
    connexity_module =  lignes[5] == "yes" ? true : false
    ajout_R3 = lignes[8] == "yes" ? true : false
    ajout_R4 = lignes[11] == "yes" ? true : false
    callback_vide  = lignes[14] == "yes" ? true : false
    ajout_R5  = lignes[17] == "yes" ? true : false
    only_on_root = lignes[20] == "yes" ? true : false
    frequency =  parse.(Int64, lignes[23])
    separation_exacte = lignes[26] == "yes" ? true : false
    solver = lignes[29]
    display_log = lignes[32]
    temps_limite =  parse.(Int64, lignes[35])

    return connexity_module, ajout_R3, ajout_R4, callback_vide, ajout_R5, only_on_root, frequency, separation_exacte, solver, display_log, temps_limite
end
########## Lecture (parsing) du contexte ##########

function lire_contexte(chemin)
    # Parsing
    lignes = Vector{String}
    open(chemin) do fichier
        lignes = readlines(fichier)
    end

    # Données élémentaires
    nT, nH, nE, nC, nU, nB = parse.(Int, split(lignes[2], ","))

    # Altitude
    A = parse(Int, lignes[5])

    # Relais
    nRc = [elt for elt in parse.(Int, split(lignes[8], ","))] # nRc[c] : nombre de relais de type c
    nR = sum(nRc)

    # Construction de l'"ensemble" des relais par type de communication
    # R[c] = [] ensemble des relais de type c
    R = Vector{Vector{Int}}(undef, nC)
    cpt_r = 1
    for c in 1:nC
        temp_r = Vector{Int}(undef, nRc[c])
        for r in 1:nRc[c]
            temp_r[r] = cpt_r
            cpt_r += 1
        end
        R[c] = temp_r
    end

    # Types de communications (par unité) : 1 unité par ligne
    # Cu[u] = [] ensemble des types de communication de l'unité u
    Cu = Vector{Vector{Int}}(undef, nU)
    ligne = 11
    for u in 1:nU
        Cu[u] = [elt for elt in parse.(Int, split(lignes[ligne], ","))]
        ligne += 1
    end

    # Portée de chaque unité pat Types de communications : 1 unité par ligne
    # Cu[u] = [] ensemble des types de communication de l'unité u
    Pu = Vector{Vector{Int}}(undef, nU)
    ligne += 2
    for u in 1:nU
        Pu[u] = [elt for elt in parse.(Int, split(lignes[ligne], ","))]
        ligne += 1
    end

    # Types de communications (par base) : 1 base par ligne
    # Cb[b] = [] ensemble des types de communication de la base b
    Cb = Vector{Vector{Int}}(undef, nB)
    ligne += 2
    for b in 1:nB
        Cb[b] = [elt for elt in parse.(Int, split(lignes[ligne], ","))]
        ligne += 1
    end

    # Positions des bases : 1 base par ligne
    # Eb[b] = (x,y) : position base b (point)
    Eb = Vector{Point}(undef, nB)
    ligne += 2
    for base in 1:nB
        x, y = parse.(Int, split(lignes[ligne], ","))
        point = Point(x, y)
        Eb[base] = point
        ligne += 1
    end

    # Positions envisageables pour le déploiement : 1 position par ligne
    # E[i] = (x,y) : position i (point)
    E = Vector{Point}(undef, nE)
    Ebool = Vector{Bool}(undef, nE) # Liste des positions actives | Ebool[i] = 0 : position inactive (non-couvrante)
    ligne += 2
    for i in 1:nE
        x, y = parse.(Int, split(lignes[ligne], ","))
        point = Point(x, y)
        E[i] = point
        Ebool[i] = true
        ligne += 1
    end

    # Poids et puissance de chacun des HAPS : 1 HAPS par ligne
    # W[h] : poids HAPS h
    # P[h] : puissance HAPS h
    W = Vector{Int}(undef, nH)
    P = Vector{Int}(undef, nH)
    ligne += 2
    for i in 1:nH
        W[i], P[i] = parse.(Int, split(lignes[ligne], ","))
        ligne += 1
    end

    # Poids, puissance et portée (seuil) de chacun des relais : 1 paragraphe par type de relais, une ligne par relais de chaque type
    # w[r] : poids relais r
    # p[r] : puissance relais r
    # s[r] : portée (seuil) relais r
    w = Vector{Int}(undef, nR)
    p = Vector{Int}(undef, nR)
    s = Vector{Int}(undef, nR)
    ligne += 2
    for c in 1:nC
        for r in R[c]
            w[r], p[r], s[r] = parse.(Int, split(lignes[ligne], ","))
            ligne += 1
        end
        ligne += 2
    end
#
# println("$Cu")
# println("$Pu")
# println("$Cb")
    contexte = Contexte(nT, nH, nE, nC, nU, nB, nR, nRc, A, R, Cu, Pu, Cb, Eb, E, Ebool, W, P, w, p, s)
    return contexte
end

########## Lecture (parsing) du scénario ##########

function lire_scenario(chemin, contexte)
    # Parsing
    lignes = Vector{String}
    open(chemin) do fichier
        lignes = readlines(fichier)
    end

    # Ensemble des déplacements des unités
    # deplacements[u] : déplacements de l'unité u
    # deplacements[u][t] : position de l'unité u au temps t
    deplacements = Vector{Vector{Point}}(undef, contexte.nU)
    ligne = 2
    for u in 1:contexte.nU
        deplacements[u] = Vector{Point}(undef, contexte.nT)
        temp = split(lignes[ligne], ")")
        temp = split.(temp[1:end-1], "(")
        popfirst!.(temp)
        for t in 1:contexte.nT
           coord = parse.(Float64, split(string(temp[t][1]), ","))
           deplacements[u][t] = Point(coord[1], coord[2])
        end
        ligne += 3
    end
    scenario = Scenario(deplacements)
    return scenario
end
