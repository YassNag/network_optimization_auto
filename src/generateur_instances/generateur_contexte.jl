########### Importations ###########

# Dépendances
using Dates
using Random

# Essentiels
include("controles_saisies.jl")


########## Générateur de contexte ##########

function generation_contexte()
    # Nettoyage de la console
    clear()

    println("\n========== Générateur de contexte ==========")
    contexte = ""
    doc = "" # Documentation de l'instance en cours de génération

    doc *= "========== Documentation =========="

    docTemp = "\n --> Nombre d'états : "
    print(docTemp)
    nT = readline()
    doc *= "\n" * docTemp * nT
    contexte *= "// Nombre d'états (#T), de HAPS (#H), de positions envisageables pour le déploiement (#E) de types de communication (#C), d'unités (#U) et de bases (#B)\n"
    contexte *= nT * ", "

    docTemp = "\n --> Nombre de HAPS : "
    print(docTemp)
    nH = readline()
    doc *= "\n" * docTemp * nH
    contexte *= nH * ", "
    nH = parse(Int, nH)

    choix_quadrillage = choix_binaire("\n --> Souhaitez-vous générer un quadrillage pour les positions de déploiement (o/n) ? ")
    if choix_quadrillage == "o"
        println("\n    - Quadrillage (1000 * 1000)")
        print("       --> Pas : ")
        pas = parse(Int, readline())
        X = [0:pas:1000;]; Y = [0:pas:1000;]
        E = [(x,y) for x in X, y in Y]
        nE = length(E)
        docTemp = "\n --> Nombre de positions envisageables pour le déploiement d'un HAPS : $nE (quadrillage, pas : $pas)"
        println(docTemp)
        doc *= "\n" * docTemp
        contexte *= string(nE) * ", "
    elseif choix_quadrillage == "n"
        docTemp = "\n --> Nombre de positions envisageables pour le déploiement d'un HAPS : "
        print(docTemp)
        nE = readline()
        doc *= "\n" * docTemp * nE
        contexte *= nE * ", "
        nE = parse(Int, nE)
    end

    docTemp = "\n --> Nombre de types de communication : "
    print(docTemp)
    nC = readline()
    doc *= "\n" * docTemp * nC
    contexte *= nC * ", "
    nC = parse(Int, nC)

    docTemp = "\n --> Nombre d'unités : "
    print(docTemp)
    nU = readline()
    doc *= "\n" * docTemp * nU
    contexte *= nU * ", "
    nU = parse(Int, nU)

    docTemp = "\n --> Nombre de bases : "
    print(docTemp)
    nB = readline()
    doc *= "\n" * docTemp * nB
    contexte *= nB * "\n\n"
    nB = parse(Int, nB)

    contexte *= "// Altitude de déploiement des HAPS (A)\n"
    docTemp = "\n --> Altitude de déploiement des HAPS : "
    print(docTemp)
    A = readline()
    doc *= "\n" * docTemp * A
    contexte *= A * "\n\n"

    nRc = []
    contexte *= "// Nombre de relais par type de communication (R_1, R_2, ...)\n"
    choix, a, b = choix_saisie("\n --> Souhaitez-vous saisir manuellement le nombre de relais par type de communication (o/n) ? ")
    for c in 1:nC
        docTemp = "\n --> Nombre de relais de type $c : "
        print(docTemp)
        doc *= "\n" * docTemp
        if choix == "o"
            val = readline()
            doc *= val
        elseif choix == "n"
            val = rand(a:b)
            println("$val")
            val = string(val)
            doc *= val
        end
        if c != nC
            contexte *= val * ", "
        else
            contexte *= val * "\n"
        end
        push!(nRc, parse(Int, val))
    end
    contexte *= "\n"

    contexte *= "// Types de communication de chaque unité (C_u)\n"
    choix, a, b = choix_saisie("\n --> Souhaitez-vous saisir manuellement les types de communication de chaque unité (o/n) ? ")
    CU = Vector{Vector{Int}}(undef, nU)
    for u in 1:nU
        CU[u]=Int[]
        if choix == "o"
            docTemp = "\n - Unité $u : "
            println(docTemp)
            doc *= "\n" * docTemp
            docTemp = "   --> Nombre de types de communication : "
            print(docTemp)
            doc *= "\n" * docTemp
            nCu = parse(Int, readline())
            docTemp *= string(nCu)
            reste = [1:nCu;]
            for c in 1:nCu
                println("      Types restants : $reste")
                docTemp = "      --> Type $c : "
                print(docTemp)
                doc *= "\n" * docTemp
                Cu_i = readline()
                filter!(el->el!=parse(Int, Cu_i),reste)
                doc *= Cu_i
                if c != nCu
                    contexte *= Cu_i * ", "
                    push!(CU[u], Cu_i)
                else
                    contexte *= Cu_i * "\n"
                    push!(CU[u], Cu_i)
                end
            end
        elseif choix == "n"
            docTemp = "\n --> Types de communication de l'unité $u : "
            print(docTemp)
            doc *= "\n" * docTemp
            nCu = rand([a:b;])
            Cu = shuffle(1:nC)[1:nCu]
            cpt = 1
            for Cu_i in Cu
                if cpt != nCu
                    print("$(Cu_i), ")
                    doc *= string(Cu_i) * ", "
                    contexte *= string(Cu_i) * ", "
                    push!(CU[u], Cu_i)
                else
                    println("$(Cu_i)")
                    doc *= string(Cu_i)
                    contexte *= string(Cu_i) * "\n"
                    push!(CU[u], Cu_i)
                end
                cpt += 1
            end
        end
    end

# println(CU)
    contexte *= "\n"

    contexte *= "// Portée de chaque unité par type de communication (P_u)\n"

    choix_portee, a_portee, b_portee = choix_saisie("\n --> Souhaitez-vous saisir manuellement la portée de communication de chaque unité (o/n) ? ")
    for u in 1:nU
        docTemp = "      --> Portée unité (seuil) : "
        print(docTemp)
        doc *= "\n" * docTemp
        cpt=1
        for c in CU[u]
            if choix_portee == "o"
                s = readline()
                doc *= s
            elseif choix_portee == "n"
                s = rand(a_portee:b_portee)
                println("$s")
                s = string(s)
                doc *= s
            end
            if cpt != length(CU[u])
                contexte *= s * ", "
            else
                contexte *= s * "\n"
            end
            cpt+=1
        end
    end

    # for u in 1:nU
    #     for c in
    #     if choix == "o"
    #         docTemp = "\n - Unité $u : "
    #         println(docTemp)
    #         doc *= "\n" * docTemp
    #         docTemp = "   --> Nombre de types de communication : "
    #         print(docTemp)
    #         doc *= "\n" * docTemp
    #         nCu = parse(Int, readline())
    #         docTemp *= string(nCu)
    #         reste = [1:nCu;]
    #         for c in 1:nCu
    #             println("      Types restants : $reste")
    #             docTemp = "      --> Type $c : "
    #             print(docTemp)
    #             doc *= "\n" * docTemp
    #             Cu_i = readline()
    #             filter!(el->el!=parse(Int, Cu_i),reste)
    #             doc *= Cu_i
    #             if c != nCu
    #                 contexte *= Cu_i * ", "
    #             else
    #                 contexte *= Cu_i * "\n"
    #             end
    #         end
    #     elseif choix == "n"
    #         docTemp = "\n --> Types de communication de l'unité $u : "
    #         print(docTemp)
    #         doc *= "\n" * docTemp
    #         nCu = rand([a:b;])
    #         Cu = shuffle(1:nC)[1:nCu]
    #         cpt = 1
    #         for Cu_i in Cu
    #             if cpt != nCu
    #                 print("$(Cu_i), ")
    #                 doc *= string(Cu_i) * ", "
    #                 contexte *= string(Cu_i) * ", "
    #             else
    #                 println("$(Cu_i)")
    #                 doc *= string(Cu_i)
    #                 contexte *= string(Cu_i) * "\n"
    #             end
    #             cpt += 1
    #         end
    #     end
    # end
    #

    contexte *= "\n"

    contexte *= "// Types de communication de chaque base (C_b)\n"
    choix, a, b = choix_saisie("\n --> Souhaitez-vous saisir manuellement les types de communication de chaque base (o/n) ? ")
    for base in 1:nB
        if choix == "o"
            docTemp = "\n - Base $base : "
            println(docTemp)
            doc *= "\n" * docTemp
            docTemp = "   --> Nombre de types de communication : "
            print(docTemp)
            doc *= "\n" * docTemp
            nCb = parse(Int, readline())
            docTemp *= string(nCb)
            reste = [1:nCb;]
            for c in 1:nCb
                println("      Types restants : $reste")
                docTemp = "      --> Type $c : "
                print(docTemp)
                doc *= "\n" * docTemp
                Cb_i = readline()
                filter!(el->el!=parse(Int, Cb_i),reste)
                doc *= Cb_i
                if c != nCb
                    contexte *= Cb_i * ", "
                else
                    contexte *= Cb_i * "\n"
                end
            end
        elseif choix == "n"
            docTemp = "\n --> Types de communication de la base $base : "
            print(docTemp)
            doc *= "\n" * docTemp
            nCb = rand([a:b;])
            Cb = shuffle(1:nC)[1:nCb]
            cpt = 1
            for Cb_i in Cb
                if cpt != nCb
                    print("$(Cb_i), ")
                    doc *= string(Cb_i) * ", "
                    contexte *= string(Cb_i) * ", "
                else
                    println("$(Cb_i)")
                    doc *= string(Cb_i)
                    contexte *= string(Cb_i) * "\n"
                end
                cpt += 1
            end
        end
    end
    contexte *= "\n"

    contexte *= "// Positions des bases (x_b, y_b)\n"
    choix, a, b = choix_saisie("\n --> Souhaitez-vous saisir manuellement les positions des bases (o/n) ? ")
    for base in 1:nB
        docTemp = "\n - Base $base : "
        println(docTemp)
        doc *= "\n" * docTemp
        docTemp = "   --> x : "
        print(docTemp)
        doc *= "\n" * docTemp
        if choix == "o"
            x = readline()
            doc *= x
        elseif choix == "n"
            x = rand(a:b)
            println("$x")
            x = string(x)
            doc *= x
        end
        docTemp = "   --> y : "
        print(docTemp)
        doc *= "\n" * docTemp
        if choix == "o"
            y = readline()
            doc *= y
        elseif choix == "n"
            y = rand(a:b)
            println("$y")
            y = string(y)
            doc *= y
        end
        contexte *= x * ", " * y * "\n"
    end
    contexte *= "\n"

    contexte *= "// Positions envisageables pour le déploiement d'un HAPS (x_e, y_e)\n"
    if choix_quadrillage == "o"
        cpt = 1
        for (x,y) in E
            docTemp = "\n - Position $cpt : "
            doc *= "\n" * docTemp
            docTemp = "   --> x : $x"
            doc *= "\n" * docTemp
            docTemp = "   --> y : $y"
            doc *= "\n" * docTemp
            contexte *= string(x) * ", " * string(y) * "\n"
            cpt += 1
        end
        contexte *= "\n"
    elseif choix_quadrillage == "n"
        choix, a, b = choix_saisie("\n --> Souhaitez-vous saisir manuellement les positions enviseables pour le déploiement des HAPS (o/n) ? ")
        for position in 1:nE
            docTemp = "\n - Position $position : "
            println(docTemp)
            doc *= "\n" * docTemp
            docTemp = "   --> x : "
            print(docTemp)
            doc *= "\n" * docTemp
            if choix == "o"
                x = readline()
                doc *= x
            elseif choix == "n"
                x = rand(a:b)
                println("$x")
                x = string(x)
                doc *= x
            end
            docTemp = "   --> y : "
            print(docTemp)
            doc *= "\n" * docTemp
            if choix == "o"
                y = readline()
                doc *= y
            elseif choix == "n"
                y = rand(a:b)
                println("$y")
                y = string(y)
                doc *= y
            end
            contexte *= x * ", " * y * "\n"
        end
        contexte *= "\n"
    end

    contexte *= "// Poids et puissance maximale de chacun des HAPS (W_h, P_h)\n"
    choix_poids, a_poids, b_poids = choix_saisie("\n --> Souhaitez-vous saisir manuellement les poids maximales de chacun des HAPS (o/n) ? ")
    choix_puissance, a_puissance, b_puissance = choix_saisie("\n --> Souhaitez-vous saisir manuellement les puissances maximales de chacun des HAPS (o/n) ? ")
    for h in 1:nH
        docTemp = "\n - HAPS $h : "
        println(docTemp)
        doc *= "\n" * docTemp
        docTemp = "   --> Poids maximal : "
        print(docTemp)
        doc *= "\n" * docTemp
        if choix_poids == "o"
            W = readline()
            doc *= W
        elseif choix_poids == "n"
            W = rand(a_poids:b_poids)
            println("$W")
            W = string(W)
            doc *= W
        end
        docTemp = "   --> Puissance maximale : "
        print(docTemp)
        doc *= "\n" * docTemp
        if choix_puissance == "o"
            P = readline()
            doc *= P
        elseif choix_puissance == "n"
            P = rand(a_puissance:b_puissance)
            println("$P")
            P = string(P)
            doc *= P
        end
        contexte *= W * ", " * P * "\n"
    end
    contexte *= "\n"

    choix_poids, a_poids, b_poids = choix_saisie("\n --> Souhaitez-vous saisir manuellement les poids de chacun des relais (o/n) ? ")
    choix_puissance, a_puissance, b_puissance = choix_saisie("\n --> Souhaitez-vous saisir manuellement les puissances de chacun des relais (o/n) ? ")
    choix_portee, a_portee, b_portee = choix_saisie("\n --> Souhaitez-vous saisir manuellement les portées (seuils) de chacun des relais (o/n) ? ")
    cpt = 1
    for c in 1:nC
        contexte *= "// Poids, puissance et portée de chacun des relais de type $c (w_r, p_r, s_r)\n"
        docTemp = "\n - Relais de type $c"
        println(docTemp)
        doc *= "\n" * docTemp
        for r in 1:nRc[c]
            docTemp = "\n   - Relais $cpt : "
            println(docTemp)
            doc *= "\n" * docTemp
            docTemp = "      --> Poids : "
            print(docTemp)
            doc *= "\n" * docTemp
            if choix_poids == "o"
                w = readline()
                doc *= w
            elseif choix_poids == "n"
                w = rand(a_poids:b_poids)
                println("$w")
                w = string(w)
                doc *= w
            end
            docTemp = "      --> Puissance : "
            print(docTemp)
            doc *= "\n" * docTemp
            if choix_puissance == "o"
                p = readline()
                doc *= p
            elseif choix_puissance == "n"
                p = rand(a_puissance:b_puissance)
                println("$p")
                p = string(p)
                doc *= p
            end
            docTemp = "      --> Portée (seuil) : "
            print(docTemp)
            doc *= "\n" * docTemp
            if choix_portee == "o"
                s = readline()
                doc *= p
            elseif choix_portee == "n"
                s = rand(a_portee:b_portee)
                println("$s")
                s = string(s)
                doc *= s
            end
            cpt += 1
            contexte *= w * ", " * p * ", " * s * "\n"
        end
        if c != nC
            contexte *= "\n"
        end
    end
    doc *= "\n\n========== Documentation ==========\n"

    nom_dossier = string(today())*"_"*Dates.format(now(), "HH:MM:SS")
    run(`mkdir contextes/$(nom_dossier)`)
    open("contextes/$(nom_dossier)/contexte.dat", "w") do fichier
        write(fichier, contexte)
    end

    open("contextes/$(nom_dossier)/documentation", "w") do fichier
        write(fichier, doc)
    end

    println("\n /!\\ Le contexte est disponible dans le dossier \"contextes/$(nom_dossier)/\" /!\\ ")
    println("\n========== Générateur de contexte ==========")
end

generation_contexte()
