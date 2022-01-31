########## Affichage du contexte de l'instance courante ##########

function affichage_contexte(instance)
    println("\n========== Affichage du contexte ==========\n")
    
    println("    --> Nombre d'états (#T) : ", instance.contexte.nT)

    println("\n    --> Nombre de HAPS (#H) : ", instance.contexte.nH)

    println("\n    --> Nombre de positions envisageables pour déploiement (#E) : ", instance.contexte.nE)

    println("\n    --> Nombre de types de communication (#C) : ", instance.contexte.nC)

    println("\n    --> Nombre d'unités terrestres (#U) : ", instance.contexte.nU)

    println("\n    --> Nombre de bases (#B) : ", instance.contexte.nB)

    println("\n    --> Altitude de déploiement (A) : ", instance.contexte.A)

    println("\n    --> Nombre de relais (#R) : ", instance.contexte.nR)

    choix = choix_binaire("\n    --> Souhaitez-vous afficher la répartition des relais par type de communication (o/n) ? ")
    if choix == "o"
        repartition_relais(instance)
    end

    choix = choix_binaire("\n    --> Souhaitez-vous afficher les types de communication de chaque unité (o/n) ? ")
    if choix == "o"
        types_communication_unites(instance)
    end

    choix = choix_binaire("\n    --> Souhaitez-vous afficher les types de communication de chaque base (o/n) ? ")
    if choix == "o"
        types_communication_bases(instance)
    end

    choix = choix_binaire("\n    --> Souhaitez-vous afficher les positions des bases (o/n) ? ")
    if choix == "o"
        liste_positions_bases(instance)
    end

    choix = choix_binaire("\n    --> Souhaitez-vous afficher les différentes positions de déploiement (o/n) ? ")
    if choix == "o"
        liste_positions_HAPS(instance)
    end

    choix = choix_binaire("\n    --> Souhaitez-vous afficher le poids et la puissance maximale de chacun des HAPS (o/n) ? ")
    if choix == "o"
        liste_poids_puissances_HAPS(instance)
    end


    choix = choix_binaire("\n    --> Souhaitez-vous afficher le poids, la puissance et la portée (seuil) de chacun des relais (o/n) ? ")
    if choix == "o"
        liste_poids_puissances_portee_relais(instance)
    end

    println("\n========== Affichage du contexte ==========")
end

function repartition_relais(instance)
    cpt = 1
    println("\n    Nombre de relais par type de communication : ")
    println("    -------------------------------------------- ")
    for i in 1:instance.contexte.nC
        println("        - Type $cpt : $(instance.contexte.nRc[i])")
        cpt += 1
    end
end

function types_communication_unites(instance)
    println("\n    Types de communication pour chaque unité : ")
    println("    ------------------------------------------ ")
    for u in 1:instance.contexte.nU
        print("        - Unité $u : ")
        for c in 1:length(instance.contexte.Cu[u])
            if c != length(instance.contexte.Cu[u])
                print(instance.contexte.Cu[u][c], ", ")
            else
                println(instance.contexte.Cu[u][c])
            end
        end
    end
end

function types_communication_bases(instance)
    println("\n    Types de communication pour chaque base : ")
    println("    ----------------------------------------- ")
    for b in 1:instance.contexte.nB
        print("        - Base $b : ")
        for c in 1:length(instance.contexte.Cb[b])
            if c != length(instance.contexte.Cb[b])
                print(instance.contexte.Cb[b][c], ", ")
            else
                println(instance.contexte.Cb[b][c])
            end
        end
    end
end

function liste_positions_bases(instance)
    cpt = 1
    println("\n    Positions des bases : ")
    println("    --------------------- ")
    println("        Format : x, y")
    for position in instance.contexte.Eb
        println("        - Position n°$cpt : ($(position.x), $(position.y))")
        cpt += 1
    end
end

function liste_positions_HAPS(instance)
    cpt = 1
    println("\n     Positions de déploiement : ")
    println("     -------------------------- ")
    println("        Format : x, y")
    for position in instance.contexte.E
        println("        - Position n°$cpt : ($(position.x), $(position.y))")
        cpt += 1
    end
end

function liste_portees_relais(instance)
    cpt = 1
    println("\n    Portée des différents types de communication : ")
    println("    --------------------------------------------- ")
    for portee in instance.contexte.D
        println("        - Type n°$cpt : $portee")
        cpt += 1
    end
end

function liste_poids_puissances_HAPS(instance)
    cpt = 1
    println("\n    Poids et puissance maximale de chacun des HAPS : ")
    println("    ------------------------------------------------ ")
    println("        Format : Wh, Ph")
    temp = [[instance.contexte.W[i], instance.contexte.P[i]] for i in 1:instance.contexte.nH]
    for poids_puissance in temp
        println("        - HAPS n°$cpt : $(poids_puissance[1]), $(poids_puissance[2])")
        cpt += 1
    end
end

function liste_poids_puissances_portee_relais(instance)
    cpt = 1
    for c in 1:instance.contexte.nC
        println("\n    Poids, puissance et portée (seuil) des relais de type $c : ")
        println("    --------------------------------------------------------- ")
        println("        Format : wr, pr, sr")
        temp = [[instance.contexte.w[r], instance.contexte.p[r], instance.contexte.s[r]] for r in instance.contexte.R[c]]
        for poids_puissance in temp
            println("        - Relais n°$cpt : $(poids_puissance[1]), $(poids_puissance[2]), $(poids_puissance[3])")
            cpt += 1
        end
    end
end

########## Divers ##########

function clear()
    Base.run(`clear`)
end


