function lecture_resultats_modele(instance, E1, E2, modele, Vertices, Arcs, Ebis, relaxed, connexity_module, expefolder)
    global nC = instance.contexte.nC
    global R = instance.contexte.R

    Ub=calcul_unites_couverte_bases(instance)
    U=[(u,t) for u in 1:nU, t in 1:nT if (u,t) ∉  Ub]
    io_sol = open("Expes/"*expefolder*"/Results/Solution.txt", "w");
    if relaxed
        io_relax = open("Expes/"*expefolder*"/Results/RelaxationSolution.txt", "w");
        x = value.(modele[:x])
        y = value.(modele[:y])
        z = value.(modele[:z])
        if connexity_module
            l = value.(modele[:l])
            f = value.(modele[:f])
        end
        println(io_relax, "x[pos,h]")
        for h in 1:nH
            for pos in 1:nE
                xx=x[h,pos]
                if xx !=0
                    println(io_relax, "x[$pos,$h] = $xx")
                end
            end
        end
        println(io_relax, "y[pos,r]")
        for r in 1:nR
            for pos in 1:nE
                yy=y[r,pos]
                if yy !=0
                    println(io_relax, "y[$pos,$r] = $yy")
                end
            end
        end
        println(io_relax, "z[u]")
        for u in U
            if z[u]!=0
                zz=z[u]
                println(io_relax, "z[$u] = $zz")
            end
        end
        if connexity_module
            println(io_relax, "l[a,c]")
            for a in Arcs
                for c in 1:nC
                    s=l[a,c]
                        if s !=0
                            println(io_relax, "l[$a,$c] = $s")
                        end
                end
            end
            println(io_relax, "f[a,v] ")
            for a in Arcs
                for v in  Vertices[2:size(Vertices,1)]
                    flow=f[a,v]
                    if flow != 0
                        println(io_relax, "f[$a,$v] = $flow")
                    end
                end
            end
        end
        close(io_relax)
        return true

    else
        x = round.(Int, value.(modele[:x]))
        y = round.(Int, value.(modele[:y]))
        if connexity_module
            f = round.(Int, value.(modele[:f]))
            l = round.(Int, value.(modele[:l]))
        end

        println("\n========== Affichage des résultats ==========")
        z1 = trunc(Int, objective_value(modele))
        println(io_sol, "\n --> Nombre total d'unités couvertes par la base: ", size(Ub,1))
        println(io_sol, "\n --> Nombre total d'unités couvertes par les HAPS/nombre d'unité non couvertes par les bases: ", z1, "/", size(U,1), " (",trunc((z1/(size(U,1)))*100, digits=2),"%)")
        println(io_sol, "\n --> Nombre total d'unités couvertes/ nombre total d'unité: ", size(Ub,1)+z1, "/", nU*nT, " (",trunc(((size(Ub,1)+z1)/(nU*nT))*100, digits=2),"%)")

        deploiement_haps = [-1 for i in 1:nH] # deploiement_haps[h] : position sur laquelle est déployé le HAPS h (si déployé)
        placement_relais = [-1 for i in 1:nR] # placement_relais[r] : HAPS sur lequel est placé le relais r (si placé)
        for h in 1:nH
            deploye = false
            for pos in 1:nE
                if x[h,pos] == 1
                    position = E[pos]
                    println(io_sol, "\n - HAPS $h déployé à la position ($(position.x), $(position.y))")
                    deploiement_haps[h] = pos
                    for c in 1:nC
                        for r in R[c]
                            if y[r,pos] == 1
                                println(io_sol, "\n     - Relais $r (type $c) placé dans ce HAPS")
                                placement_relais[r] = h
                            end
                        end
                    end
                end
            end
        end

        z = round.(Int, value.(modele[:z])) # On arrondit car Cbc renvoie parfois 0.999999... (tolérance)
        unites_couvertes = Vector{Vector{Vector{Vector{Int}}}}(undef, nU) # unites_couvertes[u][t][c] : ensemble des HAPS couvrant l'unité u au temps t pour le type de communication c (dont dispose l'unité u)
        nb_unites_non_couvertes = [0 for i in 1:nT]
        for u in 1:nU
            unites_couvertes[u] = Vector{Vector{Vector{Int}}}(undef, nT)
            for t in 1:nT
                unites_couvertes[u][t] = [Int[] for c in 1:nC]
                if (u,t) in U
                    if z[(u,t)] == 0 # L'unité u est couverte au temps t
                        for c in 1:length(Cu[u])
                            # On cherche par quel(s) types de communication
                            for h in 1:nH
                                # Par quel HAPS
                                for r in R[Cu[u][c]]
                                    # Quel relais
                                    for pos in E1[u][c][r][t]
                                        # Positions potentiellement couvrantes si relais r placé
                                        if x[h,pos] == 1 && y[r,pos] == 1
                                            # Si HAPS déployé ET relais déployé alors unité couverte par ce HAPS et pour ce type de communication
                                            push!(unites_couvertes[u][t][Cu[u][c]], h)
                                            break
                                        end
                                    end
                                end
                            end
                        end
                elseif z[(u,t)] == 1
                    nb_unites_non_couvertes[t] += 1
                end
            end
            end
        end

        bases_couvertes = Vector{Vector{Vector{Int}}}(undef, nB) # unites_couvertes[b][c] : ensemble des HAPS couvrant l'unité b pour le type de communication c (dont dispose l'unité u)
        for b in 1:nB
            bases_couvertes[b] = [Int[] for c in 1:nC]
            for c in Cb[b]
                for h in 1:nH
                # Par quel HAPS
                    for r in R[c]
                        # Quel relais
                        for pos in E2[b][r]
                            # Positions potentiellement couvrantes si relais r placé
                            if x[h,pos] == 1 && y[r,pos] == 1
                                # Si HAPS déployé ET relais déployé alors base couverte par ce HAPS et pour ce type de communication
                                push!(bases_couvertes[b][c], h)
                                break
                            end
                        end
                   end
                end
            end
        end

        println(io_sol, "\n========== Affichage des résultats ==========")
        transmission_haps = Vector{Vector{Vector{Int}}}(undef, nH) # transmission_haps[h][c] : HAPS vers lequel le HAPS h peut transmettre via le type de communication c

        if connexity_module
            #flow variables
                for a in Arcs
                    for v in  Vertices[2:size(Vertices,1)]
                            flow=f[a,v]
                            if flow !=0
                                if a[1]==1
                                    e1=pos_base
                                    e2=E[Ebis[a[2]-1]]

                                    h1=10000
                                    for h in 1:nH
                                            if x[h,Ebis[a[2]-1]] == 1
                                                h2=h
                                            end
                                    end

                                end
                                if a[2]==1
                                    e2=pos_base
                                    e1=E[Ebis[a[1]-1]]

                                    h2=10000
                                    for h in 1:nH
                                            if x[h,Ebis[a[1]-1]] == 1
                                                h1=h
                                            end
                                    end


                                end
                                if a[1]!=1 && a[2]!=1
                                    e2=E[Ebis[a[2]-1]]
                                    e1=E[Ebis[a[1]-1]]

                                    for h in 1:nH
                                        for h_bis in 1:nH
                                            if x[h,Ebis[a[1]-1]] == 1 && x[h_bis, Ebis[a[2]-1]] == 1
                                                h1=h
                                                h2=h_bis
                                            end
                                        end
                                    end

                                end
                            println(io_sol,  "\n - flow f $v arc $a, ($h1, $h2), pos ($e1, $e2) : $flow")
                            end
                    end
            end


        for h in 1:nH
            transmission_haps[h] = [Int[] for c in 1:nC]
        end

        R=R
        for a in Arcs
            for c in 1:nC
                if a[1]!=1 && a[2]!=1
                    pos=Ebis[a[1]-1]
                    pos_bis=Ebis[a[2]-1]
                    if sum(f[a,v] for v in Vertices[2:size(Vertices,1)])!=0 &&
                        sum(y[r,pos] for r in R[c]) !=0 &&
                        sum(y[r,pos_bis] for r in R[c]) !=0
                            for h in 1:nH
                                for h_bis in 1:nH
                                    if x[h,pos] == 1 && x[h_bis, pos_bis] == 1
                                        push!(transmission_haps[h][c], h_bis)
                                        println(io_sol, (h,h_bis,c))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end


        solution = Solution(1, z1, nb_unites_non_couvertes, deploiement_haps, placement_relais, unites_couvertes, bases_couvertes, transmission_haps)
    end

    close(io_sol)
    return solution
end
