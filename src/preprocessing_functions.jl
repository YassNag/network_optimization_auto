
function dominant_positions(instance, Vertices, Ebis, E3, E1, U)
    Ve=[]
    for v1 in Vertices[2:size(Vertices,1)]
        for v2 in Vertices[2:size(Vertices,1)]
            if v1 != v2
                pos1 = Ebis[v1-1]
                pos2 = Ebis[v2-1]

                condition1=true
                condition2=true
                condition3=true
                #defining condition 1
                for r in 1:nR
                    if !issubset(deleteat!(E3[pos1][r], findall(x->x==pos2,E3[pos1][r])), E3[pos2][r])
                        condition1=false
                        break
                    end
                end
                #defining condition 2
                for v in  Vertices[2:size(Vertices,1)]
                    if v !=v1 && v !=v2
                        for r in 1:nR
                            if Ebis[v1-1] in E3[Ebis[v-1]][r] && Ebis[v2-1] ∉ E3[Ebis[v-1]][r]
                                condition2=false
                                break
                            end
                        end
                    end
                end
                #defining condition 3
                for u in U
                    for c in 1:length(Cu[u[1]])
                        for r in R[Cu[u[1]][c]]
                            if Ebis[v1-1] in E1[u[1]][c][r][u[2]] && Ebis[v2-1] ∉ E1[u[1]][c][r][u[2]]
                                condition3=false
                            end
                        end
                    end
                end
                #
                if condition1 && condition2 && condition3
                    push!(Ve,(v1,v2)) # (dominé, dominante)
                end
            end
        end
    end
    return Ve
end

function calcul_positions_porte_bases(instance)
    E4 =   Vector{Vector{Vector{Int}}}(undef, nC) #ensemble de positions à portée de la base b via le type c
    for c in 1:nC
        seuil_max_relais = maximum([seuil[r] for r in R[c]])
        E4[c] =   Vector{Vector{Int}}(undef, nB)
        for b in 1:nB
            E4[c][b] = Int[]
                for index_position in 1:nE
                    position = E[index_position]
                    pos_base = Eb[b]
                    if distance(instance, pos_base, position) <= seuil_max_relais
                            push!(E4[c][b], index_position)
                    end
                end
        end
    end
        #println("E4: $E4")
return E4
end

function calcul_unites_couverte_bases(instance)
    Ub=Tuple{Int64,Int64}[]
    sr=Int[]
    for b in 1:nB
        pos_base = Eb[b]
        for t in 1:nT
             for u in 1:nU
                 for c in 1:length(Cu[u])
                    seuil_max_relais = maximum(seuil[r] for r in R[Cu[c][1]])
                    pos_unite = instance.scenario.deplacements[u][t]
                    if distance(instance, pos_unite,pos_base) <= min(seuil_max_relais, Pu[u][c])
                        push!(Ub,(u,t))
                        # d=distance(instance, pos_unite,pos_base)
                        # println("u: $u, t:  $t, distance to base: $d")
                    end
                end
            end
        end
    end
    # println(Ub)
    return Ub
end

function calcul_positions_relais(instance, Ebis)
    E3 = Vector{Vector{Vector{Int}}}(undef, nE)
    # E3[e][r] : ensemble des positions à portée de la position e pour le relais r
    for position in Ebis
        E3[position] = Vector{Vector{Int}}(undef, nR)
        for r in 1:nR
            E3[position][r] = Int[]
            for position_bis in 1:nE
                p1 = E[position]
                p2 = E[position_bis]
                position_active = Ebool[position]
                position_bis_active = Ebool[position_bis]
                if distance_air_air(instance, p1, p2) <= seuil[r] && position != position_bis && position_active == true && position_bis_active == true
                    push!(E3[position][r], position_bis)
                end
            end
        end
    end
    return E3
end


function graph_conctruction(instance, Ebis)
    Vertices=Int64[]
    Arcs=Tuple{Int64,Int64}[]

    G = complete_graph(size(Ebis,1)+1)

    push!(Vertices, 1) #position de la base
    cmp=2
    for e in Ebis
        push!(Vertices, cmp)
        cmp+=1
    end

    #remove the non essential edges
    seuil_max_relais = maximum(seuil)

    for a in edges(G)
            u, v = src(a), dst(a)
            if u==1
                e1=pos_base
                e2=E[Ebis[v-1]]
            end
            if v==1
                e1=E[Ebis[u-1]]
                e2=pos_base
            end
            if u!=1 && v!=1
                e1=E[Ebis[u-1]]
                e2=E[Ebis[v-1]]
            end

            if distance(instance, e1, e2) > seuil_max_relais
                rem_edge!(G, a)
            end
    end

    DG=DiGraph(G)

    for a in edges(DG)
            u, v = src(a), dst(a)
            push!(Arcs, (u,v))
    end

return DG, Vertices, Arcs
end


function calcul_positions_bases(instance)
    E2 = Vector{Vector{Vector{Int}}}(undef, nB)
    # E2[b][r] : ensemble des positions à portée de la base b pour le relais r
    for b in 1:nB
        E2[b] = Vector{Vector{Int}}(undef, nR)
        for r in 1:nR
            E2[b][r] = Int[]
            for index_position in 1:nE
                position = E[index_position]
                pos_base = Eb[b]
                if distance(instance, pos_base, position) <= seuil[r]
                    if !(index_position in E2[b][r])
                        push!(E2[b][r], index_position)
                    end
                end
            end
        end
    end
    return E2
end

function calcul_positions_unites(instance)
    E1 = Vector{Vector{Vector{Vector{Vector{Int}}}}}(undef, nU)
    # E1[u][r][t] : ensemble des positions à portée de l'unité u de type c pour le relais r
    for u in 1:nU
        E1[u] = Vector{Vector{Vector{Vector{Int}}}}(undef,length(Cu[u]))
            for c in 1:length(Cu[u])
                E1[u][c] = Vector{Vector{Vector{Int}}}(undef, nR)
                for r in 1:nR
                    E1[u][c][r] = Vector{Vector{Int}}(undef, nT)
                    for t in 1:nT
                        E1[u][c][r][t] = Int[]
                        for index_position in 1:nE
                            position = E[index_position]
                            pos_unite = instance.scenario.deplacements[u][t]
                            if distance(instance, pos_unite, position) <= min(seuil[r],Pu[u][c])
                                if !(index_position in E1[u][c][r][t])
                                    push!(E1[u][c][r][t], index_position)
                                end
                            end
                        end
                    end
                end
            end
    end
    return E1
end

function distance(instance, p1, p2)
    return sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2 + A^2)
end

function positions_non_couvrantes(instance)
    seuil_max_relais = maximum(seuil)
    total = 0
    for position in 1:nE
        position_active = false
        t = 1
        while t <= nT && position_active == false
             u = 1
             while u <= nU && position_active == false
                pos_unite = instance.scenario.deplacements[u][t]
                if distance(instance, pos_unite, E[position]) <= seuil_max_relais
                    position_active = true #Position à portée d'au moins une unité
                end
                u += 1
             end
            t += 1
        end
        b = 1
        while b <= nB && position_active == false
            pos_base = Eb[b]
            if distance(instance, pos_base, E[position]) <= seuil_max_relais
                position_active = true # Position à portée d'au moins une base
            end
            b += 1
        end
        if position_active == false
            total += 1
            Ebool[position] = false
        end
    end
    total_pourcentage = trunc((total/nE)*100, digits=2)
    println("\n - Nombre de positions désactivées : $total/$(nE) ($(total_pourcentage)%)")
end

function distance_air_air(instance, p1, p2)
    return sqrt((p2.x - p1.x)^2 + (p2.y - p1.y)^2)
end
