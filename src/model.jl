include("preprocessing_functions.jl")
include("resultat.jl")
include("relaxation_value.jl")

function resolution_modele(instance, connexity_module, ajout_R3, ajout_R4, callback_vide, ajout_R5, only_on_root, frequency, separation_exacte, solver, display_log, temps_limite, ResultDict, expefolder)
    global epsilone=0.1
    global nT = instance.contexte.nT
    global nH = instance.contexte.nH
    global nE = instance.contexte.nE
    global nC = instance.contexte.nC
    global nU = instance.contexte.nU
    global nB = instance.contexte.nB
    global nR = instance.contexte.nR
    global nRc = instance.contexte.nRc
    global R = instance.contexte.R
    global E=instance.contexte.E
    global Cb = instance.contexte.Cb
    global Cu = instance.contexte.Cu
    global W = instance.contexte.W
    global P = instance.contexte.P
    global w = instance.contexte.w
    global p =instance.contexte.p
    global Ebool=instance.contexte.Ebool
    global Eb=instance.contexte.Eb
    global seuil = instance.contexte.s
    global A = instance.contexte.A
    global Pu=instance.contexte.Pu
    global pos_base = Eb[1]

    # Pré-traitement
    global nb_const = 0
    temps_pretraitement_positions = @elapsed positions_non_couvrantes(instance)
    temps_precalcul_unites = @elapsed global E1 = calcul_positions_unites(instance)
    temps_precalcul_E2 = @elapsed global E2 = calcul_positions_bases(instance)
    temps_precalcul_E4 = @elapsed global E4 = calcul_positions_porte_bases(instance)

    # Positions couvrantes
    global Ebis = Int[] # On cherche les positions actives
    for position in 1:nE
        if Ebool[position] == true # position active
            push!(Ebis, position)

        end
    end

    global DG, Vertices, Arcs=graph_conctruction(instance, Ebis)

    temps_precalcul_E3 = @elapsed E3 = calcul_positions_relais(instance, Ebis) # positions à portée d'une autre position pour un relais donné

    # Création du modèle
    println("\n========== Création du modèle ==========")
    global separation_heuristique=false
    ajout_R1 = true
    ajout_R2 = true
    if callback_vide == false
        if ajout_R5 == true
            if separation_exacte==false
                println("la séparation sera effectuée avec la méthode heuristique.")
                separation_heuristique = true
            end
        else
            separation_exacte=false
            separation_heuristique=false
        end
    else
        ajout_R5=false
        separation_exacte=false
        separation_heuristique = false
    end
    relaxed=false
    temp_creation_modele = @elapsed modele = creation_modele(instance, E1, E2, E4, ajout_R1, ajout_R2, connexity_module, ajout_R3, ajout_R4, ajout_R5, relaxed,
    separation_exacte, callback_vide, only_on_root, frequency, solver, display_log, temps_limite, ResultDict, expefolder)


    println("\n - Nombre de variables : ", num_variables(modele))
    nb_contraintes = 0
    for (F, S) in list_of_constraint_types(modele)
        nb_contraintes += num_constraints(modele, F, S)
    end
    println("\n - Nombre de contraintes : ", nb_contraintes)

    ResultDict[:Var]=num_variables(modele)
    ResultDict[:Const]=nb_contraintes

    println("\n========== Création du modèle ==========")

    # Choix du solveur
    if solver == "GLPK"
        set_optimizer(modele, GLPK.Optimizer)
        gestion_parametres_GLPK(modele, display_log)
    elseif solver == "Cbc"
        set_optimizer(modele, Cbc.Optimizer)
        gestion_parametres_CBC(modele, display_log)
    elseif solver == "SCIP"
         set_optimizer(modele, SCIP.Optimizer)
         gestion_parametres_SCIP(modele, display_log)
    elseif solver == "Gurobi"
        print("\n - ")
        set_optimizer(modele, Gurobi.Optimizer)
        gestion_parametres_Gurobi(modele, display_log)
    elseif solver == "CPLEX"
        set_optimizer(modele, CPLEX.Optimizer)
        gestion_parametres_CPLEX(modele, display_log)
    elseif solver == "Mosek"
        set_optimizer(modele, Mosek.Optimizer)
        gestion_parametres_Mosek(modele, display_log)
    end

    # Limite de temps
    # temps_limite = parse(Float64, temps_limite)
    set_time_limit_sec(modele, temps_limite)

    # Résolution
    #Solving the relaxation
    relaxed=true
    relaxed_modele = creation_modele(instance, E1, E2, E4, ajout_R1, ajout_R2, connexity_module, ajout_R3, ajout_R4,ajout_R5, relaxed, separation_exacte,
    callback_vide, only_on_root, frequency, solver, display_log, temps_limite, ResultDict, expefolder)
    z2= relaxation_value(relaxed_modele, instance, separation_exacte, callback_vide, only_on_root, connexity_module, expefolder)
    relaxed=false

    if display_log == "yes"
        println("\n========== Détails de la résolution ==========\n")
        optimize!(modele)
        print("\n========== Détails de la résolution ==========")
    else
        print("\n - Résolution en cours...")
        optimize!(modele)
    end


    # Affichage des temps
    precision = 4 # Précision des résultats (nombre de décimales)
    temps_resolution = solve_time(modele)
    ResultDict[:CPU]=round(temps_resolution, digits=2)
    println("\n\n========== Affichage des temps ==========")
    println("\n - Temps suppression de positions non couvrantes : $(trunc(temps_pretraitement_positions, digits=precision)) seconde(s)")
    println("\n - Temps pré-calcul E1 (ensemble des positions à portée de l'unité u de type c pour le relais r): $(trunc(temps_precalcul_unites, digits=precision)) seconde(s)")
    println("\n - Temps pré-calcul E2 (ensemble des positions pouvant atteindre la base b pour le relais r) : $(trunc(temps_precalcul_E2, digits=precision)) seconde(s)")
    println("\n - Temps pré-calcul E3 (positions à portée d'une autre position pour un relais donné) : $(trunc(temps_precalcul_E3, digits=precision)) seconde(s)")
    println("\n - Temps pré-calcul E4 (ensemble de positions à portée de la base b via le type c) : $(trunc(temps_precalcul_E4, digits=precision)) seconde(s)")
    println("\n - Temps création modèle : $(trunc(temp_creation_modele, digits=precision)) seconde(s)")
    println("\n - Temps résolution : $(trunc(temps_resolution, digits=precision)) seconde(s)")
    println("\n========== Affichage des temps ==========")

    if termination_status(modele) == MOI.INFEASIBLE
        println("\n - Échec : problème non réalisable")
        statut = 0
    elseif termination_status(modele) == MOI.OPTIMAL # Solution optimale
        println("\n - Solution optimale trouvée")
        statut = 1
    elseif termination_status(modele) == MOI.TIME_LIMIT && has_values(modele) # Temps limite atteint mais solution réalisable sous-optimale disponible
        println("\n - Solution réalisable (possiblement sous-optimale) trouvée dans le temps imparti")
        statut = 1
    elseif termination_status(modele) == MOI.TIME_LIMIT # Temps limite atteint et pas de solution
        println("\n - Échec : temps limite de résolution atteint (> $(temps_limite/(10^3)) seconde(s)) ")
        statut = 0
    elseif termination_status(modele) == MOI.NODE_LIMIT
            println("\n -  : arret noeud 0 objective")
            statut = 1
    end

    if statut == 1 # Solution disponible
        Ebis = Int[] # On cherche les positions actives (sous-ensemble de positions)
        for position in 1:nE
            if Ebool[position] == true # position active
                push!(Ebis, position)
            end
        end

        Ub=calcul_unites_couverte_bases(instance)
        U=[(u,t) for u in 1:nU, t in 1:nT if (u,t) ∉  Ub]

        solution = lecture_resultats_modele(instance, E1, E2, modele, Vertices, Arcs, Ebis, relaxed, connexity_module, expefolder)
    elseif statut == 0 # Pas de solution
        solution = Solution(0, 0, Int[], Int[], Int[], Int[], Int[])
    end

    z1=solution.z
    if z1==0
        Gap=0
    else
        println("solopt: $z1")
        println("solFractionnaire: $z2")
        Gap= ((z2-z1)/z1)*100
    end

     println("Gap: $Gap")
     ResultDict[:Gap]=round(Gap, digits=2)
     if separation_exacte==true || separation_heuristique ==true
         println("nbre d'inégalité de couverture ajoutés: $nb_const")
         ResultDict[:CoverI]=nb_const
     end

     number_of_nodes=MOI.get(modele, MOI.NodeCount())
     println("Number of nodes in the B&C tree : $number_of_nodes")
     ResultDict[:Nodes]=number_of_nodes

    #solve the relaxation
    return solution
end



function creation_modele(instance, E1, E2, E4, ajout_R1, ajout_R2, connexity_module, ajout_R3, ajout_R4, ajout_R5, relaxed, separation_exacte,
    callback_vide, only_on_root, frequency, solver, display_log, temps_limite, ResultDict, expefolder)
    # Initialisation du modèle
    modele = Model()
    Ebis = Int[] # On cherche les positions actives (sous-ensemble de positions)
    Ebisbar = Int[] # On cherche les positions actives (sous-ensemble de positions)

    for position in 1:nE
        if Ebool[position] == true # position active
            push!(Ebis, position)
        else
            push!(Ebisbar, position)

        end
    end


    E3 = calcul_positions_relais(instance, Ebis) # positions à portée d'une autre position pour un relais donné

    #Ensemble des unités à couvrir.
    Ub=calcul_unites_couverte_bases(instance)
    U=[(u,t) for u in 1:nU, t in 1:nT if (u,t) ∉  Ub]

    # Variables
    if relaxed && connexity_module
        @variables(modele, begin
            x[1:nH, 1:nE]>= 0
            y[1:nR, 1:nE]>= 0
            z[U]>= 0
            l[Arcs, 1:nC]>= 0
            f[Arcs, Vertices[2:size(Vertices,1)]] >= 0
            end)
    end

    if relaxed && !connexity_module
        @variables(modele, begin
            x[1:nH, 1:nE]>= 0
            y[1:nR, 1:nE]>= 0
            z[U]>= 0
            end)
    end

    if !relaxed && connexity_module
        @variables(modele, begin
            x[1:nH, 1:nE], Bin
            y[1:nR, 1:nE], Bin
            z[U] >= 0
            f[Arcs, Vertices[2:size(Vertices,1)]] >= 0
            l[Arcs, 1:nC], Bin
            end)
    end

    if !relaxed && !connexity_module
        @variables(modele, begin
            x[1:nH, 1:nE], Bin
            y[1:nR, 1:nE], Bin
            z[U] >= 0
            end)
    end

    if relaxed
        @constraint(modele, c_r1[h=1:nH, pos=1:nE], x[h,pos]<=1)
        @constraint(modele, c_r2[r=1:nR, pos=1:nE], y[r,pos]<=1)
        @constraint(modele, c_r3[u in U], z[u]<=1)
        if connexity_module
            @constraint(modele, c_r4[a in Arcs,v in Vertices[2:size(Vertices,1)]], f[a,v]<=1)
            @constraint(modele, c_r5[a in Arcs, c=1:nC], l[a,c]<=1)
        end
    end

    # Contraintes
    @constraint(modele, c_01[h=1:nH, pos in Ebisbar], x[h,pos]==0)

    @constraint(modele, c_02[r=1:nR, pos in Ebisbar], y[r,pos] == 0)

    @constraint(modele, c_03[h=1:nH, pos in Ebis], x[h,pos]  <= sum(y[r,pos] for r in 1:nR))

    @constraint(modele, c_1[h=1:nH], sum(x[h,pos] for pos in Ebis) <= 1)

    @constraint(modele, c_2[pos in Ebis], sum(x[h,pos] for h in 1:nH) <= 1)

    @constraint(modele, c_3[r=1:nR], sum(y[r,pos] for pos in Ebis) <= 1)

    @constraint(modele, c_4[pos in Ebis], sum(w[r] * y[r,pos] for r in 1:nR) <= sum(W[h] * x[h,pos] for h in 1:nH))

    @constraint(modele, c_5[pos in Ebis], sum(p[r] * y[r,pos] for r in 1:nR) <= sum(P[h] * x[h,pos] for h in 1:nH))

    @constraint(modele, c_6[u in U], z[u] <= sum(y[r,pos] for c in 1:length(Cu[u[1]]), r in R[Cu[u[1]][c]], pos in E1[u[1]][c][r][u[2]]) )

    @constraint(modele, c_r33[u in U], z[u]<=1)

    @objective(modele, Max, sum(z[u] for u in U))

    if connexity_module
        # #Flow conservation f
        for v in Vertices[2:size(Vertices,1)]
                for i in 1:nv(DG)
                    if i!=v && i!=1
                        @constraint(modele, sum(f[(i,j),v] for j in outneighbors(DG, i)) -
                        sum(f[(j,i),v]  for j in inneighbors(DG, i))  ==0)
                    end
                    if i==v
                        @constraint(modele, sum(f[(i,j),v]  for j in outneighbors(DG, i)) -
                        sum( f[(j,i),v] for j in inneighbors(DG, i))  == - sum(x[h,Ebis[v-1]] for h in 1:nH))
                    end
                    if i==1
                        @constraint(modele, sum(f[(i,j),v]  for j in outneighbors(DG, i)) -
                        sum(f[(j,i),v] for j in inneighbors(DG, i))  == sum(x[h,Ebis[v-1]] for h in 1:nH))
                    end
                end
            end


        # # Defining l variables
        for a in Arcs
                for c in 1:nC
                    #finding the relays of type c that can reach j from i
                    R_cj=Int64[]
                    if a[1]!=1 && a[2]!=1
                        position = Ebis[a[1]-1]
                        for r in R[c]
                            if Ebis[a[2]-1] in E3[position][r]
                                push!(R_cj, r)
                            end
                        end
                        @constraint(modele, l[a,c] <=
                        sum(y[r,Ebis[a[1]-1]] for r in R_cj))

                        @constraint(modele, l[a,c] <=
                        sum(y[r,Ebis[a[2]-1]] for r in R[c]))

                        # @constraint(modele, sum(y[r,Ebis[a[2]-1]] for r in R[c])+sum(y[r,Ebis[a[1]-1]] for r in R_cj)-1
                        # <= l[a,c])

                    end

                    if a[2]==1
                        if c in Cb[1]
                            @constraint(modele, l[a,c] == sum(y[r,Ebis[a[1]-1]] for r in R[c]
                            if Ebis[a[1]-1] in E2[1][r]))
                        else
                            @constraint(modele, l[a,c] == 0)
                        end
                    end

                    if a[1]==1
                        if c in Cb[1]
                            if  Ebis[a[2]-1] in E4[c][1]
                                @constraint(modele, l[a,c] == sum(y[r,Ebis[a[2]-1]] for r in R[c]))
                            else
                                @constraint(modele, l[a,c] == 0)
                            end
                        else
                            @constraint(modele, l[a,c] == 0)
                        end
                    end

                end
        end

        #Symétrie forte
        @constraint(modele,c_sym[a in Arcs, c=1:nC], l[(a[2],a[1]),c] == l[(a[1],a[2]), c])

        #Flow interdiction
        for a in Arcs
            for v in Vertices[2:size(Vertices,1)]
                @constraint(modele, f[a,v] <= sum(l[a,c] for c in 1:nC))
            end
        end

        for a in Arcs
            for v in Vertices[2:size(Vertices,1)]
                @constraint(modele, f[a,v] <= sum(x[h,Ebis[v-1]] for h in 1:nH))
            end
        end

        for v in Vertices[2:size(Vertices,1)]
            @constraint(modele, sum(f[(i,1),v]  for i in inneighbors(DG, 1)) == 0)
            @constraint(modele, sum(f[(v,i),v]  for i in outneighbors(DG, v)) == 0)
        end

            if ajout_R1 == true
                @constraint(modele, c_R1[pos in Ebis, c=1:nC], sum(y[r,pos] for r in R[c]) <= sum(x[h, pos] for h in 1:nH))
            end

            if ajout_R2 == true
                @constraint(modele, c_R22[pos in Ebis, r=1:nR], y[r, pos] <= sum(x[h, pos] for h in 1:nH))
            end

        if ajout_R3 == true
            Ve=dominant_positions(instance, Vertices, Ebis, E3, E1, U)
            Lve=length(Ve)
            println("nombre de couples de positions dominés: $Lve")
            ResultDict[:DomConst]=Lve

            @constraint(modele, c_domxx[v in Ve,  h in 1:nH], x[h,Ebis[v[1]-1]] <=sum(y[r1,Ebis[v[2]-1]]  for r1 in 1:nR))
        end

        if ajout_R4 == true
            Ur=Dict()
            for v in Vertices[2:size(Vertices,1)]
                Ur[v]=[]
                for c in 1:nC
                    for r in R[c]
                        if r ∉ Ur[v]
                            exist=false
                            for u in U
                                if c in Cu[u[1]]
                                    cc=findfirst(isequal(c), Cu[u[1]])
                                    if Ebis[v-1] in E1[u[1]][cc][r][u[2]]
                                        exist=true
                                    end
                                end
                            end
                            if exist
                                push!(Ur[v],r)
                            end
                        end
                    end
                end
            end

            for v in Vertices[2:size(Vertices,1)]
                    @constraint(modele, 2 * sum(x[h, Ebis[v-1]] for h in 1:nH) - sum(y[r,Ebis[v-1]] for r in Ur[v]) <= sum(l[(v,j),c] for c in 1:nC for j in outneighbors(DG, v)))
            end
        end
        if ajout_R5
            if relaxed==false
                global pool=[]
                ###empty callback
                function my_empty_callback_function(cb_data)

                end

                ###exact callback
                function my_exact_callback_function(cb_data)
                    n = Ref{CPXLONG}()
                    CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODECOUNT, n)
                    x_vals = callback_value.(Ref(cb_data), x)
                    y_vals = callback_value.(Ref(cb_data), y)

                    if  only_on_root
                        callbackcond = n[] ==0
                    else
                        callbackcond = n[] % frequency ==0
                    end
                    if callbackcond
                        #Model
                        for h in 1:nH
                            model_sep = Model()
                            #Variables
                            @variables(model_sep, begin
                                beta[Vertices[2:size(Vertices,1)]], Bin
                                gamma[1:nR], Bin
                                delta[Vertices[2:size(Vertices,1)], 1:nR], Bin
                                end)
                            #Constraintes
                            @constraint(model_sep, sum(beta[v] for v in Vertices[2:size(Vertices,1)]) ==1)
                            @constraint(model_sep, sum(p[r]*gamma[r] for r in 1:nR ) >= sum((P[h]+1)*beta[v] for v in Vertices[2:size(Vertices,1)]))
                            for r in 1:nR
                                for v in Vertices[2:size(Vertices,1)]
                                    @constraint(model_sep, beta[v] + gamma[r] -1 <= delta[v,r])
                                end
                            end
                            #Objective
                            @objective(model_sep, Min, sum((1-y_vals[r,Ebis[v-1]])*delta[v,r] for v in Vertices[2:size(Vertices,1)] for r in 1:nR)-
                                sum(x_vals[h,Ebis[v-1]]*beta[v] for v in Vertices[2:size(Vertices,1)]))
                            # println("model created")
                            #Optimize
                            set_optimizer(model_sep, CPLEX.Optimizer)
                            set_optimizer_attribute(model_sep, "CPX_PARAM_SCRIND", 0)
                            optimize!(model_sep)
                            z_sep = objective_value(model_sep)
                            # println("objective sep: $z_sep")
                            if z_sep<0-epsilone
                                println("violated power constraint : $z_sep")
                                beta_vals = round.(Int, value.(model_sep[:beta]))
                                gamma_vals = round.(Int, value.(model_sep[:gamma]))
                                # println("beta_vals: $beta_vals")
                                # println("gamma_vals: $gamma_vals")
                                v=0
                                for v1 in Vertices[2:size(Vertices,1)]
                                    if beta_vals[v1]>epsilone
                                        v=v1
                                        break
                                    end
                                end
                                    RCover=[]
                                    for r in 1:nR
                                        if gamma_vals[r]>epsilone
                                            push!(RCover,r)
                                        end
                                    end
                                    if (h,v,RCover) ∉ pool
                                        println("power : $h, $v, $RCover, $z_sep")
                                        global nb_const+=1
                                        con = @build_constraint(
                                        sum(y[r,Ebis[v-1]] for r in RCover) <=  length(RCover)-x[h,Ebis[v-1]]
                                        )
                                        # println("Adding $(con)")
                                        MOI.submit(modele, MOI.UserCut(cb_data), con)
                                        # println("nb_const: $nb_const")
                                        push!(pool,(h,v,RCover))
                                    end
                            end

                            ## Weight separation
                            model_sep = Model()
                            #Variables
                            @variables(model_sep, begin
                                beta[Vertices[2:size(Vertices,1)]], Bin
                                gamma[1:nR], Bin
                                delta[Vertices[2:size(Vertices,1)], 1:nR], Bin
                                end)
                            #Constraintes
                            @constraint(model_sep, sum(beta[v] for v in Vertices[2:size(Vertices,1)]) ==1)
                            @constraint(model_sep, sum(w[r]*gamma[r] for r in 1:nR ) >= sum((W[h]+1)*beta[v] for v in Vertices[2:size(Vertices,1)]))
                            for r in 1:nR
                                for v in Vertices[2:size(Vertices,1)]
                                    @constraint(model_sep, beta[v] + gamma[r] -1 <= delta[v,r])
                                end
                            end
                            #Objective
                            @objective(model_sep, Min, sum((1-y_vals[r,Ebis[v-1]])*delta[v,r] for v in Vertices[2:size(Vertices,1)] for r in 1:nR)-
                                sum(x_vals[h,Ebis[v-1]]*beta[v] for v in Vertices[2:size(Vertices,1)]))
                            # println("model created")
                            #Optimize
                            set_optimizer(model_sep, CPLEX.Optimizer)
                            set_optimizer_attribute(model_sep, "CPX_PARAM_SCRIND", 0)
                            optimize!(model_sep)
                            z_sep = objective_value(model_sep)
                            # println("objective sep: $z_sep")
                            if z_sep<0-epsilone
                                println("violated weight constraint : $z_sep")
                                beta_vals = round.(Int, value.(model_sep[:beta]))
                                gamma_vals = round.(Int, value.(model_sep[:gamma]))
                                # println("beta_vals: $beta_vals")
                                # println("gamma_vals: $gamma_vals")
                                v=0
                                for v1 in Vertices[2:size(Vertices,1)]
                                    if beta_vals[v1]>epsilone
                                        v=v1
                                        break
                                    end
                                end
                                    RCover=[]
                                    for r in 1:nR
                                        if gamma_vals[r]>epsilone
                                            push!(RCover,r)
                                        end
                                    end
                                    if (h,v,RCover) ∉ pool
                                        println("weight : $h, $v, $RCover, $z_sep")
                                        global nb_const+=1
                                        con = @build_constraint(
                                        sum(y[r,Ebis[v-1]] for r in RCover) <=  length(RCover)-x[h,Ebis[v-1]]
                                        )
                                        # println("Adding $(con)")
                                        #MOI.submit(modele, MOI.UserCut(cb_data), con)
                                        MOI.submit(modele, MOI.UserCut(cb_data), con)
                                        # println("nb_const: $nb_const")
                                        push!(pool,(h,v,RCover))
                                    end
                            end
                        end
                    end
                end

                ###heuristic callback
                pool=[]
                function my_heuristic_callback_function(cb_data)
                    n = Ref{CPXLONG}()
                    CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODECOUNT, n)
                    x_vals = callback_value.(Ref(cb_data), x)
                    y_vals = callback_value.(Ref(cb_data), y)

                    if  only_on_root
                        callbackcond = n[] ==0
                    else
                        callbackcond = n[] % frequency ==0
                    end
                    if callbackcond
                        for h in 1:nH
                            vmax=0
                            max=0
                            for v in Vertices[2:size(Vertices,1)]
                                if x_vals[h,Ebis[v-1]] >= max
                                    vmax=v
                                    max= x_vals[h,Ebis[v-1]]
                                end
                            end
                            v=vmax

                            Rw=[i for i in 1:nR]
                            Rp=Rw
                            sort!(Rw, by = r -> y_vals[r,Ebis[v-1]]*p[r], rev=true)
                            sort!(Rp, by = r -> y_vals[r,Ebis[v-1]]*w[r], rev=true)

                            #Power cover inequalities
                            weight=0
                            Rcover=Int[]
                            cpt=1
                            while weight <= P[h] && cpt <= nR
                                weight+=p[Rp[cpt]]
                                push!(Rcover, Rp[cpt])
                                # println(weight)
                                cpt+=1
                            end
                            if  weight > P[h] && Rcover!=[]
                                cut_val=sum(y_vals[r,Ebis[v-1]] for r in Rcover)
                                Rhs=length(Rcover)-x_vals[h,Ebis[v-1]]
                                # println("cut power value: $cut_val, Rhs: $Rhs")
                                if cut_val > Rhs+epsilone && (h,v,Rcover) ∉ pool
                                        global nb_const+=1
                                        println("Power : $h, $v, $Rcover, $cut_val, $Rhs")
                                        con = @build_constraint(
                                        sum(y[r,Ebis[v-1]] for r in Rcover) <=  length(Rcover)-x[h,Ebis[v-1]]
                                        )
                                        # println("Adding $(con)")
                                        #MOI.submit(modele, MOI.UserCut(cb_data), con)
                                        MOI.submit(modele, MOI.UserCut(cb_data), con)
                                        println("nb_const: $nb_const")
                                        push!(pool,(h,v,Rcover))
                                        # println(pool)
                                end
                            end

                            #weight cover inequalities
                            weight=0
                            Rcover=Int[]
                            Rr=[[r for r in R[c]] for c in 1:nC]
                            for c in 1:nC
                                sort!(Rr[c], by = r -> y_vals[r,Ebis[v-1]]/w[r], rev=true)
                            end
                            cpt=1
                            while weight <= W[h] && cpt <= nR
                                weight+=w[Rw[cpt]]
                                push!(Rcover, Rw[cpt])
                                # println(weight)
                                cpt+=1
                            end
                            # println("Cover : $Rcover")
                            # println(weight, W[h])
                            if  weight > W[h] && Rcover!=[]
                                cut_val=sum(y_vals[r,Ebis[v-1]] for r in Rcover)
                                Rhs=length(Rcover)-x_vals[h,Ebis[v-1]]
                                # println("cut power value: $cut_val, Rhs: $Rhs")
                                if cut_val > Rhs + epsilone && (h,v,Rcover) ∉ pool
                                    println("Weight : $h, $v, $Rcover, $cut_val, $Rhs")
                                    global nb_const+=1
                                    con = @build_constraint(
                                    sum(y[r,Ebis[v-1]] for r in Rcover) <=  length(Rcover)-x[h,Ebis[v-1]]
                                    )
                                    # println("Adding $(con)")
                                    #MOI.submit(modele, MOI.UserCut(cb_data), con)
                                    MOI.submit(modele, MOI.UserCut(cb_data), con)
                                    # println("nb_const: $nb_const")
                                    push!(pool,(h,v,Rcover))
                                end
                            end
                        end
                    end
                end

                if separation_exacte==true
                    MOI.set(modele, MOI.UserCutCallback(), my_exact_callback_function)
                end
                if separation_heuristique==true
                    MOI.set(modele, MOI.UserCutCallback(), my_heuristic_callback_function)
                end
                if callback_vide==true
                    MOI.set(modele, MOI.UserCutCallback(), my_empty_callback_function)
                end
            end
        end
    end
    return modele
end
