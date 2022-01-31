function relaxation_value(relaxed_modele, instance, separation_exacte, callback_vide, only_on_root, connexity_module, expefolder)
    println("solving the relaxation")
    set_optimizer(relaxed_modele, CPLEX.Optimizer)
    relax_integrality(relaxed_modele)

    if separation_exacte==false && separation_heuristique==false
        optimize!(relaxed_modele)
        z2 = objective_value(relaxed_modele)
        lecture_resultats_modele(instance, E1, E2, relaxed_modele, Vertices, Arcs, Ebis, true, connexity_module, expefolder)
        println("relxation solution: $z2")
    end
    if separation_exacte==true
        z2=exact_cutting_plane(relaxed_modele, instance, connexity_module, expefolder)
    end
    if separation_heuristique==true
        z2=heuristic_cutting_plane(relaxed_modele, instance, connexity_module, expefolder)
    end
    return z2
end


function find_heuristic_cuts(x_vals, y_vals, pool, instance)
    res=[]
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
                    # nb_const+=1
                    println("Power : $h, $v, $Rcover, $cut_val, $Rhs")
                    push!(pool,(h,v,Rcover))
                    push!(res,(h,v,Rcover))
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
                # println("nb_const: $nb_const")
                push!(pool,(h,v,Rcover))
                push!(res,(h,v,Rcover))
            end
        end
    end
    return res
end

function find_exact_cuts(x_vals, y_vals, pool, instance)
    res=[]
    for h in 1:nH
        ##Power capacity separation
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
        #Optimize
        set_optimizer(model_sep, CPLEX.Optimizer)
        set_optimizer_attribute(model_sep, "CPX_PARAM_SCRIND", 0)
        optimize!(model_sep)
        z_sep = objective_value(model_sep)
        # println("objective sep: $z_sep")
        println("z_sep: $z_sep")
        if z_sep<0-epsilone
            println("violated power constraint : $z_sep")
            beta_vals = round.(Int, value.(model_sep[:beta]))
            gamma_vals = round.(Int, value.(model_sep[:gamma]))

            v=0
            for v1 in Vertices[2:size(Vertices,1)]
                if beta_vals[v1]>epsilone
                    v=v1
                    break
                end
            end
            #FInd the cover
            RCover=[]
            for r in 1:nR
                if gamma_vals[r]>epsilone
                    push!(RCover,r)
                end
            end

            if (h,v,RCover) ∉ pool
                println("power : $h, $v, $RCover, $z_sep")
                push!(pool,(h,v,RCover))
                push!(res,(h,v,RCover))
            end

        end

        ## Weight capacity separation
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
        #Optimize
        set_optimizer(model_sep, CPLEX.Optimizer)
        set_optimizer_attribute(model_sep, "CPX_PARAM_SCRIND", 0)
        optimize!(model_sep)
        z_sep = objective_value(model_sep)
        println("z_sep: $z_sep")
        if z_sep<0-epsilone
            println("violated weight constraint : $z_sep")
            beta_vals = round.(Int, value.(model_sep[:beta]))
            gamma_vals = round.(Int, value.(model_sep[:gamma]))
            #Find v
            v=0
            for v1 in Vertices[2:size(Vertices,1)]
                if beta_vals[v1]>epsilone
                    v=v1
                    break
                end
            end
            #Find the cover
            RCover=[]
            for r in 1:nR
                if gamma_vals[r]>epsilone
                    push!(RCover,r)
                end
            end
            if (h,v,RCover) ∉ pool
                println("weight : $h, $v, $RCover, $z_sep")
                push!(pool,(h,v,RCover))
                push!(res,(h,v,RCover))
            end
        end
    end
    return res # toutes les contraintes violée, contrainte violée localement
end

function exact_cutting_plane(relaxed_modele, instance, connexity_module, expefolder)
    is_over = false
    pool=[]
    y = relaxed_modele[:y]
    x = relaxed_modele[:x]
    # nb_contraintes = 0
    # for (F, S) in list_of_constraint_types(relaxed_modele)
    #     nb_contraintes += num_constraints(relaxed_modele, F, S)
    # end
    # println("nb_contraintes before cutting plane: $nb_contraintes")
    while !is_over
        optimize!(relaxed_modele)
        x_vals = value.(relaxed_modele[:x])
        y_vals = value.(relaxed_modele[:y])
        res = find_exact_cuts(x_vals, y_vals, pool, instance)
        if length(res)==0
            is_over = true
        else
            for p in res
              @constraint(relaxed_modele, sum(y[r,Ebis[p[2]-1]] for r in p[3]) <=  length(p[3])-x[p[1],Ebis[p[2]-1]])
            end
        end
    end
   z2 = objective_value(relaxed_modele)
   lecture_resultats_modele(instance, E1, E2, relaxed_modele, Vertices, Arcs, Ebis, true, connexity_module, expefolder)
   println("relxation solution: $z2")
   return z2
end


function heuristic_cutting_plane(relaxed_modele, instance, connexity_module, expefolder)
    is_over = false
    pool=[]
    y = relaxed_modele[:y]
    x = relaxed_modele[:x]
    # nb_contraintes = 0
    # for (F, S) in list_of_constraint_types(relaxed_modele)
    #     nb_contraintes += num_constraints(relaxed_modele, F, S)
    # end
    # println("nb_contraintes before cutting plane: $nb_contraintes")
    while !is_over
        optimize!(relaxed_modele)
        x_vals = value.(relaxed_modele[:x])
        y_vals = value.(relaxed_modele[:y])
        res = find_heuristic_cuts(x_vals, y_vals, pool, instance)
        if length(res)==0
            is_over = true
        else
            for p in res
              @constraint(relaxed_modele, sum(y[r,Ebis[p[2]-1]] for r in p[3]) <=  length(p[3])-x[p[1],Ebis[p[2]-1]])
            end
        end
    end
   z2 = objective_value(relaxed_modele)
   lecture_resultats_modele(instance, E1, E2, relaxed_modele, Vertices, Arcs, Ebis, true, connexity_module, expefolder)
   println("relxation solution: $z2")
   return z2
end
