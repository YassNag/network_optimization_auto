function gestion_parametres_GLPK(modele, display_log)
    if choix_verbosite == "o"
        set_optimizer_attribute(modele, "msg_lev", 3)
    end
    return 
end

function gestion_parametres_CBC(modele, display_log)
    choix_verbosite = choix_binaire("\n --> Souhaitez-vous obtenir les détails de la résolution (o/n) ? ")
    if choix_verbosite == "n"
        # CBC verbeux par défaut
        set_optimizer_attribute(modele, "logLevel", 0)
    end
    return
end

function gestion_parametres_SCIP(modele, display_log)
    if display_log == "no"
        # SCIP verbeux par défaut
        set_optimizer_attribute(modele, "display/verblevel", 0)
    end
    return
end

function gestion_parametres_CPLEX(modele, display_log)
    if display_log == "no"
        # CPLEX verbeux par défaut
        set_optimizer_attribute(modele, "CPX_PARAM_SCRIND", 0)
    end
    return
end

function gestion_parametres_Gurobi(modele, display_log)
    if display_log == "no"
        # Gurobi verbeux par défaut
        set_optimizer_attribute(modele, "OutputFlag", 0)
    end
    return
end

function gestion_parametres_Mosek(modele, display_log)
    if display_log == "no"
        # Mosek verbeux par défaut
        set_optimizer_attribute(modele, "LOG", 0)
    end
    return
end
