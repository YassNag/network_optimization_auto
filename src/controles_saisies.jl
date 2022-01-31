
########## Vérification des saisies utilisateur ##########

function choix_multiple(message, max, min=1)
    choix = "a"
    while tryparse(Int, choix) == nothing || tryparse(Int, choix) < min || tryparse(Int, choix) > max
        print(message)
        choix = readline()
    end
    return parse(Int, choix)
end

function choix_binaire(message)
    choix = "a"
    while choix != "o" && choix != "n"
        print(message)
        choix = readline()
    end
    return choix
end
