########## Utilitaires ##########

function choix_binaire(message)
    choix = "a"
    while choix != "o" && choix != "n"
        print(message)
        choix = readline()
    end
    return choix
end

function choix_multiple(message, max, min=1)
    choix = "a"
    while tryparse(Int, choix) == nothing || tryparse(Int, choix) < min || tryparse(Int, choix) > max
        print(message)
        choix = readline()
    end
    return parse(Int, choix)
end

function clear()
    Base.run(`clear`)
end

function choix_saisie(message)
    choix = choix_binaire(message)
    if choix == "o"
        a = -1
        b = -1
    elseif choix == "n"
        println(" \n - Tirage uniforme entre [a,b]")
        print("    --> a : ")
        a = parse(Int, readline())
        b = a - 1
        while b < a
            print("    --> b (≥ a) : ")
            b = parse(Int, readline())
        end
    end
    return choix, a, b
end

