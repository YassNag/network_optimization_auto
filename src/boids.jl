using PyPlot

function initialisation(nC, nUc, trace)
    liste_boids = Vector{Vector{Vector{Float64}}}(undef, nC)
    cpt_b = 1
    for c in 1:nC
        liste_boids[c] = Vector{Vector{Float64}}(undef, nUc[c])
        for u in 1:nUc[c]
            # Positions
            x = rand(0:1000)
            y = rand(0:1000)
            # Vélocités
            vx = 0.0 
            vy = 0.0
            liste_boids[c][u] = [x, y, vx, vy]
            trace[cpt_b][1] = [x, y]
            cpt_b += 1
        end
    end 
    return liste_boids
end

function affichage(nC, nUc, liste_couleurs, liste_plot, liste_boids, trace, t)
    cpt_b = 1
    for c in 1:nC
        X = [b[1] for b in liste_boids[c]]
        Y = [b[2] for b in liste_boids[c]]
        push!(liste_plot, plot(X, Y, liste_couleurs[c], marker="^", linestyle="", markersize=3, zorder=2))
        if t >= 3 && trace_active == true
            for b in 1:nUc[c]
                push!(liste_plot, plot([trace[cpt_b][i][1] for i in 1:3], [trace[cpt_b][i][2] for i in 1:3], liste_couleurs[c], marker="", alpha=0.2, zorder=1))
                cpt_b += 1
            end
        elseif t < 3 && trace_active
            for b in 1:nUc[c]
                push!(liste_plot, plot([trace[cpt_b][i][1] for i in 1:t], [trace[cpt_b][i][2] for i in 1:t], liste_couleurs[c], marker="", alpha=0.2, zorder=1))
                cpt_b += 1
            end
        end
    end
end

function reset_affichage(liste_plot)
    for p in liste_plot
        p[1].remove()
    end
end

function regle_1(nC, nUc, b, liste_boids) 
    # Cohésion : les boids se dirigent vers le centre de masse du groupe (perçu par le boid courant)
    centre_masse = [0.0, 0.0]
    for c in 1:nC
        for b2 in liste_boids[c]
            if b2 != b
                centre_masse[1] += b2[1]
                centre_masse[2] += b2[2]
            end
        end
    end
    centre_masse /= (sum(nUc) - 1)
    #plot([b[1], centre_masse[1]], [b[2], centre_masse[2]], color="orange")
    #centre_masse = [500, 500]
    return [(centre_masse[1] - b[1])*facteur_cohesion/100, (centre_masse[2] - b[2])*facteur_cohesion/100]
end

function dist(b1, b2)
    return sqrt((b2[1] - b1[1])^2 + (b2[2] - b1[2])^2)
end

function regle_2(nC, nUc, b, liste_boids) 
    # Séparation : les boids s'évitent
    repulsion = [0.0, 0.0]
    for c in 1:nC
        for b2 in liste_boids[c]
            if b2 != b && dist(b, b2) < distance_proximite
                repulsion[1] -= (b2[1] - b[1])
                repulsion[2] -= (b2[2] - b[2])
            end
        end
    end
    return repulsion*facteur_repulsion
end

function regle_3(nC, nUc, b, liste_boids) 
    # Alignement : les boids s'alignent sur la même trajectoire
    alignement = [0.0, 0.0]
    for c in 1:nC
        for b2 in liste_boids[c]
            if b2 != b
                alignement[1] += (b2[3])
                alignement[2] += (b2[4])
            end
        end
    end
    alignement /= (sum(nUc) - 1)
    return [alignement[1] + b[3], alignement[2] + b[4]]/facteur_alignement
end 

function regle_4(nC, nUc, b)
    # On force les boids à rester dans le terrain
    repulsion = [0.0, 0.0]
    if b[1] > 900
        repulsion[1] -= 20
    elseif b[1] < 100
        repulsion[1] += 20
    end
    
    if b[2] > 900
        repulsion[2] -= 20
    elseif b[2] < 100
        repulsion[2] += 20
    end
    return repulsion
end

function limiter_vitesse(b)
 vitesse_boid = sqrt(b[3]^2 + b[4]^2) 
 if vitesse_boid > vitesse_limite # trop rapide
     b[3] = (b[3] / vitesse_boid) * vitesse_limite
     b[4] = (b[4] / vitesse_boid) * vitesse_limite     
 end
end

function mouvoir_boids(nC, nUc, liste_boids, trace) 
    cpt_b = 1
    for c in 1:nC
        for b in liste_boids[c]
            # Calcul des vélocités (vecteurs vitesse)
            v1 = regle_1(nC, nUc, b, liste_boids) 
            v2 = regle_2(nC, nUc, b, liste_boids)
            v3 = regle_3(nC, nUc, b, liste_boids)
            v4 = regle_4(nC, nUc, b)
               
            # Mise à jour des vélocités en fonction des règles
            b[3] += v1[1] + v2[1] + v3[1] + v4[1]
            b[4] += v1[2] + v2[2] + v3[1] + v4[2]
            limiter_vitesse(b)
            # Mise à jours des positions en fonction des vélocités
            b[1] += b[3]
            b[2] += b[4]
            trace[cpt_b][3] =  trace[cpt_b][2];  trace[cpt_b][2] =  trace[cpt_b][1];  trace[cpt_b][1] = [b[1], b[2]]
            cpt_b += 1
        end
    end
end

function main()
   global trace_active = true
   global vitesse_limite = 30
   global distance_proximite = 50
   global facteur_repulsion = 1
   global facteur_cohesion = 10
   global facteur_alignement = 10
   nC = 3
   nUc = [5, 5, 5]
   nT = 200
   x_limite = 1000
   y_limite = 1000
   liste_couleurs = ["red", "yellow", "blue", "cyan", "lime", "darkorchid", "fuchsia", "snow"]
   # Initialisation graphique
    ion()
    fig = figure("Boids")
    ax = fig.add_subplot(111)
    grid(false)
    xlim(0, x_limite)
    ylim(0, y_limite)
    # Création du sol
    sol = matplotlib.patches.Rectangle((0, 0), 1000, 1000, linewidth=1, edgecolor="black", facecolor="#40916c")
    ax.add_patch(sol)
    # Grille
    pas = 50
    i = pas
    while i < 1000
        ax.plot([i, i], [0, 1000], color="white", linewidth=0.2) # Longueur
        ax.plot([0, 1000], [i, i], color="white", linewidth=0.2) # Largeur
        i += pas
    end 
   
   # Initialisation boids
   trace = [Vector{Vector{Float64}}(undef, 3) for i in 1:sum(nUc)]
   for i in 1:length(trace)
       trace[i] = [Vector{Float64}(undef, 2) for j in 1:3]
   end
   liste_boids = initialisation(nC, nUc, trace)

   # Boucle
   liste_plot = []
   for t in 1:nT
       affichage(nC, nUc, liste_couleurs, liste_plot, liste_boids, trace, t)
       savefig("images/image_"*string(t))
       mouvoir_boids(nC, nUc, liste_boids, trace) 
       sleep(0.1) 
       reset_affichage(liste_plot)
       liste_plot = []
   end
   close()
end

main()
