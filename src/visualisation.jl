function customisation_pyplot()
    # On désactive certaines touches de base
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["keymap.pan"] = ""
    rcParams["toolbar"] = "None"
    rcParams["keymap.all_axes"] = ""
end

function fullscreen()
    manager = get_current_fig_manager()
    manager.window.showMaximized()
end

function on_press(event, demarrer, continuer, pause, nb_images_seconde, trace_unite)
    if demarrer[1] == false
        if event.key == "d" # Démarrer la visualisation
            demarrer[1] = true
        end
    else # Les autres commandes sont bloquées tant que la visualisation n'a pas été lancée
        if event.key == "i" # Interrompre la visualisation
            pause[1] = false
            continuer[1] = false
        elseif event.key == "p" # Mettre en pause ou reprendre la visualisation
            pause[1] = !pause[1]
        elseif event.key == "a" # Accélerer la visualisation (on régule le temps de pause entre chaque frame)
            nb_images_seconde[1] += 1
        elseif event.key == "r" # Ralentir la visualisation
            if nb_images_seconde[1] > 1
                nb_images_seconde[1] -= 1
            end
        elseif event.key == "t" # Afficher la trace des unités
            trace_unite[1] = !trace_unite[1]
        end
    end
end

function initialisation_suite(ax, instance, solution, liste_couleurs)
    affichage_nb_unites_non_couvertes = ax.text2D(-0.25, 0.75, "", transform=ax.transAxes, fontsize=12, bbox=Dict("facecolor"=>"whitesmoke", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8),  zorder=1000)
    total = instance.contexte.nT * instance.contexte.nU
    ax.text2D(-0.25, 0.85, "Nombre total d'unités non couvertes : $(solution.z)/$(total) ($(trunc((solution.z/total)*100, digits=2))%)", transform=ax.transAxes, fontsize=12,  bbox=Dict("facecolor"=>"whitesmoke", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=1000)
    ax.plot([], [], [], color="darkorange", marker="o", markeredgecolor="black", linestyle="", markersize=20,  label="HAPS")
    ax.plot([], [], [], color="darkorange", marker="s", linestyle="", markersize=10, label="Déploiement possible")
    for position in 1:instance.contexte.nE
        if instance.contexte.Ebool[position] == false # position inactive
            ax.plot([], [], [], color="red", marker="s", linestyle="", markersize=10, label="Déploiement impossible")
            break
        end
    end
    ax.plot([], [], [], color="lightgray", marker="s", markeredgecolor="black", linestyle="", markersize=15, label="Base")
    ax.legend(loc="upper left", bbox_to_anchor=(0.90, 0.95), fontsize=15, framealpha=0.8, facecolor="whitesmoke", edgecolor="black").set_zorder(1000)
    return affichage_nb_unites_non_couvertes
end

function initialisation_base(instance, liste_couleurs, mode)
    ion()
    using3D()
    fig = figure("Outil de visualisation") # figsize=(16,9), tight_layout=true
    #fullscreen()
    ax = fig.add_subplot(111, projection="3d")
    ax.margins(0)
    ax.view_init(elev=15, azim=-5)
    ax.dist = 10
    fig.set_facecolor("#14213d")
    ax.set_facecolor("#14213d")
    axis("off")
    ax.grid(false)
    ax.set_zticks([0, 15, 20, 50, 85, 500])
    xlim(0, 1000)
    ylim(0, 1000)
    zlim(0, 50)
    tight_layout(pad=0)
    if mode == 1
        suptitle("Visualisation de la solution",  y=0.95, bbox=Dict("facecolor"=>"whitesmoke", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), fontsize=15)
    elseif mode == 2
        suptitle("Visualisation du scénario", y=0.95, bbox=Dict("facecolor"=>"whitesmoke", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), fontsize=15)
    end
    affichage_nb_images_seconde = ax.text2D(-0.25, 0.90, "", transform=ax.transAxes, fontsize=12,  bbox=Dict("facecolor"=>"whitesmoke", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=1000)
    affichage_etat_courant = ax.text2D(-0.25, 0.80, "", transform=ax.transAxes, fontsize=12,  bbox=Dict("facecolor"=>"whitesmoke", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=1000)
    ax.plot([], [], [], color="white", marker="^", linestyle="", markersize=10, markeredgecolor="black", label="Unités terrestres")
    for i in 1:instance.contexte.nC
        if mode == 1
            ax.plot([], [], [], color=liste_couleurs[i], marker=".", markersize=10, label="Relais de type $i")
        end
    end
    generation_sol(ax)
    generation_ciel(ax)
    generation_etoiles(ax)
    ax.legend(loc="upper left", bbox_to_anchor=(0.90, 0.95), fontsize=15, framealpha=0.8, facecolor="whitesmoke", edgecolor="black").set_zorder(1000)
    return fig, ax, affichage_nb_images_seconde, affichage_etat_courant
end

function generation_sol(ax)
    # Meshgrid
    # Ancienne méthode, génère des ralentissements
    #xIt = 0:50:1000
    #yIt = 0:50:1000
    #x = xIt' .* ones(length(xIt))
    #y = ones(length(yIt))' .* yIt
    #z = x .* 0
    #plot_surface(x, y, z, color="#239b5e", edgecolor="none", alpha=1, shade=false) # 2d6a4f
    # Surface
    x = [-10 -10; 1010 1010]
    y = [-10 1010; -10 1010]
    z = [0 0; 0 0]
    plot_surface(x,y,z, color="#40916c", edgecolor="black", alpha=1, shade=false, zorder=-5, clip_on=false)
    # Grille
    pas = 50
    i = 0
    while i < 1010
        ax.plot([i, i], [-10, 1010], [0], color="white", linewidth=0.1, zorder=15, clip_on=false) # Longueur
        ax.plot([-10, 1010], [i, i], [0], color="white", linewidth=0.1, zorder=15, clip_on=false) # Largeur
        i += pas
    end

    # Plan tombant 1
    x = [-10 1010; -10 1010]
    y = [-10 -10; -10 -10]
    z = [0 0; -15 -15]
    plot_surface(x, y, z, color="#774936", edgecolor="black", alpha=1, shade=false, zorder=-5, clip_on=false)

    # Plan tombant 2
    x = [1010 1010; 1010 1010]
    y = [-10 1010; -10 1010]
    z = [0 0; -15 -15]
    plot_surface(x, y, z, color="#774936", edgecolor="black", alpha=1, shade=false, zorder=-5, clip_on=false)

    # Plan tombant 3
    x = [-10 1010; -10 1010]
    y = [1010 1010; 1010 1010]
    z = [0 0; -15 -15]
    plot_surface(x, y, z, color="#774936", edgecolor="black", alpha=1, shade=false, zorder=-5, clip_on=false)

    # Plan tombant 4
    x = [-10 -10; -10 -10]
    y = [-10 1010; -10 1010]
    z = [0 0; -15 -15]
    plot_surface(x, y, z, color="#774936", edgecolor="black", alpha=1, shade=false, zorder=-5, clip_on=false)

    # Fermeture
    x = [-10 -10; 1010 1010]
    y = [-10 1010; -10 1010]
    z = [-15 -15; -15 -15]
    plot_surface(x, y, z, color="#774936", edgecolor="black", alpha=1, shade=false, zorder=-5, clip_on=false)

    # Démarcations (arrêtes)
    #plot([0, 0], [0, 0], [0, -500], color="#081c15", linewidth=1, zorder=10)
    #plot([1000, 1000], [0, 0], [0, -500], color="#081c15", linewidth=1, zorder=10)
    #plot([1000, 1000], [1000, 1000], [0, -500], color="#081c15", linewidth=1, zorder=10)

    # Affichage dimensions
    #ax.plot([0,1000], [-50,-50], [0,0], color="white", marker="", markersize=5, zorder=10)
    #ax.plot([1050,1050], [0,1000], [0,0], color="white", marker="", markersize=5, zorder=10)
end

function generation_ciel(ax)
    ###### Fond n°1 ######

    # Troposphère (15 km)
    x = [-10 1010; -10 1010]
    y = [1010 1010; 1010 1010]
    z = [0 0; 15 15]
    plot_surface(x,y,z, color="#48cae4", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    # Stratosphère (15-50 km)
    x = [-10 1010; -10 1010]
    y = [1010 1010; 1010 1010]
    z = [15 15; 50 50]
    plot_surface(x,y,z, color="#00b4d8", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    # Mésosphère (50-85 km)
    x = [-10 1010; -10 1010]
    y = [1010 1010; 1010 1010]
    z = [50 50; 85 85]
    plot_surface(x,y,z, color="#0096c7", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    # Thermosphère (85-500 km)
    x = [-10 1010; -10 1010]
    y = [1010 1010; 1010 1010]
    z = [85 85; 500 500]
    plot_surface(x,y,z, color="#0077b6", alpha=0.9, shade=false, zorder=-10, clip_on=false)

   ###### Fond n°2 ######

    # Troposphère (15 km)
    x = [-10 -10; -10 -10]
    y = [-10 1010; -10 1010]
    z = [0 0; 15 15]
    plot_surface(x,y,z, color="#48cae4", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    # Stratosphère (15-50 km)
    x = [-10 -10; -10 -10]
    y = [-10 1010; -10 1010]
    z = [15 15; 50 50]
    plot_surface(x,y,z, color="#00b4d8", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    # Mésosphère (50-85 km)
    x = [-10 -10; -10 -10]
    y = [-10 1010; -10 1010]
    z = [50 50; 85 85]
    plot_surface(x,y,z, color="#0096c7", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    # Thermosphère (85-500 km)
    x = [-10 -10; -10 -10]
    y = [-10 1010; -10 1010]
    z = [85 85; 500 500]
    plot_surface(x,y,z, color="#0077b6", alpha=0.9, shade=false, zorder=-10, clip_on=false)

    ax.plot([-10,-10], [-10,1010], [20,20], color="white", marker="", markersize=5, linestyle="--", alpha=0.7, zorder=15, clip_on=false)
    ax.plot([-10,1010], [1010,1010], [20,20], color="white", marker="", markersize=5, linestyle="--", alpha=0.7, zorder=15, clip_on=false)

    text3D(x=0, y=-90, z=15, s="15 km.", zdir="y", color="white", bbox=Dict("facecolor"=>"#6096ba", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=15, clip_on=false)
    text3D(x=0, y=-90, z=20, s="20 km.", zdir="y", color="white", bbox=Dict("facecolor"=>"#6096ba", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=15, clip_on=false)
    text3D(x=0, y=10, z=7.5, s="Troposphère", zdir="y", color="black", zorder=15, clip_on=false)
    text3D(x=0, y=10, z=35, s="Stratosphère", zdir="y", color="black", zorder=15, clip_on=false)
    text3D(x=0, y=-90, z=50, s="50 km.", zdir="y", color="white", bbox=Dict("facecolor"=>"#6096ba", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=15, clip_on=false)
    text3D(x=0, y=10, z=67.5, s="Mésosphère", zdir="y", color="black", zorder=15, clip_on=false)
    text3D(x=0, y=-90, z=85, s="85 km.", zdir="y", color="white", bbox=Dict("facecolor"=>"#6096ba", "boxstyle"=>"round,pad=0.5", "alpha"=>0.8), zorder=15, clip_on=false)
    text3D(x=0, y=10, z=100, s="Thermosphère", zdir="y", color="black", zorder=15, clip_on=false)
end

function generation_etoiles(ax)
    # Seed : chaîne de caractères vers suite d'entiers (conversion via code ASCII)
    chaine = "orion"
    conversion = string.([Int(c) for c in chaine])
    graine = ""
    for c in conversion
        graine *= c
    end
    graine = parse(Int, graine)
    Random.seed!(graine)

    # Fond 1
    x = rand([-1000:-1000;], 800)
    y = rand([-2000:2000;], 800)
    z = rand([-200:400;], 800)
    for g in 0:15
        ax.plot(x[1+50*g:50*(g+1)], y[1+50*g:50*(g+1)], z[1+50*g:50*(g+1)], color="white", marker=".", linestyle="", alpha=rand([0.4:0.1:1;]), markersize=alpha=rand([0.5:0.1:1;]), zorder=0, clip_on=false)
    end

    # Fond 2
    x = rand([-1000:2000;], 800)
    y = rand([2000:2000;], 800)
    z = rand([-200:400;], 800)
    for g in 0:15
        ax.plot(x[1+50*g:50*(g+1)], y[1+50*g:50*(g+1)], z[1+50*g:50*(g+1)], color="white", marker=".", linestyle="", alpha=rand([0.4:0.1:1;]), markersize=alpha=rand([0.5:0.1:1;]), zorder=0, clip_on=false)
    end

    # Soleil
    #ax.plot([-800], [-100], [70], color="#a10702", marker="o", markersize=105, linestyle="", alpha=0.5, zorder=-1)
    #ax.plot([-800], [-100], [70], color="#ff6600", marker="o", markersize=100, linestyle="", alpha=1, zorder=0)
end

function plot_base(coord_base)
    x = coord_base.x
    y = coord_base.y
    # Murs
    Z = [trunc.([0.05*i, 0.05*i], digits=2) for i in 0:20]
    cpt = 1
    for z in Z
        if z == [0.0,0.0]
            # marker="", markersize=0.2, markeredgecolor="black",
            plot([x-10,x+10], [y+10,y+10], z, color="dimgray", linewidth=1, zorder=16, clip_on=false) # mur haut
            plot([x-10,x+10], [y-10,y-10], z, color="dimgray", linewidth=1, zorder=16, clip_on=false) # mur bas
            plot([x-10,x-10], [y-10,y+10], z, color="dimgray", linewidth=1, zorder=16, clip_on=false) # mur gauche
            plot([x+10,x+10], [y-10,y+10], z, color="dimgray", linewidth=1, zorder=16, clip_on=false) # mur droite
        elseif z == [1.0, 1.0]
            plot([x-10,x+10], [y+10,y+10], z, color="dimgray", linewidth=1, zorder=18, clip_on=false) # mur haut
            plot([x-10,x+10], [y-10,y-10], z, color="dimgray", linewidth=1, zorder=18, clip_on=false) # mur bas
            plot([x-10,x-10], [y-10,y+10], z, color="dimgray", linewidth=1, zorder=18, clip_on=false) # mur gauche
            plot([x+10,x+10], [y-10,y+10], z, color="dimgray", linewidth=1, zorder=18, clip_on=false) # mur droite
        else
            if cpt % 2 == 0
                couleur = "silver"
            else
                couleur = "snow"
            end
            plot([x-10,x+10], [y+10,y+10], z, color=couleur, linewidth=1, zorder=16, clip_on=false) # mur haut
            plot([x-10,x+10], [y-10,y-10], z, color=couleur, linewidth=1, zorder=16, clip_on=false) # mur bas
            plot([x-10,x-10], [y-10,y+10], z, color=couleur, linewidth=1, zorder=16, clip_on=false) # mur gauche
            plot([x+10,x+10], [y-10,y+10], z, color=couleur, linewidth=1, zorder=16, clip_on=false) # mur droite
            cpt += 1
        end
    end
    # Toit
    cpt = 1
    for i in x-9.9:0.1:x+9.9
        if cpt % 2 == 0
            couleur = "silver"
        else
            couleur = "snow"
        end
        plot([i,i], [y-10,y+10], [1.0,1.0], color=couleur, marker="", linewidth=1, zorder=17, clip_on=false) # toit
        cpt += 1
    end
end

function plot_position_valide(x,y)
    for i in x-2:1:x+2
        plot([i,i], [y-3,y+3], [0,0], color="darkorange", marker="", alpha=0.2, linewidth=2, zorder=16, clip_on=false) # toit
    end
end

function plot_position_invalide(x,y)
    for i in x-2:1:x+2
        plot([i,i], [y-3,y+3], [0,0], color="red", marker="", alpha=0.2, linewidth=2, zorder=16, clip_on=false) # toit
    end
end

########### Visualisation d'une solution complète ###########

function visualisation(instance, solution)
    # Données (contexte & scénario)
    nT = instance.contexte.nT
    nH = instance.contexte.nH
    nE = instance.contexte.nE
    nC = instance.contexte.nC
    nU = instance.contexte.nU
    nB = instance.contexte.nB
    nR = instance.contexte.nR
    nRc = instance.contexte.nRc
    A = instance.contexte.A
    R = instance.contexte.R
    Cu = instance.contexte.Cu
    Cb = instance.contexte.Cb
    Eb = instance.contexte.Eb
    E = instance.contexte.E
    Ebool = instance.contexte.Ebool
    W = instance.contexte.W
    P = instance.contexte.P
    w = instance.contexte.w
    p = instance.contexte.p
    s = instance.contexte.s

    deplacements = instance.scenario.deplacements

    # Solution
    z = solution.z
    nb_unites_non_couvertes = solution.nb_unites_non_couvertes
    deploiement_haps = solution.deploiement_haps
    placement_relais = solution.placement_relais
    unites_couvertes = solution.unites_couvertes
    bases_couvertes = solution.bases_couvertes
    transmission_haps = solution.transmission_haps

    println("\n========== Visualisation de la solution ==========")

    # Liste de couleurs prédéfinies pour les différents types d'unités
    liste_couleurs = ["red", "yellow", "blue", "lime", "darkorchid", "fuchsia", "snow", "cyan"]

    customisation_pyplot()

    # Terrain (3-dimensions)
    fig, ax, affichage_nb_images_seconde, affichage_etat_courant = initialisation_base(instance, liste_couleurs, 1)
    affichage_nb_unites_non_couvertes = initialisation_suite(ax, instance, solution, liste_couleurs)
    show()

    # Comète (trajectoire)
    XC = [i for i in range(-200, -200, length=nT)]
    YC = [i for i in range(-200, 1200, length=nT)]
    ZC = [i for i in range(40, 10, length=nT)]
    trace_comete = trunc(Int, 0.025*nT)

    # On affiche les positions de déploiement envisageables
    for position in 1:nE
        coord = E[position]
        if Ebool[position] == true # position active
            ax.plot([coord.x], [coord.y], [0], color="darkorange", marker="s", linestyle="", markersize=2.5, alpha=0.5, zorder=16, clip_on=false)
            #plot_position_valide(coord.x, coord.y)
        else
            ax.plot([coord.x], [coord.y], [0], color="red", marker="s", linestyle="", markersize=2.5, alpha=0.5, zorder=16, clip_on=false)
            #plot_position_invalide(coord.x, coord.y)
        end
    end

    # On affiche les bases et les liens de couverture
    for b in 1:nB
        coord_base = Eb[b]
        plot_base(coord_base)
        #ax.plot([coord_base.x], [coord_base.y], [0], color="black", marker="p", linestyle="", markersize=8, alpha=1, zorder=16, clip_on=false)
        types_communication = Dict{Int, Vector{Int}}()
        for c in Cb[b]
            if length(bases_couvertes[b][c]) != 0 # base b couverte par le type c (dont elle dispose)
                for haps in bases_couvertes[b][c]
                    #ax.plot([coord_base.x, E[deploiement_haps[haps]].x], [coord_base.y, E[deploiement_haps[haps]].y], [1, A], color=liste_couleurs[c], linestyle="--",  linewidth=0.6, zorder=18, clip_on=false)
                    if haskey(types_communication, haps)
                        push!(types_communication[haps], c)
                    else
                        types_communication[haps] = [c]
                    end
                end
            end
        end
        # Segmentation des liens selon le nombre de types de communication couvrant la base b (dans un même HAPS)
        for (haps, liste_c) in types_communication
            vect = [E[deploiement_haps[haps]].x - coord_base.x, E[deploiement_haps[haps]].y - coord_base.y, A - 1]
            point_1 = [coord_base.x, coord_base.y, 1]
            for i in 1:length(liste_c)
                point_2 = point_1 .+ (1/length(liste_c)) .* vect
                ax.plot([point_1[1], point_2[1]], [point_1[2], point_2[2]], [point_1[3], point_2[3]], color=liste_couleurs[liste_c[i]], linestyle="--",  linewidth=1, zorder=18, clip_on=false)
                point_1 = point_2
            end
        end
    end

    # On affiche les HAPS déployés
    for position_haps in deploiement_haps
        if position_haps != -1
            ax.plot([E[position_haps].x], [E[position_haps].y], [A], color="darkorange", marker="o", markeredgecolor="black", linestyle="", markersize=14, zorder=25, clip_on=false)
            ax.plot([E[position_haps].x, E[position_haps].x], [E[position_haps].y, E[position_haps].y], [0, A], color="black", linestyle="dotted", zorder=19, clip_on=false, linewidth=1)
        end
    end

    # On affiche les arcs de transmission
    transmissions_haps_haps = Dict{Tuple{Int, Int}, Vector{Int}}()
    if length(transmission_haps) != 0
        for haps in 1:nH
            position_haps = deploiement_haps[haps]
            if position_haps != -1
                for c in 1:nC
                    for haps_bis in transmission_haps[haps][c]
                        position_haps_bis = deploiement_haps[haps_bis]
                        #ax.plot([E[position_haps].x, E[position_haps_bis].x], [E[position_haps].y, E[position_haps_bis].y], [A, A], color=liste_couleurs[c], linestyle="dotted", zorder=20, clip_on=false, linewidth=1)
                        if haskey(transmissions_haps_haps, (haps,haps_bis))
                            push!(transmissions_haps_haps[(haps, haps_bis)], c)
                        else
                            transmissions_haps_haps[(haps, haps_bis)] = [c]
                        end
                    end
                end
            end
        end
        couples = [(haps,haps_bis) for haps in 1:nH, haps_bis in 1:nH if haps < haps_bis]
        for (haps, haps_bis) in couples
            position_haps = deploiement_haps[haps]
            position_haps_bis = deploiement_haps[haps_bis]
            if haskey(transmissions_haps_haps, (haps,haps_bis))
                #println(" - $haps --> $(haps_bis) $(transmissions_haps_haps[(haps, haps_bis)])")
                vect = [E[position_haps_bis].x - E[position_haps].x, E[position_haps_bis].y - E[position_haps].y]
                point_1 = [E[position_haps].x, E[position_haps].y]
                for c in transmissions_haps_haps[(haps,haps_bis)]
                    point_2 = point_1 .+ (1/length(transmissions_haps_haps[(haps,haps_bis)])) .* vect
                    ax.plot([point_1[1], point_2[1]], [point_1[2], point_2[2]], [A + 0.5, A + 0.5], color=liste_couleurs[c], linestyle="dotted", zorder=18, clip_on=false, linewidth=1.5)
                    ax.plot([(point_1[1]+point_2[1])/2], [(point_1[2]+point_2[2])/2], [A + 0.5], color=liste_couleurs[c], zorder=19, marker=">", markersize=4)
                    point_1 = point_2
                    #ax.quiver([point_1[1]], [point_1[2]], [A + 0.5], [vect[1]], [vect[2]], [0], color=liste_couleurs[c], linestyle="dotted", zorder=18, clip_on=false, headwidth=1)
                end
            end
            if haskey(transmissions_haps_haps, (haps_bis,haps))
                #println(" - $(haps_bis) --> $haps $(transmissions_haps_haps[(haps_bis, haps)])")
                vect = [E[position_haps].x - E[position_haps_bis].x, E[position_haps].y - E[position_haps_bis].y]
                point_1 = [E[position_haps_bis].x, E[position_haps_bis].y]
                for c in transmissions_haps_haps[(haps_bis,haps)]
                    point_2 = point_1 .+ (1/length(transmissions_haps_haps[(haps_bis,haps)])) .* vect
                    ax.plot([point_1[1], point_2[1]], [point_1[2], point_2[2]], [A - 0.5, A - 0.5], color=liste_couleurs[c], linestyle="dotted", zorder=18, clip_on=false, linewidth=1.5)
                    ax.plot([(point_1[1]+point_2[1])/2], [(point_1[2]+point_2[2])/2], [A - 0.5], color=liste_couleurs[c], zorder=19, marker="<", markersize=4)
                    point_1 = point_2
                end
            end
        end
    end

    # On affiche la portée de chaque relais déployés dans chaque HAPS (cercle sur le plan xy)
    for c in 1:nC
        for r in R[c]
            if placement_relais[r] != -1
                nb_points = 40
                theta = range(0, 2 * pi, length=nb_points)
                rayon = sqrt((s[r])^2 - (A)^2)
                X = E[deploiement_haps[placement_relais[r]]].x .+ rayon*cos.(theta)
                Y = E[deploiement_haps[placement_relais[r]]].y .+ rayon*sin.(theta)
                points_a_conserver = [[]]
                for i in 1:nb_points
                    if (X[i] >= -10 && X[i] <= 1010) && (Y[i] >= -10 && Y[i] <= 1010)
                        push!(points_a_conserver[end], i)
                    else
                        push!(points_a_conserver, [])
                    end
                end
                for k in 1:length(points_a_conserver)
                    ax.plot([X[i] for i in points_a_conserver[k]], [Y[i] for i in points_a_conserver[k]], [0 for i in points_a_conserver[k]], color=liste_couleurs[c], marker=".", zorder=15, markersize=4, alpha=0.3, clip_on=false)
                end
            end
        end
    end

    sauvegarde_animation = choix_binaire("\n --> Souhaitez-vous enregistrer les images (o/n) ? ")
    if sauvegarde_animation == "o"
        nom_dossier = string(today())*"_"*Dates.format(now(), "HH:MM:SS")
        run(`mkdir animations/$(nom_dossier)`)
    end

    println("\n Liste des actions disponibles : ")
    println(" ------------------------------- ")
    println("    'd' : démarrer la visualisation")
    println("    'i' : interrompre la visualisation")
    println("    'p' : mettre en pause/reprendre la visualisation")
    println("    'a' : accélerer la visualisation")
    println("    'r' : ralentir la visualisation")
    println("    't' : activer/désactiver la trace des unités")
    println("\n /!\\ Cliquez sur le graphique, puis sur les touches  /!\\ ")

    # On affiche les différentes unités terrestres au fil du temps
    demarrer = [false]
    continuer = [true]
    pause = [false]
    trace_unite = [true]
    nb_images_seconde = [10] # Une image = un état
    cid_bis = fig.canvas.mpl_connect("key_press_event", event -> on_press(event, demarrer, continuer, pause, nb_images_seconde, trace_unite))
    t = 1
    taille_trace = 6
    trace = Vector{Vector{Point}}(undef, nU)
    trace = [Vector{Point}(undef, taille_trace) for i in 1:length(trace)]
    while t <= nT && continuer[1] == true
        affichage_nb_images_seconde.set_text("Nombre d'images par seconde = $(nb_images_seconde[1])")
        affichage_etat_courant.set_text("État : $t/$(nT)")
        affichage_nb_unites_non_couvertes.set_text("Nombre d'unités non couvertes = $(nb_unites_non_couvertes[t])/$(nU) ($(trunc(((nb_unites_non_couvertes[t])/(nU)*100), digits=2))%)")
        liste_points = []
        # Comète
        push!(liste_points, ax.plot([XC[t]], [YC[t]], [ZC[t]], color="#caf0f8", alpha=1, marker="o", markersize=5, zorder=0, clip_on=false))
        t <= trace_comete ? push!(liste_points, ax.plot(XC[1:t], YC[1:t], ZC[1:t], color="#caf0f8", alpha=0.4, marker="", linewidth=1.5, zorder=-10, clip_on=false)) : push!(liste_points, ax.plot(XC[t-trace_comete:t], YC[t-trace_comete:t], ZC[t-trace_comete:t], color="#caf0f8", alpha=0.4, marker="", linewidth=1.5, zorder=-10, clip_on=false))
        # Unités
        for u in 1:nU
            coordUnite = deplacements[u][t]
            # On utilise ici plot et non ax.scatter3D car ce dernier posait des soucis de superposition avec le terrain
            push!(liste_points, ax.plot([coordUnite.x], [coordUnite.y], [0], color="#ffe8d6", marker="^", markeredgecolor="black", markersize=7, alpha=1.0, linestyle="", zorder=21, clip_on=false))
            for i in 0:taille_trace-2
                trace[u][taille_trace - i] = trace[u][taille_trace - i - 1]
            end
            trace[u][1] = coordUnite
            if trace_unite[1] == true
                if t >= taille_trace
                    push!(liste_points, ax.plot([trace[u][i].x for i in 1:taille_trace], [trace[u][i].y for i in 1:taille_trace], [0], color="#ffe8d6", linewidth=0.9, marker="", alpha=0.6, zorder=20, clip_on=false))
                else
                    push!(liste_points, ax.plot([trace[u][i].x for i in 1:t], [trace[u][i].y for i in 1:t], [0], color="#ffe8d6", marker="", linewidth=0.9, alpha=0.6, zorder=20, clip_on=false))
                end
             end
             types_communication = Dict{Int, Vector{Int}}()
             for c in Cu[u]
                 if length(unites_couvertes[u][t][c]) != 0 # unité u couverte au temps t par le type c (dont elle dispose)
                     for haps in unites_couvertes[u][t][c]
                         #push!(liste_points, ax.plot([coordUnite.x, E[deploiement_haps[haps]].x], [coordUnite.y, E[deploiement_haps[haps]].y], [0, A], color=liste_couleurs[c], linestyle="--",  linewidth=0.6, zorder=19, clip_on=false))
                        if haskey(types_communication, haps)
                           push!(types_communication[haps], c)
                        else
                            types_communication[haps] = [c]
                        end
                    end
                 end
            end
            # Segmentation des liens selon le nombre de types de communication couvrant la base b (dans un même HAPS)
            for (haps, liste_c) in types_communication
                vect = [E[deploiement_haps[haps]].x - coordUnite.x, E[deploiement_haps[haps]].y - coordUnite.y, A]
                point_1 = [coordUnite.x, coordUnite.y, 0]
                for i in 1:length(liste_c)
                    point_2 = point_1 .+ (1/length(liste_c)) .* vect
                    push!(liste_points, ax.plot([point_1[1], point_2[1]], [point_1[2], point_2[2]], [point_1[3], point_2[3]], color=liste_couleurs[liste_c[i]], linestyle="--",  linewidth=1, zorder=18, clip_on=false))
                    point_1 = point_2
                end
            end
        end
        if t == 1 # Affichage des unités à l'état 1 et on attend avant de démarrer la visualisation
            while demarrer[1] == false
                sleep(0.1)
            end
        end
        if sauvegarde_animation == "o" # Enregistrement de chaque frame (on laisse un temps pour quitter l'animation si nécessaire)
            fig.savefig("animations/$(nom_dossier)/image_"*string(t), format="png")
            sleep(1)
        end
        while pause[1] == true # Pause à durée indéterminée
            sleep(0.1)
        end
        sleep(1/nb_images_seconde[1])
        if t != nT
            for point in liste_points
                point[1].remove()
            end
        else
            print("\n --> Appuyez sur la touche \"Entrée\" pour terminer ")
            readline()
        end
        t += 1
    end
    fig.canvas.mpl_disconnect(cid_bis)
    close()
    if sauvegarde_animation == "o"
        println("\n /!\\ Les images sont disponibles dans le dossier \"animations/$(nom_dossier)/\" /!\\ ")
    end
    println("\n========== Visualisation de la solution ==========")
end

########### Visualisation du scénario uniquement ###########

function visualisation_scenario(instance)
    # Données (contexte & scénario)
    nT = instance.contexte.nT
    nH = instance.contexte.nH
    nE = instance.contexte.nE
    nC = instance.contexte.nC
    nU = instance.contexte.nU
    nB = instance.contexte.nB
    nR = instance.contexte.nR
    nRc = instance.contexte.nRc
    A = instance.contexte.A
    R = instance.contexte.R
    Cu = instance.contexte.Cu
    Cb = instance.contexte.Cb
    Eb = instance.contexte.Eb
    E = instance.contexte.E
    Ebool = instance.contexte.Ebool
    W = instance.contexte.W
    P = instance.contexte.P
    w = instance.contexte.w
    p = instance.contexte.p
    s = instance.contexte.s

    deplacements = instance.scenario.deplacements

    println("\n========== Visualisation du scénario ==========")

    # Liste de couleurs prédéfinies pour les différents types d'unités
    #liste_couleurs = ["red", "navy", "deeppink", "darkorange", "limegreen", "coral", "saddlebrown", "gold", "darkkhaki", "olive", "darkolivegreen", "lightseagreen", "paleturquoise", "deepskyblue",  "blue", "darkviolet", "magenta", "darkseagreen", "mediumpurple"]
    liste_couleurs = ["red", "yellow", "blue", "cyan", "lime", "darkviolet", "fuchsia", "slategrey", "snow"]

    customisation_pyplot()

    # Terrain (3-dimensions)
    fig, ax, affichage_nb_images_seconde, affichage_etat_courant = initialisation_base(instance, liste_couleurs, 2)

    # Comète (trajectoire)
    XC = [i for i in range(-200, -200, length=nT)]
    YC = [i for i in range(-200, 1200, length=nT)]
    ZC = [i for i in range(40, 10, length=nT)]
    trace_comete = trunc(Int, 0.03*nT)

    # Affichage
    show()

    sauvegarde_animation = choix_binaire("\n --> Souhaitez-vous enregistrer les images (o/n) ? ")
    if sauvegarde_animation == "o"
        nom_dossier = string(today())*"_"*Dates.format(now(), "HH:MM:SS")
        run(`mkdir animations/$(nom_dossier)`)
    end

    println("\n Liste des actions disponibles : ")
    println(" ------------------------------- ")
    println("    'd' : démarrer la visualisation")
    println("    'i' : interrompre la visualisation")
    println("    'p' : mettre en pause/reprendre la visualisation")
    println("    'a' : accélerer la visualisation")
    println("    'r' : ralentir la visualisation")
    println("    't' : activer/désactiver la trace des unités")
    println("\n /!\\ Cliquez sur le graphique, puis sur les touches  /!\\ ")

    # On affiche les différentes unités terrestres au fil du temps
    demarrer = [false]
    continuer = [true]
    pause = [false]
    trace_unite = [true]
    nb_images_seconde = [10] # Une image = un état
    cid_bis = fig.canvas.mpl_connect("key_press_event", event -> on_press(event, demarrer, continuer, pause, nb_images_seconde, trace_unite))
    t = 1
    taille_trace = 6
    trace = Vector{Vector{Point}}(undef, nU)
    trace = [Vector{Point}(undef, taille_trace) for i in 1:length(trace)]
    while t <= nT && continuer[1] == true
        affichage_nb_images_seconde.set_text("Nombre d'images par seconde = $(nb_images_seconde[1])")
        affichage_etat_courant.set_text("État = $t/$(nT)")
        liste_points = []
        # Comète
        push!(liste_points, ax.plot([XC[t]], [YC[t]], [ZC[t]], color="#caf0f8", alpha=1, marker="o", markersize=7, zorder=0, clip_on=false))
        t <= trace_comete ? push!(liste_points, ax.plot(XC[1:t], YC[1:t], ZC[1:t], color="#caf0f8", alpha=0.4, marker="", linewidth=4, zorder=-10, clip_on=false)) : push!(liste_points, ax.plot(XC[t-trace_comete:t], YC[t-trace_comete:t], ZC[t-trace_comete:t], color="#caf0f8", alpha=0.4, marker="", linewidth=4, zorder=-10, clip_on=false))
        # Unités
        for u in 1:nU
            coordUnite = deplacements[u][t]
            # On utilise ici plot et non scatter3D car ce dernier posait des soucis de superposition avec le terrain
            push!(liste_points, ax.plot([coordUnite.x], [coordUnite.y], [0], color="#ffe8d6", marker="^", markeredgecolor="black", markersize=7, alpha=1.0, linestyle="", zorder=15, clip_on=false))
            for i in 0:taille_trace-2
                trace[u][taille_trace - i] = trace[u][taille_trace - i - 1]
            end
            trace[u][1] = coordUnite
            if trace_unite[1] == true
                if t >= taille_trace
                    push!(liste_points, ax.plot([trace[u][i].x for i in 1:taille_trace], [trace[u][i].y for i in 1:taille_trace], [0], color="#ffe8d6", marker="", alpha=0.5, zorder=20, linewidth=0.9, clip_on=false))
                else
                    push!(liste_points, ax.plot([trace[u][i].x for i in 1:t], [trace[u][i].y for i in 1:t], [0], color="#ffe8d6", marker="", alpha=0.5, zorder=20, linewidth=0.9, clip_on=false))
                end
            end
        end
        if t == 1 # Affichage les unités à l'état 1 et on attend avant de démarrer la visualisation
            while demarrer[1] == false
                sleep(0.1)
            end
        end
        if sauvegarde_animation == "o" # Enregistrement de chaque frame (on laisse un temps pour quitter l'animation si nécessaire)
            fig.savefig("animations/$(nom_dossier)/image_"*string(t), format="png")
            sleep(1)
        end
        while pause[1] == true # Pause à durée indéterminée
            sleep(0.1)
        end
        sleep(1/nb_images_seconde[1])
        if t != nT
            for point in liste_points
                point[1].remove()
            end
        else
            print("\n --> Appuyez sur la touche \"Entrée\" pour terminer ")
            readline()
        end
        t += 1
    end
    fig.canvas.mpl_disconnect(cid_bis)
    close()
    if sauvegarde_animation == "o"
        println("\n /!\\ Les images sont disponibles dans le dossier \"animations/$(nom_dossier)/\" /!\\ ")
    end
    println("\n========== Visualisation du scénario ==========")
end
