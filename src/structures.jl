struct Point
    x::Float64
    y::Float64
end

struct Contexte
    nT::Int # Nombre d'états
    nH::Int # Nombre de HAPS
    nE::Int # Nombre de positions envisageables pour le déploiement d'un HAPS
    nC::Int # Nombre de types de communication
    nU::Int # Nombre d'unités
    nB::Int # Nombre de bases
    nR::Int # Nombre de relais
    nRc::Vector{Int} # Nombre de relais par type | nRc[c] : nombre de relais de type c
    A::Int # Altitude
    R::Vector{Vector{Int}} # Ensemble des relais par type de communication | R[c] : relais de type c
    Cu::Vector{Vector{Int}} # Ensemble des types de communication par unité | Cu[u] : types de communication de l'unité u
    Pu::Vector{Vector{Int}} # Ensemble des portée par types de communication et par unité | Pu[u] : portées de u, i ème portée --> i ème moyen de com dans Cu[u]
    Cb::Vector{Vector{Int}} # Ensemble des types de communication par base | Cb[b] : types de communication de la base b
    Eb::Vector{Point} # Positions des bases | E[b] : position de la base b (point)
    E::Vector{Point} # Liste des positions envisageables pour le déploiement | E[i] : position i (point)
    Ebool::Vector{Bool} # Liste des positions actives | Ebool[i] = 0 : position inactive (non-couvrante)
    W::Vector{Int} # Poids maximal de chacun des HAPS | W[h] : poids HAPS h
    P::Vector{Int} # Puissance maximale de chacun des HAPS | P[h] : puissance HAPS h
    w::Vector{Int} # Poids des différents relais | w[r] : poids relais r
    p::Vector{Int} # Puissance des différents relais | p[r] : puissance relais r
    s::Vector{Int} # Portée (seuil) des différents relais | s[r] : portée (seuil) relais r
end

struct Scenario
    deplacements::Vector{Vector{Point}} # Déplacements des unités terrestres sur un horizon temporel
    # deplacements[u] : déplacements de l'unité u
    # deplacements[u][t] : position de l'unité u au temps t
end

struct Instance
    contexte::Contexte
    scenario::Scenario
end

struct Solution
    statut::Int # 0 : pas de solution, 1 : solution disponible
    z::Int # Nombre d'unités non couvertes
    nb_unites_non_couvertes::Vector{Int} # Nombre d'unités non couvertes à chaque instant
    deploiement_haps::Vector{Int} # Liste des positions sur lesquelles sont déployées les HAPS (si déployés) | deploiement_haps[h] : position sur laquelle est déployé le HAPS h (si déployé)
    placement_relais::Vector{Int} # Liste des HAPS sur lesquels sont placés les relais (si placés) | placement_relais[r] : HAPS sur lequel est placé le relais r (si placé)
    unites_couvertes::Vector{Vector{Vector{Vector{Int}}}} # Liste des HAPS couvrant les unités à chaque instant (si couvertes) | unites_couvertes[u][t][c] : ensemble des HAPS couvrant l'unité u au temps t pour le type de communication c (dont dispose l'unité u)
    bases_couvertes::Vector{Vector{Vector{Int}}} # Liste des HAPS couvrant les bases (si couvertes) | bases_couvertes[b][c] : ensemble des HAPS couvrant la base b pour le type de communication c (dont dispose la base b)
    transmission_haps::Vector{Vector{Vector{Int}}} # Liste des HAPS vers lesquels peut transmettre un HAPS via un certain type de communication | transmission_haps[h][c] : HAPS vers lequel le HAPS h peut transmettre via le type de communication c
end
