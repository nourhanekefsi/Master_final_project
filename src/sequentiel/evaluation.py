import numpy as np
from collections import defaultdict

# Fonctions pour évaluer un seul cluster
def cohesiveness(cluster, graph):
    """
    Calcule la cohésion d'un cluster.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI sous forme de dictionnaire {nœud: {voisins: poids}}
    :return: Score de cohésion
    """
    W_in = 0
    W_out = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor in cluster:
                W_in += weight
            else:
                W_out += weight
    return W_in / (W_in + W_out) if (W_in + W_out) > 0 else 0

def density(cluster, graph):
    """
    Calcule la densité d'un cluster.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score de densité
    """
    W_in = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor in cluster:
                W_in += weight
    n = len(cluster)
    return (2 * W_in) / (n * (n - 1)) if n > 1 else 0

def AIEW(cluster, graph):
    """
    Calcule le poids moyen des arêtes internes (Average Inner Edge Weight).
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score AIEW
    """
    W_in = 0
    E_C = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor in cluster:
                W_in += weight
                E_C += 1
    return W_in / E_C if E_C > 0 else 0

def ABEW(cluster, graph):
    """
    Calcule le poids moyen des arêtes frontalières (Average Border Edge Weight).
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score ABEW
    """
    W_out = 0
    BE_C = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor not in cluster:
                W_out += weight
                BE_C += 1
    return W_out / BE_C if BE_C > 0 else 0

def AWM(cluster, graph):
    """
    Calcule la modularité pondérée moyenne (Average Weighted Modularity).
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score AWM
    """
    aiew = AIEW(cluster, graph)
    abew = ABEW(cluster, graph)
    return aiew / (aiew + abew) if (aiew + abew) > 0 else 0

def FF(cluster, graph):
    """
    Fonction de fitness pour un seul cluster.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score FF combiné
    """
    return (
        cohesiveness(cluster, graph) +
        density(cluster, graph) +
        AIEW(cluster, graph) -
        ABEW(cluster, graph) +
        AWM(cluster, graph)
    )

# Fonction pour évaluer un individu (ensemble de clusters)
def FS_fitness(individual, graph):
    """
    Calcule le fitness total d'un individu (ensemble de clusters).
    :param individual: Liste de clusters
    :param graph: Graphe PPI
    :return: Score FS_fitness
    """
    return sum(FF(cluster, graph) for cluster in individual)

# Fonction d'optimisation locale pour un cluster
def local_optimization(cluster, graph, max_iter=20):
    """
    Optimise un cluster localement en ajoutant/supprimant des nœuds.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :param max_iter: Nombre maximum d'itérations
    :return: Cluster optimisé
    """
    optimized_cluster = set(cluster)
    changed = True
    iteration = 0

    while changed and iteration < max_iter:
        changed = False
        current_ff = FF(optimized_cluster, graph)

        # Étape 1: Suppression des nœuds internes
        inner_nodes = {
            node for node in optimized_cluster 
            if any(neighbor not in optimized_cluster for neighbor in graph.get(node, {}))
        }
        if len(inner_nodes) > 0:
            worst_node = min(inner_nodes, key=lambda x: FF(optimized_cluster - {x}, graph))
            new_cluster = optimized_cluster - {worst_node}
            if FF(new_cluster, graph) > current_ff:
                optimized_cluster = new_cluster
                changed = True
                continue

        # Étape 2: Ajout des nœuds frontaliers
        boundary_nodes = set()
        for node in optimized_cluster:
            for neighbor in graph.get(node, {}):
                if neighbor not in optimized_cluster:
                    boundary_nodes.add(neighbor)
        
        if boundary_nodes:
            best_node = max(boundary_nodes, key=lambda x: FF(optimized_cluster | {x}, graph))
            new_cluster = optimized_cluster | {best_node}
            if FF(new_cluster, graph) > current_ff:
                optimized_cluster = new_cluster
                changed = True

        iteration += 1

    return optimized_cluster

# Métriques de performance
def overlap_score(detected, known):
    """
    Calcule le score de chevauchement entre un cluster détecté et un cluster connu.
    :param detected: Cluster détecté (ensemble de protéines)
    :param known: Cluster connu (ensemble de protéines)
    :return: Score de chevauchement
    """
    intersection = len(detected & known)
    return (intersection ** 2) / (len(detected) * len(known)) if len(detected) * len(known) > 0 else 0

def precision_recall_fmeasure(detected_complexes, known_complexes, threshold=0.2):
    """
    Calcule la précision, le rappel et le F-measure.
    :param detected_complexes: Liste des clusters détectés
    :param known_complexes: Liste des clusters connus
    :param threshold: Seuil pour considérer un match
    :return: (precision, recall, fmeasure)
    """
    TP = 0
    for detected in detected_complexes:
        for known in known_complexes:
            if overlap_score(set(detected), set(known)) >= threshold:
                TP += 1
                break
    
    precision = TP / len(detected_complexes) if len(detected_complexes) > 0 else 0
    recall = TP / len(known_complexes) if len(known_complexes) > 0 else 0
    fmeasure = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    return precision, recall, fmeasure

def coverage_rate(detected_complexes, known_complexes):
    """
    Calcule le taux de couverture.
    :param detected_complexes: Liste des clusters détectés
    :param known_complexes: Liste des clusters connus
    :return: Taux de couverture
    """
    covered_proteins = set()
    for known in known_complexes:
        for detected in detected_complexes:
            if len(set(detected) & set(known)) > 0:
                covered_proteins.update(known)
                break
    
    total_proteins = set().union(*known_complexes)
    return len(covered_proteins) / len(total_proteins) if len(total_proteins) > 0 else 0

def accuracy(detected_complexes, known_complexes, threshold=0.2):
    """
    Calcule l'accuracy (moyenne géométrique de la précision et du rappel).
    :param detected_complexes: Liste des clusters détectés
    :param known_complexes: Liste des clusters connus
    :return: Score d'accuracy
    """
    precision, recall, _ = precision_recall_fmeasure(detected_complexes, known_complexes, threshold)
    return np.sqrt(precision * recall)

def MMR(detected_complexes, known_complexes):
    """
    Calcule le Maximum Matching Ratio.
    :param detected_complexes: Liste des clusters détectés
    :param known_complexes: Liste des clusters connus
    :return: Score MMR
    """
    matched_pairs = set()
    total = 0
    
    for known in known_complexes:
        max_score = 0
        best_match = None
        for detected in detected_complexes:
            score = overlap_score(set(detected), set(known))
            if score > max_score:
                max_score = score
                best_match = tuple(detected)
        if best_match and best_match not in matched_pairs:
            matched_pairs.add(best_match)
            total += max_score
    
    return total / len(known_complexes) if len(known_complexes) > 0 else 0

def jaccard_index(detected_complexes, known_complexes):
    """
    Calcule l'indice de Jaccard.
    :param detected_complexes: Liste des clusters détectés
    :param known_complexes: Liste des clusters connus
    :return: Score de Jaccard
    """
    # Jaccard pour les clusters détectés
    jaccard_C = 0
    for detected in detected_complexes:
        max_overlap = 0
        for known in known_complexes:
            intersection = len(set(detected) & set(known))
            union = len(set(detected) | set(known))
            overlap = intersection / union if union > 0 else 0
            if overlap > max_overlap:
                max_overlap = overlap
        jaccard_C += max_overlap
    jaccard_C /= len(detected_complexes) if len(detected_complexes) > 0 else 0
    
    # Jaccard pour les clusters connus
    jaccard_G = 0
    for known in known_complexes:
        max_overlap = 0
        for detected in detected_complexes:
            intersection = len(set(detected) & set(known))
            union = len(set(detected) | set(known))
            overlap = intersection / union if union > 0 else 0
            if overlap > max_overlap:
                max_overlap = overlap
        jaccard_G += max_overlap
    jaccard_G /= len(known_complexes) if len(known_complexes) > 0 else 0
    
    # Jaccard combiné
    return (2 * jaccard_C * jaccard_G) / (jaccard_C + jaccard_G) if (jaccard_C + jaccard_G) > 0 else 0

def total_score(detected_complexes, known_complexes, threshold=0.2):
    """
    Calcule le score total combinant plusieurs métriques.
    :param detected_complexes: Liste des clusters détectés
    :param known_complexes: Liste des clusters connus
    :return: Score total
    """
    _, _, fmeasure = precision_recall_fmeasure(detected_complexes, known_complexes, threshold)
    cr = coverage_rate(detected_complexes, known_complexes)
    acc = accuracy(detected_complexes, known_complexes, threshold)
    mmr = MMR(detected_complexes, known_complexes)
    jaccard = jaccard_index(detected_complexes, known_complexes)
    
    return fmeasure + cr + acc + mmr + jaccard