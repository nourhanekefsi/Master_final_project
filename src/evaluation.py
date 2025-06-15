import numpy as np
from collections import defaultdict
import numpy as np
from math import sqrt

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
#=============================================================================================================================
                                 # Metriques d'evaluation
#================================================================================================================================
import numpy as np
from math import sqrt

def overlap_score(pred, real):
    return len(set(pred) & set(real)) / sqrt(len(pred) * len(real)) if pred and real else 0

def jaccard_index(pred, real):
    return len(set(pred) & set(real)) / len(set(pred) | set(real)) if set(pred) | set(real) else 0

def compute_metrics(predicted_complexes, real_complexes, threshold=0.2):
    m = len(predicted_complexes)
    n = len(real_complexes)

    # Intersection max for PPV
    ppv_numerator = sum(max(len(set(pred) & set(real)) for real in real_complexes) for pred in predicted_complexes)
    ppv_denominator = sum(len(pred) for pred in predicted_complexes)
    ppv = ppv_numerator / ppv_denominator if ppv_denominator > 0 else 0

    # Intersection max for Sn
    sn_numerator = sum(max(len(set(pred) & set(real)) for pred in predicted_complexes) for real in real_complexes)
    sn_denominator = sum(len(real) for real in real_complexes)
    sn = sn_numerator / sn_denominator if sn_denominator > 0 else 0

    # F-measure
    f_measure = (2 * ppv * sn) / (ppv + sn) if (ppv + sn) > 0 else 0

    # Accuracy (as sqrt(Sn * PPV))
    accuracy = sqrt(sn * ppv)

    # MMR
    mmr = np.mean([max(overlap_score(pred, real) for real in real_complexes) for pred in predicted_complexes]) if predicted_complexes else 0

    # Jaccard
    jaccard = np.mean([max(jaccard_index(pred, real) for real in real_complexes) for pred in predicted_complexes]) if predicted_complexes else 0

    # Covered Rate
    matched_reals = sum(
        any(overlap_score(pred, real) >= threshold for pred in predicted_complexes)
        for real in real_complexes
    )
    covered_rate = matched_reals / n if n > 0 else 0

    # Score total
    score_total = f_measure + covered_rate + mmr + jaccard + accuracy

    return {
        "PPV": ppv,
        "Recall (Sn)": sn,
        "F-mesure": f_measure,
        "Accuracy": accuracy,
        "MMR": mmr,
        "Jaccard": jaccard,
        "Covered Rate": covered_rate,
        "Score Total": score_total
    }

