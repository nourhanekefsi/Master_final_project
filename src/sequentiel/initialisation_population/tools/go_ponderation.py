import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import combinations
from typing import Dict, Set, Tuple

def load_ppi_network(ppi_path: str) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """Charge le réseau PPI et retourne les protéines + voisins."""
    ppi = pd.read_csv(ppi_path, sep='\t', usecols=['protein1', 'protein2'])
    proteins = set(ppi['protein1']).union(set(ppi['protein2']))
    
    # Dictionnaire de voisins
    neighbors = defaultdict(set)
    for p1, p2 in ppi.itertuples(index=False):
        neighbors[p1].add(p2)
        neighbors[p2].add(p1)
    
    return proteins, neighbors

def load_go_annotations(go_path: str, proteins: Set[str]) -> Dict[str, Set[str]]:
    """Charge les annotations GO pour les protéines du réseau."""
    go_data = pd.read_csv(
        go_path, 
        sep='\t', 
        header=None, 
        usecols=[0, 5], 
        names=['protein', 'go_term']
    )
    go_data = go_data[go_data['protein'].isin(proteins)]
    
    # Groupement par protéine
    go_annotations = defaultdict(set)
    for protein, go_term in go_data.itertuples(index=False):
        if pd.notna(go_term) and go_term.startswith('GO:'):
            go_annotations[protein].add(go_term)
    
    return go_annotations

def calculate_combined_weights(
    proteins: Set[str],
    neighbors: Dict[str, Set[str]],
    go_annotations: Dict[str, Set[str]]
) -> Dict[Tuple[str, str], float]:
    """Calcule les poids combinés GO + CN."""
    weights = {}
    protein_list = list(proteins)
    
    # Pré-calculer les annotations GO et le degré des protéines
    go_sets = {p: go_annotations.get(p, set()) for p in protein_list}
    degrees = {p: len(neighbors.get(p, set())) for p in protein_list}
    
    # Moyenne des annotations GO pour normalisation
    avg_go = np.mean([len(go) for go in go_sets.values()]) if go_sets else 1.0
    
    for i, (p1, p2) in enumerate(combinations(protein_list, 2)):
        if p2 not in neighbors.get(p1, set()):
            continue  # Ignorer les paires non interactives
        
        # 1. Calcul du score CN (Common Neighbors)
        cn_neighbors = neighbors.get(p1, set()) & neighbors.get(p2, set())
        cn_score = 0.0
        if degrees[p1] > 0 and degrees[p2] > 0:
            cn_score = len(cn_neighbors) ** 2 / (degrees[p1] * degrees[p2])
        
        # 2. Calcul du score GO
        go_p1 = go_sets.get(p1, set())
        go_p2 = go_sets.get(p2, set())
        go_intersection = go_p1 & go_p2
        go_score = 0.0
        if go_intersection:
            min_go = min(len(go_p1), len(go_p2))
            denom = max(min_go, avg_go)
            go_score = len(go_intersection) / denom
        
        # 3. Combinaison des scores (moyenne)
        combined_score = (go_score + cn_score) / 2.0
        weights[(p1, p2)] = combined_score
        weights[(p2, p1)] = combined_score  # Non-directed
    
    return weights

def save_weighted_network(
    output_path: str, 
    weights: Dict[Tuple[str, str], float],
    original_ppi: str
):
    """Sauvegarde le réseau avec combined_weight en remplacement de mean_weight."""
    ppi = pd.read_csv(original_ppi, sep='\t')
    
    # Remplacer mean_weight par le combined_weight
    ppi['weight'] = ppi.apply(
        lambda row: weights.get((row['protein1'], row['protein2']), 0.0),
        axis=1
    )
    
    # Sauvegarder avec les mêmes colonnes que l'entrée
    ppi.to_csv(output_path, sep='\t', index=False)

def main():
    # Chemins des fichiers
    ppi_file = "/Users/ryham/Documents/mémoireM2/Master_final_project/Data/clean data/weighted_networks/weighted_BIOGRID_levure.txt"  # Format: protein1 protein2 mean_weight
    go_file = "/Users/ryham/Documents/mémoireM2/Master_final_project/Data/go_slim_mapping.tab"          # Format SGD
    output_file = "/Users/ryham/Documents/mémoireM2/Master_final_project/Data/clean data/weighted_networks/tmp/GO_weighted_BIOGRID_levure.txt"
    
    # 1. Chargement du réseau PPI
    proteins, neighbors = load_ppi_network(ppi_file)
    print(f"Protéines dans le réseau: {len(proteins)}")
    
    # 2. Chargement des annotations GO
    go_annotations = load_go_annotations(go_file, proteins)
    print(f"Protéines avec annotations GO: {len(go_annotations)}")
    
    # 3. Calcul des poids combinés
    weights = calculate_combined_weights(proteins, neighbors, go_annotations)
    print(f"Paires pondérées: {len(weights)//2}")
    
    # 4. Sauvegarde (écrase mean_weight par combined_weight)
    save_weighted_network(output_file, weights, ppi_file)
    print(f"Réseau sauvegardé dans {output_file}")

if __name__ == "__main__":
    main()