import pandas as pd
import numpy as np
import os
from pathlib import Path
from collections import defaultdict

# Configuration des chemins
BASE_DIR = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\Github_CODE\Data")
CLEAN_DATA_DIR = BASE_DIR / "clean data"
INTERACTIONS_DIR = CLEAN_DATA_DIR / "interactions"
OUTPUT_DIR = CLEAN_DATA_DIR / "weighted_ppi"
SIMILARITY_DIR = CLEAN_DATA_DIR / "autres"

# Créer le dossier de sortie
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_similarity_scores(network_name):
    """Charge tous les scores de similarité disponibles pour un réseau"""
    scores = defaultdict(dict)
    
    # 1. Charger les scores HCN
    hcn_file = SIMILARITY_DIR / f"HCN_scores_{network_name}.txt"
    if hcn_file.exists():
        hcn_df = pd.read_csv(hcn_file, sep='\t')
        for _, row in hcn_df.iterrows():
            pair = tuple(sorted([row['protein1'], row['protein2']]))
            scores[pair]['HCN'] = row['HCN_score']
    
    # 2. Charger les scores de co-expression
    coexpr_file = SIMILARITY_DIR / f"coexpression_{network_name}.txt"
    if coexpr_file.exists():
        coexpr_df = pd.read_csv(coexpr_file, sep='\t')
        for _, row in coexpr_df.iterrows():
            pair = tuple(sorted([row['protein1'], row['protein2']]))
            scores[pair]['PCC'] = row['PCC']
    
    # 3. Charger les scores fonctionnels
    func_file = SIMILARITY_DIR / f"matrice_S_{network_name}.txt"
    if func_file.exists():
        func_df = pd.read_csv(func_file, sep='\t', index_col=0)
        proteins = func_df.index
        for i in range(len(proteins)):
            for j in range(i+1, len(proteins)):
                pair = tuple(sorted([proteins[i], proteins[j]]))
                scores[pair]['Func'] = func_df.iloc[i, j]
    
    # 4. Charger les scores de localisation
    sl_file = SIMILARITY_DIR / f"SL_matrix_{network_name}.tsv"
    if sl_file.exists():
        sl_df = pd.read_csv(sl_file, sep='\t', index_col=0)
        proteins = sl_df.index
        for i in range(len(proteins)):
            for j in range(i+1, len(proteins)):
                pair = tuple(sorted([proteins[i], proteins[j]]))
                scores[pair]['SL'] = sl_df.iloc[i, j]
    
    return scores

def calculate_weighted_ppi(network_name):
    """Calcule la moyenne pondérée pour un réseau PPI"""
    print(f"\nTraitement du réseau {network_name}")
    
    # Charger les scores de similarité
    similarity_scores = load_similarity_scores(network_name)
    if not similarity_scores:
        print("Aucun score de similarité trouvé")
        return None
    
    # Calculer la moyenne pondérée
    weighted_interactions = []
    for pair, scores in similarity_scores.items():
        # Calculer la moyenne des scores disponibles
        valid_scores = [s for s in scores.values() if not np.isnan(s)]
        if valid_scores:
            weight = np.mean(valid_scores)
            weighted_interactions.append({
                'protein1': pair[0],
                'protein2': pair[1],
                'weight': weight
            })
    
    if not weighted_interactions:
        print("Aucune interaction valide avec scores")
        return None
    
    # Créer le DataFrame final
    result_df = pd.DataFrame(weighted_interactions)
    return result_df[['protein1', 'protein2', 'weight']]

def process_all_networks():
    """Traite tous les réseaux PPI disponibles"""
    # Lister tous les fichiers PPI originaux
    ppi_files = list(INTERACTIONS_DIR.glob("*.txt"))
    
    for ppi_file in ppi_files:
        network_name = ppi_file.stem
        weighted_ppi = calculate_weighted_ppi(network_name)
        
        if weighted_ppi is not None:
            output_file = OUTPUT_DIR / f"weighted_{network_name}.txt"
            weighted_ppi.to_csv(output_file, sep='\t', index=False, header=False)
            print(f"Fichier sauvegardé: {output_file} ({len(weighted_ppi)} interactions)")

if __name__ == "__main__":
    process_all_networks()