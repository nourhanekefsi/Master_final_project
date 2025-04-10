import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import os
import glob

def load_ppi_network(file_path):
    """Charge un réseau PPI à partir d'un fichier"""
    ppi = pd.read_csv(file_path, sep='\t', header=None, names=['protein1', 'protein2'])
    return ppi

def build_graph(ppi_df):
    """Construit un graphe NetworkX à partir d'un dataframe PPI"""
    G = nx.Graph()
    for _, row in ppi_df.iterrows():
        G.add_edge(row['protein1'], row['protein2'])
    return G

def calculate_hcn_similarity(G, output_file):
    """Calcule la similarité HCN pour toutes les paires de protéines connectées"""
    results = []
    
    # Pré-calculer les voisins pour tous les nœuds
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes()}
    
    for edge in G.edges():
        v, u = edge
        N_v = neighbors[v]
        N_u = neighbors[u]
        
        # Calculer les composantes de la formule HCN
        NCN = N_v & N_u  # Intersection des voisins (common neighbors)
        N_union = N_v | N_u  # Union des voisins
        
        # Calculer chaque terme de la formule
        numerator = len(NCN)**2
        denominator = (len(N_v) * len(N_u) * len(N_union))
        
        # Éviter la division par zéro
        if denominator == 0:
            hcn = 0.0
        else:
            hcn = numerator / denominator
        
        results.append({
            'protein1': v,
            'protein2': u,
            'HCN_score': hcn,
            'common_neighbors': len(NCN),
            'degree_v': len(N_v),
            'degree_u': len(N_u)
        })
    
    # Sauvegarder les résultats
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"Résultats HCN sauvegardés dans {output_file}")
    return result_df

def process_all_ppi_networks(base_dir, output_dir):
    """Traite tous les fichiers PPI dans le dossier spécifié"""
    # Créer le dossier de sortie s'il n'existe pas
    os.makedirs(output_dir, exist_ok=True)
    
    # Trouver tous les fichiers PPI
    ppi_files = glob.glob(os.path.join(base_dir, "*i*"))[:5]  # Prendre les 5 premiers fichiers
    
    for ppi_file in ppi_files:
        print(f"\nTraitement du fichier: {os.path.basename(ppi_file)}")
        
        # Charger le réseau PPI
        ppi_df = load_ppi_network(ppi_file)
        print(f"Nombre d'interactions chargées: {len(ppi_df)}")
        
        # Construire le graphe
        G = build_graph(ppi_df)
        print(f"Nombre de nœuds: {G.number_of_nodes()}")
        print(f"Nombre d'arêtes: {G.number_of_edges()}")
        
        # Calculer les scores HCN
        output_name = os.path.splitext(os.path.basename(ppi_file))[0]
        output_file = os.path.join(output_dir, f"HCN_scores_{output_name}.txt")
        
        hcn_results = calculate_hcn_similarity(G, output_file)
        
        # Afficher quelques statistiques
        print(f"Score HCN moyen: {hcn_results['HCN_score'].mean():.4f}")
        print(f"Nombre moyen de voisins communs: {hcn_results['common_neighbors'].mean():.2f}")

if __name__ == "__main__":
    # Configuration des chemins
    base_dir = "C:/Users/PC/Documents/M2 HPC/PFE/PFE_CODE/Data/clean data/interactions"
    output_dir = "C:/Users/PC/Documents/M2 HPC/PFE/PFE_CODE/Data/clean data/autres"
    
    # Exécuter le traitement pour tous les réseaux PPI
    process_all_ppi_networks(base_dir, output_dir)