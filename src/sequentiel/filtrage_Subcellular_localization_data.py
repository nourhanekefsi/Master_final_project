from typing import Dict, List, Set
import pandas as pd
from pathlib import Path
from collections import defaultdict
import numpy as np
from tqdm import tqdm

# Configuration des chemins
BASE_DIR = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\Github_CODE\Data")
CLEAN_DATA_DIR = BASE_DIR / "clean data"
RAW_DATA_DIR = BASE_DIR / "raw data" / "autres"
OUTPUT_DIR = BASE_DIR / "clean data" / "autres"
INTERACTIONS_DIR = CLEAN_DATA_DIR / "interactions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Fichiers d'entrée
HUMAN_COMPARTMENT = RAW_DATA_DIR / "human_compartment_integrated_full.tsv"
YEAST_COMPARTMENT = RAW_DATA_DIR / "yeast_compartment_integrated_full.tsv"

def load_compartment_data(compartment_file: Path, threshold: float = 1.0) -> Dict[str, Set[str]]:
    """
    Charge les données de localisation subcellulaire et filtre par score de confiance.
    
    Args:
        compartment_file: Chemin vers le fichier de localisation
        threshold: Seuil de score minimum pour considérer une localisation
        
    Returns:
        Dictionnaire {protéine: set(localisations)}
    """
    try:
        # Chargement des données
        df = pd.read_csv(compartment_file, sep='\t', header=None,
                        names=['protein_id', 'gene_name', 'go_id', 'go_term', 'score'])
        
        # Filtrage par score et nettoyage
        df = df[df['score'] >= threshold]
        df['protein_id'] = df['protein_id'].str.upper().str.strip()
        df['go_term'] = df['go_term'].str.strip()
        
        # Groupement par protéine
        protein_locations = defaultdict(set)
        for _, row in df.iterrows():
            protein_locations[row['protein_id']].add(row['go_term'])
            # On garde aussi le mapping par nom de gène si différent
            if row['gene_name'] != row['protein_id']:
                protein_locations[row['gene_name'].upper()].add(row['go_term'])
        
        return dict(protein_locations)
    except Exception as e:
        print(f"Erreur lors du chargement du fichier {compartment_file}: {e}")
        return {}

def calculate_sl_similarity(protein1: str, protein2: str, compartment_data: Dict[str, Set[str]]) -> float:
    """
    Calcule la similarité de localisation subcellulaire selon l'équation 5.
    
    Args:
        protein1: Première protéine
        protein2: Seconde protéine
        compartment_data: Dictionnaire des localisations
        
    Returns:
        Score de similarité SL entre 0 et 1
    """
    loc1 = compartment_data.get(protein1, set())
    loc2 = compartment_data.get(protein2, set())
    
    if not loc1 or not loc2:
        return 0.0
    
    intersection = len(loc1 & loc2)
    return (intersection ** 2) / (len(loc1) * len(loc2))

def create_sl_matrix(proteins: List[str], compartment_data: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    Crée une matrice de similarité SL pour toutes les paires de protéines.
    
    Args:
        proteins: Liste des protéines à inclure
        compartment_data: Dictionnaire des localisations
        
    Returns:
        DataFrame de similarité (matrice carrée)
    """
    n = len(proteins)
    sl_matrix = np.zeros((n, n))
    
    # Création d'un index pour un accès rapide
    protein_index = {prot: i for i, prot in enumerate(proteins)}
    
    # Calcul pour toutes les paires uniques
    for i in tqdm(range(n), desc="Calcul SL"):
        for j in range(i, n):
            similarity = calculate_sl_similarity(
                proteins[i], proteins[j], compartment_data)
            sl_matrix[i][j] = similarity
            sl_matrix[j][i] = similarity  # Matrice symétrique
    
    return pd.DataFrame(sl_matrix, index=proteins, columns=proteins)

def main():
    # Chargement des données de localisation
    print("Chargement des données de localisation subcellulaire...")
    human_comp = load_compartment_data(HUMAN_COMPARTMENT)
    yeast_comp = load_compartment_data(YEAST_COMPARTMENT)
    
    # Exemple d'utilisation avec un réseau PPI
    # (À adapter selon vos besoins spécifiques)
    networks = {
        "human": {
            "compartment_data": human_comp,
            "ppi_files": ["BIOGRID_humain", "STRING_humain"]
        },
        "yeast": {
            "compartment_data": yeast_comp,
            "ppi_files": ["BIOGRID_levure", "DIP_levure", "STRING_levure"]
        }
    }
    
    for organism, data in networks.items():
        print(f"\nTraitement des réseaux {organism}...")
        compartment_data = data["compartment_data"]
        
        for network in data["ppi_files"]:
            print(f"  Traitement du réseau {network}...")
            
            # Charger le PPI (à adapter selon votre structure)
            ppi_file = INTERACTIONS_DIR / f"{network}.txt"
            if not ppi_file.exists():
                print(f"Fichier PPI {ppi_file} non trouvé")
                continue
                
            # Chargement des interactions
            ppi = pd.read_csv(ppi_file, sep='\t', names=['protein1', 'protein2'])
            proteins = sorted(set(ppi['protein1']).union(set(ppi['protein2'])))
            
            # Création de la matrice SL
            sl_matrix = create_sl_matrix(proteins, compartment_data)
            
            # Sauvegarde
            output_file = OUTPUT_DIR / f"SL_matrix_{network}.tsv"
            sl_matrix.to_csv(output_file, sep='\t', float_format='%.5f')
            print(f"Matrice SL sauvegardée dans {output_file}")

if __name__ == "__main__":
    main()