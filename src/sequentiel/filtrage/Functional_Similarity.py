import pandas as pd
import numpy as np
from pathlib import Path
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity
import re
from typing import Dict, Tuple, Optional
import logging

# Configuration des chemins
BASE_DIR = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data")
RAW_DATA_DIR = BASE_DIR / "raw data" / "autres"
CLEAN_DATA_DIR = BASE_DIR / "clean data"
INTERACTIONS_DIR = CLEAN_DATA_DIR / "interactions"
OUTPUT_DIR = CLEAN_DATA_DIR / "autres"
GO_SLIM_YEAST = RAW_DATA_DIR / "go_slim_mapping.tab"
GO_SLIM_HUMAN = RAW_DATA_DIR / "uniprotkb_Homo_sapiens_Human_AND_model_2025_04_03.tsv"
MAPPING_FILE = RAW_DATA_DIR / "YEAST_559292_idmapping.dat"

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_mapping(mapping_file: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Charge le mapping UniProt vers SGD et inversement"""
    try:
        dtype = {'UniProt': 'string', 'DB': 'string', 'ID': 'string'}
        df = pd.read_csv(mapping_file, sep='\t', header=None,
                        names=['UniProt', 'DB', 'ID'], dtype=dtype)
        
        sgd_mapping = df[df['DB'] == 'SGD']
        uniprot_to_sgd = dict(zip(sgd_mapping['UniProt'], sgd_mapping['ID']))
        sgd_to_uniprot = dict(zip(sgd_mapping['ID'], sgd_mapping['UniProt']))
        
        logger.info(f"Loaded {len(uniprot_to_sgd)} UniProt to SGD mappings")
        return uniprot_to_sgd, sgd_to_uniprot
    except Exception as e:
        logger.error(f"Erreur lors du chargement du mapping: {e}")
        return {}, {}

def load_go_slim(go_file: Path, is_human: bool = True, 
                sgd_to_uniprot: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Charge les annotations GO et convertit les identifiants en UniProt si nécessaire"""
    try:
        if is_human:
            # Chargement des données humaines (déjà en UniProt)
            df = pd.read_csv(go_file, sep='\t', usecols=['Entry', 'Gene Ontology (GO)'], 
                            dtype={'Entry': 'string'})
            
            # Extraction des termes GO
            go_data = []
            for entry, go_str in df.dropna(subset=['Gene Ontology (GO)']).itertuples(index=False):
                terms = re.findall(r'GO:\d+', go_str)
                go_data.extend([(entry, term) for term in terms])
            
            df_go = pd.DataFrame(go_data, columns=['protein', 'go_term'])
        else:
            # Chargement des données de levure (SGD -> besoin de conversion)
            df = pd.read_csv(go_file, sep='\t', header=None, 
                           usecols=[2, 5], names=['protein', 'go_term'],
                           dtype={'protein': 'string', 'go_term': 'string'})
            
            # Filtrage des termes GO valides
            df_go = df[df['go_term'].str.contains(r'^GO:\d+$', na=False)].copy()
            
            # Conversion SGD -> UniProt si le mapping est fourni
            if sgd_to_uniprot:
                logger.info("Converting SGD IDs to UniProt IDs in GO annotations")
                df_go['protein'] = df_go['protein'].map(sgd_to_uniprot)
                df_go = df_go.dropna(subset=['protein'])  # Enlever les non mappés
        
        logger.info(f"Loaded {len(df_go)} GO annotations")
        return df_go.drop_duplicates().reset_index(drop=True)
    except Exception as e:
        logger.error(f"Erreur lors du chargement des annotations GO: {e}")
        return pd.DataFrame(columns=['protein', 'go_term'])

def calculate_functional_similarity(proteins: list, go_annotations: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Calcule la similarité fonctionnelle entre les protéines"""
    try:
        proteins_set = sorted(set(proteins))
        logger.info(f"Calculating functional similarity for {len(proteins_set)} unique proteins")
        
        # Filtrer les annotations pour les protéines d'intérêt
        go_filtered = go_annotations[go_annotations['protein'].isin(proteins_set)].copy()
        
        if len(go_filtered) == 0:
            logger.warning("No matching GO annotations found for the proteins")
            return None
        
        # Créer une matrice protéine-terme GO
        protein_cat = pd.CategoricalDtype(categories=proteins_set, ordered=True)
        go_filtered['protein'] = go_filtered['protein'].astype(protein_cat)
        
        # Construction de la matrice creuse
        row_ind = go_filtered['protein'].cat.codes
        col_ind = pd.Categorical(go_filtered['go_term']).codes
        data = np.ones(len(go_filtered), dtype=np.int8)
        
        sparse_matrix = sparse.csr_matrix(
            (data, (row_ind, col_ind)),
            shape=(len(proteins_set), len(go_filtered['go_term'].unique()))
        )
        # Calcul de similarité cosinus
        logger.info("Computing cosine similarity...")
        similarity_matrix = cosine_similarity(sparse_matrix)
        
        return pd.DataFrame(similarity_matrix, 
                          index=proteins_set, 
                          columns=proteins_set)
    except Exception as e:
        logger.error(f"Error calculating functional similarity: {e}")
        return None

def load_ppi_network(ppi_file: Path) -> pd.DataFrame:
    """Charge un réseau PPI (supposé utiliser des UniProt IDs)"""
    try:
        ppi = pd.read_csv(ppi_file, sep='\t', header=None,
                         names=['protein1', 'protein2'], dtype='string')
        
        # Nettoyage des identifiants
        ppi['protein1'] = ppi['protein1'].str.upper().str.strip()
        ppi['protein2'] = ppi['protein2'].str.upper().str.strip()
        
        # Suppression des doublons et des valeurs manquantes
        ppi = ppi.dropna().drop_duplicates()
        
        logger.info(f"Loaded PPI network with {len(ppi)} interactions")
        return ppi
    except Exception as e:
        logger.error(f"Error loading PPI network: {e}")
        return pd.DataFrame(columns=['protein1', 'protein2'])

def main():
    # Chargement des mappings
    uniprot_to_sgd, sgd_to_uniprot = load_mapping(MAPPING_FILE)
    
    # Configuration des réseaux à traiter
    networks = {
        "STRING_levure": (GO_SLIM_YEAST, False),
        "BIOGRID_humain": (GO_SLIM_HUMAN, True),
        "STRING_humain": (GO_SLIM_HUMAN, True),
        "BIOGRID_levure": (GO_SLIM_YEAST, False),
        "DIP_levure": (GO_SLIM_YEAST, False)   
    }

    for network_name, (go_file, is_human) in networks.items():
        logger.info(f"\nProcessing {network_name} network")
        
        # Chargement du réseau PPI
        ppi_file = INTERACTIONS_DIR / f"{network_name}.txt"
        ppi = load_ppi_network(ppi_file)
        
        if ppi.empty:
            logger.warning(f"No valid interactions found for {network_name}")
            continue
        
        # Chargement des annotations GO avec conversion si nécessaire
        go_annotations = load_go_slim(
            go_file, 
            is_human=is_human,
            sgd_to_uniprot=sgd_to_uniprot if not is_human else None
        )
        
        if go_annotations.empty:
            logger.warning(f"No GO annotations found for {network_name}")
            continue
        
        # Liste de toutes les protéines du réseau
        all_proteins = list(set(ppi['protein1']).union(set(ppi['protein2'])))
        
        # Calcul de la similarité fonctionnelle
        similarity_matrix = calculate_functional_similarity(all_proteins, go_annotations)
        
        if similarity_matrix is None:
            logger.warning(f"Could not compute similarity matrix for {network_name}")
            continue
        
        # Ajout des poids aux interactions
        ppi['weight'] = ppi.apply(
            lambda row: similarity_matrix.loc[row['protein1'], row['protein2']]
            if row['protein1'] in similarity_matrix.index and row['protein2'] in similarity_matrix.columns
            else 0.0,
            axis=1
        )
        
        # Filtrage des interactions avec poids > 0
        ppi = ppi[ppi['weight'] > 0]
        
        # Sauvegarde
        output_file = OUTPUT_DIR / f"FS_{network_name}.txt"
        ppi.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved weighted network to {output_file} with {len(ppi)} interactions")

if __name__ == "__main__":
    main()