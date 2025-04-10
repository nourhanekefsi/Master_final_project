import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import os
from collections import defaultdict
import glob

def load_mapping(file_path):
    """Charge les mappings d'identifiants avec tous les types pertinents"""
    mapping = {
        'uniprot_to_gene': defaultdict(list),
        'uniprot_to_name': defaultdict(list),
        'gene_to_uniprot': defaultdict(list),
        'name_to_uniprot': defaultdict(list),
        'orf_to_uniprot': defaultdict(list),
        'uniprot_to_uniprot': defaultdict(list),
        'ensembl_to_uniprot': defaultdict(list),
        'string_to_uniprot': defaultdict(list)
    }
    
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                uniprot_id = parts[0]
                db = parts[1]
                db_id = parts[2]
                
                # Mapper Uniprot à lui-même
                mapping['uniprot_to_uniprot'][uniprot_id].append(uniprot_id)
                
                if db == 'Gene_OrderedLocusName':
                    mapping['uniprot_to_gene'][uniprot_id].append(db_id)
                    mapping['gene_to_uniprot'][db_id].append(uniprot_id)
                elif db == 'Gene_Name':
                    mapping['uniprot_to_name'][uniprot_id].append(db_id)
                    mapping['name_to_uniprot'][db_id.upper()].append(uniprot_id)
                elif db == 'ORF':
                    mapping['orf_to_uniprot'][db_id].append(uniprot_id)
                elif db == 'Ensembl':
                    mapping['ensembl_to_uniprot'][db_id].append(uniprot_id)
                elif db == 'STRING':
                    mapping['string_to_uniprot'][db_id].append(uniprot_id)
    
    # Simplifier les mappings en gardant le premier élément
    simple_mappings = {}
    for map_type in mapping:
        simple_mappings[map_type] = {k: v[0] for k, v in mapping[map_type].items() if v}
    
    return simple_mappings

def load_expression_data(file_path, mappings, species):
    """Charge les données d'expression en fonction de l'espèce"""
    try:
        # Trouver le début des données
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('!dataset_table_begin'):
                    break
            
            # Lire les données avec toutes les colonnes
            df = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
        
        # Extraire les colonnes d'expression (commençant par GSM)
        expression_cols = [col for col in df.columns if col.startswith('GSM')]
        
        # Vérifier que toutes les colonnes GSM ont le même nombre de valeurs
        gsm_lengths = {col: len(df[col].dropna()) for col in expression_cols}
        if len(set(gsm_lengths.values())) > 1:
            print("Avertissement: Les colonnes GSM ont des longueurs différentes!")
            # Prendre le nombre minimal de valeurs
            min_length = min(gsm_lengths.values())
            for col in expression_cols:
                df[col] = df[col].dropna().head(min_length)
        
        # Garder seulement les colonnes nécessaires
        df = df[['IDENTIFIER'] + expression_cols].dropna()
        df = df.set_index('IDENTIFIER')
        
        # Convertir en numérique en gérant les erreurs
        for col in expression_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df = df.dropna(how='all')
        
        # Mapping différent selon l'espèce
        if species == 'human':
            # Pour les données humaines, essayer plusieurs stratégies
            def human_mapper(x):
                try:
                    x = str(x).upper()
                    if x in mappings['name_to_uniprot']:
                        return mappings['name_to_uniprot'][x]
                    return None
                except:
                    return None
            
            df.index = df.index.map(human_mapper)
        else:
            # Pour la levure, utiliser les ORFs/noms de gènes
            def yeast_mapper(x):
                try:
                    x = str(x)
                    if x in mappings['gene_to_uniprot']:
                        return mappings['gene_to_uniprot'][x]
                    if x in mappings['orf_to_uniprot']:
                        return mappings['orf_to_uniprot'][x]
                    if x.upper() in mappings['name_to_uniprot']:
                        return mappings['name_to_uniprot'][x.upper()]
                    return None
                except:
                    return None
            
            df.index = df.index.map(yeast_mapper)
        
        return df.dropna()
    except Exception as e:
        print(f"Erreur lors du chargement des données d'expression: {str(e)}")
        return pd.DataFrame()

def calculate_pcc(v, u):
    """Calcule le coefficient de corrélation de Pearson avec vérification de longueur"""
    try:
        # Convertir en arrays numpy et s'assurer qu'ils sont 1D
        v = np.asarray(v, dtype=float).flatten()
        u = np.asarray(u, dtype=float).flatten()
        
        # Trouver la longueur minimale
        min_len = min(len(v), len(u))
        if min_len < 2:
            return np.nan
            
        # Tronquer les vecteurs à la même longueur
        v = v[:min_len]
        u = u[:min_len]
        
        # Masque pour les valeurs manquantes
        mask = ~np.isnan(v) & ~np.isnan(u)
        v_clean = v[mask]
        u_clean = u[mask]
        
        if len(v_clean) < 2:
            return np.nan
        
        # Vérifier si l'un des vecteurs est constant
        if np.all(v_clean == v_clean[0]) or np.all(u_clean == u_clean[0]):
            return np.nan
            
        pcc, _ = pearsonr(v_clean, u_clean)
        return (pcc + 1) / 2  # Transformation vers [0,1]
    except Exception as e:
        print(f"Error calculating PCC between shapes {np.shape(v)} and {np.shape(u)}: {str(e)}")
        return np.nan

def calculate_coexpression(ppi, expr_data, output_file):
    """Calcule et sauvegarde les scores de co-expression"""
    results = []
    proteins_in_ppi = set(ppi['protein1']).union(set(ppi['protein2']))
    proteins_in_expr = set(expr_data.index)
    common_proteins = proteins_in_ppi.intersection(proteins_in_expr)
    
    print(f"\nRecouvrement entre PPI et expression:")
    print(f"- Protéines dans PPI: {len(proteins_in_ppi)}")
    print(f"- Protéines dans expression: {len(proteins_in_expr)}")
    print(f"- Protéines communes: {len(common_proteins)}")
    
    if len(common_proteins) < 2:
        print("Avertissement: Pas assez de protéines communes pour calculer la co-expression!")
        return
    
    # Pré-traitement: s'assurer que toutes les protéines ont le même nombre de points de données
    min_length = min([len(expr_data.loc[p].values) for p in common_proteins])
    expr_data = expr_data.apply(lambda x: x[:min_length], axis=1)
    
    for _, row in ppi.iterrows():
        p1, p2 = row['protein1'], row['protein2']
        if p1 in expr_data.index and p2 in expr_data.index:
            # Get expression values and ensure they're 1D arrays
            v = expr_data.loc[p1].values
            u = expr_data.loc[p2].values
            
            # Flatten arrays in case they're 2D
            v = np.asarray(v).flatten()
            u = np.asarray(u).flatten()
            
            # Skip if one array is empty
            if len(v) == 0 or len(u) == 0:
                continue
                
            pcc = calculate_pcc(v, u)
            if not np.isnan(pcc):
                results.append({'protein1': p1, 'protein2': p2, 'PCC': pcc})
    
    if results:
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"\nRésultats sauvegardés dans {output_file} ({len(result_df)} paires valides)")
    else:
        print("\nAucun résultat valide à sauvegarder.")

def process_dataset(base_dir, species):
    """Traite un ensemble de données complet"""
    print(f"\n=== Traitement des données {species} ===")
    
    # Déterminer les noms de fichiers en fonction de l'espèce
    if species == 'human':
        mapping_file = os.path.join(base_dir, "raw data/autres/HUMAN_9606_idmapping.dat")
        expr_file = os.path.join(base_dir, "raw data/autres/human_co-expression.soft")
        ppi_pattern = os.path.join(base_dir, "clean data/interactions/*humain*")
    else:  # levure
        mapping_file = os.path.join(base_dir, "raw data/autres/YEAST_559292_idmapping.dat")
        expr_file = os.path.join(base_dir, "raw data/autres/levure_co-expression.soft")
        ppi_pattern = os.path.join(base_dir, "clean data/interactions/*levure*")
    
    # Charger les mappings
    if not os.path.exists(mapping_file):
        print(f"Fichier de mapping introuvable: {mapping_file}")
        return
    
    mappings = load_mapping(mapping_file)
    print(f"\nStatistiques de mapping:")
    print(f"- Uniprot IDs: {len(mappings['uniprot_to_uniprot'])}")
    print(f"- Gene names: {len(mappings['name_to_uniprot'])}")
    print(f"- Gene IDs: {len(mappings['gene_to_uniprot'])}")
    print(f"- ORFs: {len(mappings['orf_to_uniprot'])}")
    print(f"- Ensembl IDs: {len(mappings['ensembl_to_uniprot'])}")
    print(f"- STRING IDs: {len(mappings['string_to_uniprot'])}")
    
    # Charger les données d'expression
    if not os.path.exists(expr_file):
        print(f"Fichier d'expression introuvable: {expr_file}")
        return
    
    expr_data = load_expression_data(expr_file, mappings, species)
    if expr_data.empty:
        print("Avertissement: Aucune donnée d'expression valide après mapping!")
        return
    
    print(f"\nDonnées d'expression chargées: {len(expr_data)} protéines")
    
    # Traiter les fichiers PPI
    ppi_files = glob.glob(ppi_pattern)
    if not ppi_files:
        print(f"Aucun fichier PPI trouvé avec le pattern: {ppi_pattern}")
    
    for ppi_file in ppi_files:
        print(f"\nTraitement de {os.path.basename(ppi_file)}")
        try:
            ppi = pd.read_csv(ppi_file, sep='\t', header=None, names=['protein1', 'protein2'], dtype=str)
            
            # Mapper les identifiants
            ppi['protein1'] = ppi['protein1'].apply(lambda x: mappings['uniprot_to_uniprot'].get(x, None))
            ppi['protein2'] = ppi['protein2'].apply(lambda x: mappings['uniprot_to_uniprot'].get(x, None))
            ppi = ppi.dropna()
            
            if ppi.empty:
                print("Avertissement: Aucune interaction valide après mapping!")
                continue
            
            print(f"Interactions chargées: {len(ppi)}")
            
            # Calculer la co-expression
            output_name = os.path.splitext(os.path.basename(ppi_file))[0]
            output_file = os.path.join(base_dir, f"clean data/autres/coexpression_{output_name}.txt")
            
            calculate_coexpression(ppi, expr_data, output_file)
        except Exception as e:
            print(f"Erreur lors du traitement: {str(e)}")

def main():
    base_dir = "C:/Users/PC/Documents/M2 HPC/PFE/PFE_CODE/Data"
    
    # Créer le dossier de résultats s'il n'existe pas
    os.makedirs(os.path.join(base_dir, "clean data/autres"), exist_ok=True)
    
    # Traiter les deux espèces
    process_dataset(base_dir, 'human')
    process_dataset(base_dir, 'levure')

if __name__ == "__main__":
    main()