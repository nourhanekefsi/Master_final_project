import pandas as pd
import re
from unidecode import unidecode

# Chemins des fichiers
complex_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\molecularComplexes.tab"
mapping_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\autres\YEAST_559292_idmapping.dat"
manual_mapping_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\autres\manual_mapping.tsv"
output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_SGD.txt"

# 1. Charger les mappings
def load_mappings():
    # Mapping automatique
    auto_map = {}
    with open(mapping_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[1] == "Gene_Name":
                clean_name = unidecode(parts[2].upper())
                auto_map[clean_name] = parts[0]
    
    # Mapping manuel
    manual_map = pd.read_csv(manual_mapping_file, sep='\t')
    manual_map = manual_map[manual_map['UniProtID'] != '-']
    for _, row in manual_map.iterrows():
        clean_name = unidecode(row['Original_Name'].upper())
        auto_map[clean_name] = row['UniProtID']
    
    return auto_map

# 2. Nettoyage des noms
def clean_name(raw_name):
    name = re.sub(r'^\d+X', '', raw_name, flags=re.IGNORECASE)
    name = re.sub(r'(DIMER|TRIMER|TETRAMER|OLIGOMER|COMPLEX)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'[^A-Z0-9,-]', '', name)
    return unidecode(name.strip().upper())

# 3. Traitement principal
def process_complexes():
    mapping = load_mappings()
    results = []
    complex_counter = 1  # Compteur pour les IDs numériques
    
    with open(complex_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                original_id = parts[0]
                subunits = parts[3]
                protein_list = []
                
                for subunit in subunits.split(':'):
                    cleaned = clean_name(subunit)
                    
                    # Gestion des multi-protéines
                    if ',' in cleaned:
                        for prot in cleaned.split(','):
                            prot = prot.strip()
                            if prot in mapping:
                                protein_list.append(mapping[prot])
                    elif cleaned in mapping:
                        protein_list.append(mapping[cleaned])
                
                if len(protein_list)>2:
                    results.append({
                        'id': complex_counter,  # ID numérique séquentiel
                        'original_id': original_id,  # ID original du complexe
                        'proteins': ';'.join(sorted(protein_list))
                    })
                    complex_counter += 1
    
    # Sauvegarde
    result_df = pd.DataFrame(results)
    result_df = result_df[['id', 'proteins']] 
    
    result_df['proteins'] = result_df['proteins'].drop_duplicates()
    
    result_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Fichier généré : {output_file}")
    print(f"Nombre total de complexes : {len(results)}")
    print(f"Nombre total de protéines uniques : {len(set(result_df['proteins'].str.cat(sep=';').split(';')))}")

if __name__ == "__main__":
    process_complexes()