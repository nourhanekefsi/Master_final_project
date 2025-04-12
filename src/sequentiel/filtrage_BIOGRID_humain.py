import pandas as pd
import re
from pathlib import Path

def process_biogrid_max_coverage():
    # Chemins des fichiers
    input_file = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Protein_Interactions\BIOGRID-MV-Physical.txt")
    output_file = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\BIOGRID_humain.txt")
    
    # Charger les données en forçant le type string
    df = pd.read_csv(input_file, sep='\t', comment='#', header=None, dtype=str)
    df.columns = [
        'ID_A', 'ID_B', 'Alt_IDs_A', 'Alt_IDs_B', 
        'Aliases_A', 'Aliases_B', 'Method', 'Author',
        'PubIDs', 'TaxID_A', 'TaxID_B', 'IntType',
        'SourceDB', 'IntIDs', 'Confidence'
    ]
    
    # 1. Filtrer pour Homo sapiens uniquement
    human_df = df[(df['TaxID_A'] == 'taxid:9606') & (df['TaxID_B'] == 'taxid:9606')].copy()
    
    # 2. Extraction robuste des UniProt IDs
    def extract_uniprot(alt_ids):
        if pd.isna(alt_ids): 
            return None
        patterns = [
            r'uniprot/swiss[\W-]?prot:([A-Z0-9]{6,8})',
            r'uniprot:([A-Z0-9]{6,8})',
            r'([A-Z0-9]{6,8})\.\d'
        ]
        for pattern in patterns:
            match = re.search(pattern, str(alt_ids), re.IGNORECASE)
            if match:
                return match.group(1)
        return None
    
    human_df['Protein1'] = human_df['Alt_IDs_A'].apply(extract_uniprot)
    human_df['Protein2'] = human_df['Alt_IDs_B'].apply(extract_uniprot)
    
    # 3. Nettoyage de base
    clean_df = human_df.dropna(subset=['Protein1', 'Protein2'])
    clean_df = clean_df[clean_df['Protein1'] != clean_df['Protein2']].copy()
    
    # 4. Système de scoring des interactions (conservé pour information)
    method_scores = {
        'x-ray crystallography': 4,
        'electron microscopy': 3,
        'affinity chromatography': 3,
        'coimmunoprecipitation': 2,
        'two hybrid': 2,
        'pull down': 2,
        'mass spectrometry': 1
    }
    
    type_scores = {
        'direct interaction': 3,
        'physical association': 1,
        'complex': 2
    }
    
    clean_df['Pub_Count'] = clean_df['PubIDs'].str.count(r'\|') + 1
    
    clean_df['Method_Score'] = clean_df['Method'].map(
        lambda x: max([method_scores.get(method.lower(), 1) 
                      for method in str(x).split('|')])
    )
    clean_df['Type_Score'] = clean_df['IntType'].map(
        lambda x: max([type_scores.get(typ.lower(), 1) 
                      for typ in str(x).split('|')])
    )
    clean_df['Total_Score'] = (
        clean_df['Method_Score'] * 3 + 
        clean_df['Type_Score'] * 2 + 
        clean_df['Pub_Count']
    )
    
    # 5. Supprimer les doublons (même paire de protéines dans un ordre différent)
    # Créer une colonne avec la paire triée pour identifier les doublons
    clean_df['SortedPair'] = clean_df.apply(lambda x: tuple(sorted([x['Protein1'], x['Protein2']])), axis=1)
    
    # Garder toutes les interactions uniques (sans filtrage par score)
    non_redundant = clean_df.drop_duplicates(subset=['SortedPair'])
    
    # 6. Statistiques finales
    unique_proteins = pd.unique(
        non_redundant[['Protein1', 'Protein2']].values.ravel('K')
    )
    
    print(f"\nRésultats finaux :")
    print(f"- Interactions humaines totales : {len(human_df)}")
    print(f"- Interactions avec UniProt valides : {len(clean_df)}")
    print(f"- Interactions non-redondantes : {len(non_redundant)}")
    print(f"- Protéines uniques : {len(unique_proteins)}")
    print(f"- Score moyen (pour information) : {non_redundant['Total_Score'].mean():.1f}")
    
    # 7. Sauvegarde (seulement les deux colonnes Protein1 et Protein2)
    non_redundant[['Protein1', 'Protein2']].to_csv(
        output_file, 
        sep='\t', 
        index=False, 
        header=False  # Pas d'en-tête dans le fichier de sortie
    )

if __name__ == "__main__":
    process_biogrid_max_coverage()