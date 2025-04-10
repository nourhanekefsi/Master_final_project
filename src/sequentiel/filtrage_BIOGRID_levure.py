import pandas as pd
import re
from pathlib import Path

def extract_uniprot(alt_ids):
    """Extract UniProt IDs from alternative IDs field"""
    if pd.isna(alt_ids): 
        return None
    alt_ids = str(alt_ids)
    patterns = [
        r'uniprot/swiss[\W-]?prot:([A-Z0-9]{6,10})',
        r'uniprot:([A-Z0-9]{6,10})',
        r'swiss[\W-]?prot:([A-Z0-9]{6,10})',
        r'([A-Z0-9]{6,10})\.\d'
    ]
    for pattern in patterns:
        match = re.search(pattern, alt_ids, re.IGNORECASE)
        if match:
            return match.group(1).upper()
    return None

def process_biogrid_high_confidence():
    """Process BioGRID data to get ~50,000 high-confidence interactions"""
    # File paths
    input_file = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Protein_Interactions\BIOGRID-ORGANISM-Saccharomyces_cerevisiae.txt")
    output_file = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\BIOGRID_levure.txt")
    
    # Load data
    try:
        df = pd.read_csv(input_file, sep='\t', comment='#', header=None, dtype=str)
        df.columns = [
            'ID_A', 'ID_B', 'Alt_IDs_A', 'Alt_IDs_B', 
            'Aliases_A', 'Aliases_B', 'Method', 'Author',
            'PubIDs', 'TaxID_A', 'TaxID_B', 'IntType',
            'SourceDB', 'IntIDs', 'Confidence'
        ]
    except Exception as e:
        print(f"Error loading file: {e}")
        return
    
    # 1. Filter for Saccharomyces cerevisiae
    yeast_df = df[(df['TaxID_A'] == 'taxid:559292') & (df['TaxID_B'] == 'taxid:559292')].copy()
    
    # 2. Extract UniProt IDs
    yeast_df['Protein1'] = yeast_df['Alt_IDs_A'].apply(extract_uniprot)
    yeast_df['Protein2'] = yeast_df['Alt_IDs_B'].apply(extract_uniprot)
    
    # 3. Clean data
    clean_df = yeast_df.dropna(subset=['Protein1', 'Protein2'])
    clean_df = clean_df[clean_df['Protein1'] != clean_df['Protein2']].copy()
    
    # 4. High-confidence scoring system (less strict than ultra-strict version)
    method_scores = {
        'x-ray crystallography': 10,
        'electron microscopy': 8,
        'affinity chromatography': 6,
        'coimmunoprecipitation': 4,
        'two hybrid': 2,
        'pull down': 3,
        'mass spectrometry': 2
    }
    
    type_scores = {
        'direct interaction': 6,
        'physical association': 3,
        'complex': 4
    }
    
    clean_df['Pub_Count'] = clean_df['PubIDs'].str.count(r'\|').add(1).fillna(1)
    
    clean_df['Method_Score'] = clean_df['Method'].apply(
        lambda x: max([method_scores.get(method.lower(), 1)
                      for method in str(x).split('|') 
                      if 'psi-mi' not in method.lower()], default=1)
    )
    
    clean_df['Type_Score'] = clean_df['IntType'].apply(
        lambda x: max([type_scores.get(typ.lower(), 1)
                      for typ in str(x).split('|') 
                      if 'psi-mi' not in typ.lower()], default=1)
    )
    
    clean_df['Total_Score'] = (
        clean_df['Method_Score'] * 3 +
        clean_df['Type_Score'] * 2 +
        clean_df['Pub_Count']
    )
    
    # 5. Remove duplicates keeping highest scores
    clean_df['SortedPair'] = clean_df.apply(
        lambda x: tuple(sorted([x['Protein1'], x['Protein2']])), axis=1
    )
    
    non_redundant = clean_df.sort_values('Total_Score', ascending=False).drop_duplicates(subset=['SortedPair'])
    
    # 6. Dynamic threshold to get ~50,000 interactions
    target_count = 50000
    if len(non_redundant) > target_count:
        # Find score threshold that gives us ~50,000 interactions
        thresholds = sorted(non_redundant['Total_Score'].unique(), reverse=True)
        for threshold in thresholds:
            filtered = non_redundant[non_redundant['Total_Score'] >= threshold]
            if len(filtered) <= target_count:
                break
        
        # If we're still too far from target, take top N
        if len(filtered) < target_count * 0.8 or len(filtered) > target_count * 1.2:
            filtered = non_redundant.head(target_count)
    else:
        filtered = non_redundant
        print(f"Warning: Only {len(non_redundant)} interactions available")
    
    # 7. Statistics
    unique_proteins = pd.unique(filtered[['Protein1', 'Protein2']].values.ravel('K'))
    
    print(f"\nRésultats finaux (haute confiance):")
    print(f"- Interactions levures totales : {len(yeast_df):,}")
    print(f"- Interactions avec UniProt valides : {len(clean_df):,}")
    print(f"- Interactions non-redondantes : {len(non_redundant):,}")
    print(f"- Interactions sélectionnées : {len(filtered):,}")
    print(f"- Protéines uniques : {len(unique_proteins):,}")
    print(f"- Score minimum retenu : {filtered['Total_Score'].min():.1f}")
    print(f"- Score maximum : {filtered['Total_Score'].max():.1f}")
    print(f"- Score moyen : {filtered['Total_Score'].mean():.1f}")
    
    # 8. Save
    try:
        filtered[['Protein1', 'Protein2']].to_csv(
            output_file, 
            sep='\t', 
            index=False, 
            header=False
        )
        print(f"\nDonnées sauvegardées dans : {output_file}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde : {e}")

if __name__ == "__main__":
    process_biogrid_high_confidence()