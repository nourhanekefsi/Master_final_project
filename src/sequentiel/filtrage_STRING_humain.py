import pandas as pd
from pathlib import Path

# Chemins des fichiers
input_path = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Protein_Interactions\9606.protein.links.detailed.v12.0.txt")
output_path = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\STRING_humain_filtered_interactions.txt")

# 1. Chargement et filtrage initial
df = pd.read_csv(input_path, sep=" ")

# Créer un masque de filtrage
mask = (df["combined_score"] >= 700) & (
    df[["experimental", "coexpression", "database", "textmining"]].max(axis=1) > 90
)

# Appliquer le filtre et faire une copie explicite
filtered_df = df.loc[mask].copy()

# 2. Nettoyage des données
# Supprimer "4932." des noms de protéines
for col in ['protein1', 'protein2']:
    filtered_df.loc[:, col] = filtered_df[col].str.replace('9606.', '')

# 3. Gestion des interactions uniques
# Créer des identifiants d'interaction canoniques
filtered_df.loc[:, 'interaction_key'] = filtered_df.apply(
    lambda x: frozenset({x['protein1'], x['protein2']}), axis=1
)

# Supprimer les doublons en gardant la première occurrence
unique_interactions_df = filtered_df.drop_duplicates(subset='interaction_key')

# 4. Sauvegarde des résultats
unique_interactions_df[['protein1', 'protein2']].to_csv(
    output_path, sep="\t", index=False, header=False
)

# 5. Calcul des statistiques
unique_proteins = pd.unique(
    unique_interactions_df[['protein1', 'protein2']].values.ravel('K')
)
num_interactions = len(unique_interactions_df)

# Affichage des résultats
print("\n STRING ----- humain")
print("Traitement terminé avec succès.")
print(f"Nombre de protéines uniques : {len(unique_proteins)}")
print(f"Nombre d'interactions uniques : {num_interactions}")