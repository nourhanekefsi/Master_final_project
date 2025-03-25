import pandas as pd
import re

# Charger le fichier
input_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Protein_Interactions\BIOGRID-MV-Physical.txt"
output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\BIOGRID_filtered_interactions.txt"

# Charger les données
df = pd.read_csv(input_file, sep='\t', low_memory=False)

# Renommer la colonne pour éviter le problème du "#"
df.rename(columns={"#ID Interactor A": "ID Interactor A"}, inplace=True)

# Colonnes importantes
col_alt_a = "Alt IDs Interactor A"
col_alt_b = "Alt IDs Interactor B"
col_taxid_a = "Taxid Interactor A"
col_taxid_b = "Taxid Interactor B"
col_interaction_type = "Interaction Types"

# 1. Filtrer pour garder uniquement Homo sapiens (taxid:9606)
df = df[(df[col_taxid_a] == "taxid:9606") & (df[col_taxid_b] == "taxid:9606")]

# 2. Sélectionner uniquement les interactions protéine-protéine
df = df[df[col_interaction_type].str.contains("MI:0407", na=False)]  # "MI:0407" = direct interaction

# 3. Extraire les identifiants UniProt
def extract_uniprot(alt_ids):
    """Extrait l'identifiant UniProt depuis la colonne 'Alt IDs'."""
    if isinstance(alt_ids, str):
        match = re.search(r"uniprot/swiss-prot:([\w-]+)", alt_ids)
        if match:
            return match.group(1)  # Retourne uniquement l'ID UniProt
    return None  # Si pas d'UniProt, retourner None

df["Proteine1"] = df[col_alt_a].apply(extract_uniprot)
df["Proteine2"] = df[col_alt_b].apply(extract_uniprot)

# 4. Supprimer les lignes où l'un des deux identifiants est manquant
df_filtered = df.dropna(subset=["Proteine1", "Proteine2"])

# 5. Garder uniquement les colonnes nécessaires et sauvegarder
df_final = df_filtered[["Proteine1", "Proteine2"]]
df_final.to_csv(output_file, sep='\t', index=False, header=False)

print(f"Fichier filtré enregistré sous : {output_file}")
