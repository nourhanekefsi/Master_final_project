import os
import panda as pd
from collections import defaultdict

# Chemins des fichiers d'entrée
input_files = [
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_SGD.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\Portal_complexes_levure.txt"
    
]

# Chemin du fichier de sortie
output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_levure.txt"

def process_complexes(input_files, output_file):
    unique_complexes = set()
    
    # Lire tous les fichiers et collecter les complexes uniques
    for file_path in input_files:
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    # Supprimer les numéros au début de la ligne et les espaces multiples
                    line = line.split()
                    if line:
                        # Supprimer le numéro au début si présent
                        parts = line[1].split(';')
                        proteins = [p for p in parts]
                        if len(proteins)>2:
                            complex_str = ';'.join(sorted(proteins))  # Tri pour l'uniformité
                            unique_complexes.add(complex_str)
        except FileNotFoundError:
            print(f"Attention: Fichier non trouvé - {file_path}")
            continue
    # Create a DataFrame with proper structure
    df_complexes = pd.DataFrame(
        [(idx, complex_str) for idx, complex_str in enumerate(sorted(unique_complexes), 1)],
        columns=['id', 'Complexes']
    )
    df_complexes = df_complexes.drop_duplicates(subset=['Complexes'])
    df_complexes.to_csv(output_file, sep='\t', index=False, header=False)
    # Écrire les résultats dans le fichier de sortie
    
    print(f"Terminé! {len(df_complexes)} complexes uniques ont été écrits dans {output_file}")

# Exécuter la fonction
process_complexes(input_files, output_file)