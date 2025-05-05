import os
import pandas as pd
import numpy as np

# Chemins des fichiers d'entrée
input_files = [
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_SGD.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\Portal_complexes_levure.txt"    
]

# Chemin du fichier de sortie
output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_levure.txt"

def process_complexes(input_files, output_file):
    unique_complexes = set([])
    unique_proteins = set()
    complexes = set([])
    # Lire tous les fichiers et collecter les complexes uniques
    for input_file in input_files:
        df = pd.read_csv(input_file, sep= '\t', header=0, names=[ 'id','proteins'])
        for complex in df['proteins']:
            complex = complex.split(';')
            complex = sorted(set(complex))
            complexes.add(tuple(complex))
            for p in complex:
                unique_proteins.add(p)
            complex = " ".join([c for c in complex])
            unique_complexes.add(complex)
            
    mean_length = []
    for complex in complexes:
        mean_length.append(len(complex))       
    mean_length = np.mean(mean_length)        
    df_complexes = pd.DataFrame(
        [(idx, str(complex_str)) for idx, complex_str in enumerate(sorted(unique_complexes), 1)],
        columns=['id', 'Complexes']
    )
    df_complexes = df_complexes.drop_duplicates(subset=['Complexes'])
    df_complexes.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"La taille moyenne est {mean_length} | Le nombre de proteines uniques {len(unique_proteins)}")
    # Écrire les résultats dans le fichier de sortie
    
    print(f"Terminé! {len(df_complexes)} complexes uniques ont été écrits dans {output_file}")

# Exécuter la fonction
process_complexes(input_files, output_file)