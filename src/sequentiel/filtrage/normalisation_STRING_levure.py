import os
from collections import defaultdict

# Chemins des fichiers
interactions_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\STRING_levure.temp"
mapping_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\autres\YEAST_559292_idmapping.dat"
output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\STRING_levure.txt"

def create_mapping_dict(mapping_file):
    """Crée un dictionnaire de mapping STRING -> UniProt et inversement"""
    str_to_uniprot = defaultdict(list)
    uniprot_to_str = defaultdict(list)
    
    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                uniprot_id = parts[0]
                db_type = parts[1]
                external_id = parts[2]

                if db_type == "STRING":
                    # Supprimer le préfixe '4932.' si présent
                    if external_id.startswith("4932."):
                        external_id = external_id.replace("4932.", "")
                    str_to_uniprot[external_id].append(uniprot_id)
                    uniprot_to_str[uniprot_id].append(external_id)
    
    return str_to_uniprot

def process_interactions(interactions_file, str_to_uniprot):
    """Traite le fichier d'interactions"""
    interactions = []
    missing_mappings = set()
    
    with open(interactions_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                prot1, prot2 = parts[0], parts[1]
                
                # Vérifier si les protéines existent dans le mapping
                mapped1 = str_to_uniprot.get(prot1, [None])[0]
                mapped2 = str_to_uniprot.get(prot2, [None])[0]
                
                if mapped1 and mapped2:
                    interactions.append((mapped1, mapped2))
                else:
                    if not mapped1:
                        missing_mappings.add(prot1)
                    if not mapped2:
                        missing_mappings.add(prot2)
    
    return interactions

def main():
    # 1. Créer le dictionnaire de mapping
    str_to_uniprot = create_mapping_dict(mapping_file)
    
    # 2. Traiter les interactions
    filtered_interactions = process_interactions(interactions_file, str_to_uniprot)
    
    # 3. Sauvegarder les résultats
    with open(output_file, 'w') as f:
        for prot1, prot2 in filtered_interactions:
            f.write(f"{prot1}\t{prot2}\n")

    # 4. Calcul des statistiques
    unique_proteins = set()
    for p1, p2 in filtered_interactions:
        unique_proteins.update([p1, p2])

    # 5. Affichage simplifié
    print("\nTraitement terminé avec succès.")
    print(f"Nombre total d'interactions conservées : {len(filtered_interactions)}")
    print(f"Nombre total de protéines uniques : {len(unique_proteins)}")


if __name__ == "__main__":
    main()
