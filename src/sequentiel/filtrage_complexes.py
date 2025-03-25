import csv

def extract_protein_complexes(input_file, output_file, max_complexes=300, max_proteins=2000):
    unique_proteins = set()
    selected_complexes = []
    
    with open(input_file, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        
        for row in reader:
            if len(selected_complexes) >= max_complexes:
                break
            
            proteins = set(row['subunits_uniprot_id'].split(';')) if row['subunits_uniprot_id'] else set()
            if len(unique_proteins.union(proteins)) <= max_proteins:
                selected_complexes.append((len(selected_complexes) + 1, list(proteins)))
                unique_proteins.update(proteins)
    
    with open(output_file, 'w', encoding='utf-8') as out_file:
        for idx, proteins in selected_complexes:
            out_file.write(f"{idx}\t{';'.join(proteins)}\n")
    
    print(f"Extraction terminée : {len(selected_complexes)} complexes écrits dans {output_file}")

# Exemple d'utilisation
extract_protein_complexes(r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\corum_humanComplexes.txt', r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes.txt')
