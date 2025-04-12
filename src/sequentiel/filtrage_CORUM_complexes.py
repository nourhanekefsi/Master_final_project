import csv

def process_complexes(input_file, output_file):
    # Dictionnaire pour stocker les complexes (id -> set de protéines)
    complexes = {}
    
    # Lire le fichier d'entrée
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            complex_id = row['complex_id']
            proteins = row['subunits_uniprot_id']
            
            if proteins:
                # Séparer les protéines et supprimer les doublons
                protein_list = [p.strip() for p in proteins.split(';') if p.strip()]
                unique_proteins = list(set(protein_list))
                
                # Stocker dans le dictionnaire
                complexes[complex_id] = unique_proteins
    
    # Écrire le fichier de sortie
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['complex_id', 'proteins'])
        
        for complex_id, proteins in complexes.items():
            # Joindre les protéines avec des points-virgules
            protein_str = ';'.join(proteins)
            writer.writerow([complex_id, protein_str])

if __name__ == '__main__':
    input_filename = r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\corum_humanComplexes.txt'  # Remplacez par votre fichier d'entrée
    output_filename = r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\CORUM_complexes_humain.txt'  # Fichier de sortie
    
    process_complexes(input_filename, output_filename)
    print(f"Le fichier de sortie a été généré : {output_filename}")

#############################################################################################################################


import csv
from pathlib import Path

def load_ppi_network(ppi_file):
    """Charge le réseau PPI et retourne un set de protéines uniques"""
    proteins = set()
    try:
        with open(ppi_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    proteins.update({p.strip() for p in row[:2] if p.strip()})
        return proteins
    except Exception as e:
        print(f"Erreur lecture {ppi_file}: {str(e)}")
        return set()

def filter_complexes(complexes_file, ppi_proteins, output_file):
    """Filtre les complexes conservant seulement ceux avec toutes les protéines dans le PPI"""
    stats = {'total': 0, 'kept': 0}
    
    try:
        with open(complexes_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8', newline='') as f_out:
            
            writer = csv.writer(f_out, delimiter='\t')
            writer.writerow(['complex_id', 'proteins'])
            
            for row in csv.reader(f_in, delimiter='\t'):
                if len(row) < 2:
                    continue
                
                stats['total'] += 1
                complex_id, proteins_str = row[0].strip(), row[1].strip()
                proteins = [p.strip() for p in proteins_str.split(';') if p.strip()]
                
                if all(p in ppi_proteins for p in proteins):
                    writer.writerow([complex_id, ';'.join(proteins)])
                    stats['kept'] += 1
        
        return stats
    except Exception as e:
        print(f"Erreur traitement {complexes_file}: {str(e)}")
        return {'total': 0, 'kept': 0}

def process_ppi_network(ppi_file, complexes_file, output_file, network_name):
    """Processus complet pour un réseau PPI"""
    print(f"\nTraitement du réseau {network_name}...")
    
    # Charger les protéines PPI
    ppi_proteins = load_ppi_network(ppi_file)
    if not ppi_proteins:
        print(f"Échec: Aucune protéine chargée depuis {ppi_file}")
        return
    
    print(f"- Protéines uniques dans PPI: {len(ppi_proteins):,}")

    # Filtrer les complexes
    stats = filter_complexes(complexes_file, ppi_proteins, output_file)
    
    if stats['total'] == 0:
        print("Échec: Aucun complexe traité")
        return
    
    # Afficher les statistiques
    print(f"- Complexes analysés: {stats['total']:,}")
    print(f"- Complexes conservés: {stats['kept']:,} ({stats['kept']/stats['total']*100:.1f}%)")
    print(f"- Fichier généré: {output_file}")

if __name__ == '__main__':
    # Configuration des chemins
    DATA_DIR = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data")
    
    # Fichier de complexes source
    CORUM_FILE = DATA_DIR / "complexes" / "CORUM_complexes_humain.txt"
    
    # Traitement BIOGRID
    process_ppi_network(
        ppi_file=DATA_DIR / "weighted_networks" / "weighted_BIOGRID_humain.txt",
        complexes_file=CORUM_FILE,
        output_file=DATA_DIR / "complexes" / "BIOGRID_complexes_humain.txt",
        network_name="BIOGRID"
    )
    
    # Traitement STRING
    process_ppi_network(
        ppi_file=DATA_DIR / "weighted_networks" / "weighted_STRING_humain.txt",
        complexes_file=CORUM_FILE,
        output_file=DATA_DIR / "complexes" / "STRING_complexes_humain.txt",
        network_name="STRING"
    )