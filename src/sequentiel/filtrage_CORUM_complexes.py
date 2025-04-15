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
from collections import defaultdict, deque

def load_ppi_network(ppi_file):
    """Charge le réseau PPI et retourne un set de protéines uniques et le graphe PPI"""
    proteins = set()
    graph = defaultdict(set)
    try:
        with open(ppi_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) >= 2:
                    p1, p2 = row[0].strip(), row[1].strip()
                    if p1 and p2:
                        proteins.add(p1)
                        proteins.add(p2)
                        graph[p1].add(p2)
                        graph[p2].add(p1)
        return proteins, graph
    except Exception as e:
        print(f"Erreur lecture {ppi_file}: {str(e)}")
        return set(), defaultdict(set)

def is_single_connected_component(proteins, ppi_graph):
    """Vérifie si les protéines forment un seul composant connecté dans le réseau PPI"""
    if not proteins:
        return False
    
    visited = set()
    queue = deque()
    
    # Prendre une protéine quelconque comme point de départ
    start_protein = next(iter(proteins))
    queue.append(start_protein)
    visited.add(start_protein)
    
    while queue:
        current = queue.popleft()
        for neighbor in ppi_graph[current]:
            if neighbor in proteins and neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    
    return visited == proteins

def filter_complexes(complexes_file, ppi_proteins, ppi_graph, output_file):
    """Filtre les complexes conservant seulement ceux avec toutes les protéines dans le PPI et formant un seul composant connecté"""
    stats = {'total': 0, 'kept': 0, 'missing_proteins': 0, 'disconnected': 0}
    
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
                proteins = set(p.strip() for p in proteins_str.split(';') if p.strip())
                
                # Vérifier que toutes les protéines sont dans le PPI
                if not all(p in ppi_proteins for p in proteins):
                    stats['missing_proteins'] += 1
                    continue
                
                # Vérifier que les protéines forment un seul composant connecté
                if not is_single_connected_component(proteins, ppi_graph):
                    stats['disconnected'] += 1
                    continue
                
                # Écrire le complexe valide
                writer.writerow([complex_id, ';'.join(proteins)])
                stats['kept'] += 1
        
        return stats
    except Exception as e:
        print(f"Erreur traitement {complexes_file}: {str(e)}")
        return {'total': 0, 'kept': 0, 'missing_proteins': 0, 'disconnected': 0}

def process_ppi_network(ppi_file, complexes_file, output_file, network_name):
    """Processus complet pour un réseau PPI"""
    print(f"\nTraitement du réseau {network_name}...")
    
    # Charger les protéines PPI et le graphe
    ppi_proteins, ppi_graph = load_ppi_network(ppi_file)
    if not ppi_proteins:
        print(f"Échec: Aucune protéine chargée depuis {ppi_file}")
        return
    
    print(f"- Protéines uniques dans PPI: {len(ppi_proteins):,}")
    print(f"- Interactions dans PPI: {sum(len(v) for v in ppi_graph.values())//2:,}")

    # Filtrer les complexes
    stats = filter_complexes(complexes_file, ppi_proteins, ppi_graph, output_file)
    
    if stats['total'] == 0:
        print("Échec: Aucun complexe traité")
        return
    
    # Afficher les statistiques
    print(f"- Complexes analysés: {stats['total']:,}")
    print(f"- Complexes conservés: {stats['kept']:,} ({stats['kept']/stats['total']*100:.1f}%)")
    print(f"  - Rejetés (protéines manquantes): {stats['missing_proteins']:,}")
    print(f"  - Rejetés (non connectés): {stats['disconnected']:,}")
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