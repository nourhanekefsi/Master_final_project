import csv
import pandas as pd
import numpy as np

def process_corum(file_path):
    complexes = {}
    all_proteins = set()
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        mean_length = []
        for row in reader:
            complex_id = int(row['complex_id'])
            proteins = row['subunits_uniprot_id']
            if isinstance(proteins, str) and proteins:
                protein_list = sorted(set(p.strip() for p in proteins.split(';') if p.strip()))
                if len(protein_list)>2:
                    complexes[complex_id] = protein_list
                    all_proteins.update(protein_list)
                    mean_length.append(len(protein_list))
    mean_length = np.mean(mean_length)
    print(f"[CORUM] Complexes: {len(complexes)} | Protéines uniques: {len(all_proteins)} | taille moyenne {mean_length}")
    return complexes, all_proteins

def process_humap(file_path):
    complexes = {}
    all_proteins = set()
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        complex_id = 0
        mean_length = []
        for row in reader:
            proteins = row['Uniprot_ACCs']
            if proteins:
                protein_list = sorted(set(p.strip() for p in proteins.split(' ') if p.strip()))
                if len(protein_list)>2:
                    complexes[complex_id] = protein_list
                    all_proteins.update(protein_list)
                    complex_id += 1
                    mean_length.append(len(protein_list))
    mean_length = np.mean(mean_length)
    print(f"[Hu.MAP] Complexes: {len(complexes)} | Protéines uniques: {len(all_proteins)} | taille moyenne {mean_length}")
    return complexes, all_proteins

def export_unique_complexes(complexes_dict, output_path):
    seen = set()
    unique_complexes = []
    complex_id = 0
    mean_length = []
    for id, proteins in complexes_dict.items():
        protein_tuple = tuple(sorted(proteins))
        if protein_tuple not in seen:
            seen.add(protein_tuple)
            unique_complexes.append((complex_id, protein_tuple))
            complex_id += 1
            mean_length.append(len(protein_tuple))
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['complex_id', 'proteins'])
        for complex_id, protein_tuple in unique_complexes:
            writer.writerow([complex_id, ' '.join(protein_tuple)])
        
    mean_length = np.mean(mean_length)
    all_proteins = set(p for _, proteins in unique_complexes for p in proteins)
    print(f"\n[Union] Complexes non redondants : {len(unique_complexes)}")
    print(f"[Union] Protéines uniques totales : {len(all_proteins)}")
    print(f"la taille moyenne des complexes humain totale {mean_length}")
    print(f"Fichier sauvegardé : {output_path}")

if __name__ == '__main__':
    corum_file = r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\corum_humanComplexes.txt'
    humap_file = r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\hu.MAP_complexes.txt'
    output_file = r'C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\merged_unique_complexes.txt'

    corum_complexes, corum_proteins = process_corum(corum_file)
    humap_complexes, humap_proteins = process_humap(humap_file)

    all_complexes = {**corum_complexes, **humap_complexes}
    export_unique_complexes(all_complexes, output_file)

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
    if len(proteins) < 2:  # Un complexe d'une seule protéine est toujours connecté
        return True
    
    visited = set()
    queue = deque()
    
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
    """Filtre les complexes selon les critères spécifiés"""
    stats = {
        'total': 0,
        'kept': 0,
        'small': 0,
        'missing_proteins': 0,
        'disconnected': 0
    }
    
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
                proteins_set = set(proteins)
                
                # Vérifier la taille du complexe
                if len(proteins_set) < 3:
                    stats['small'] += 1
                    continue
                
                # Vérifier que toutes les protéines sont dans le PPI
                if not all(p in ppi_proteins for p in proteins_set):
                    stats['missing_proteins'] += 1
                    continue
                
                # Vérifier la connectivité
                if not is_single_connected_component(proteins_set, ppi_graph):
                    stats['disconnected'] += 1
                    continue
                
                # Écrire le complexe valide
                writer.writerow([complex_id, ';'.join(proteins)])
                stats['kept'] += 1
        
        return stats
    except Exception as e:
        print(f"Erreur traitement {complexes_file}: {str(e)}")
        return {'total': 0, 'kept': 0, 'small': 0, 'missing_proteins': 0, 'disconnected': 0}

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
    print(f"  - Rejetés (taille < 3): {stats['small']:,}")
    print(f"  - Rejetés (protéines manquantes): {stats['missing_proteins']:,}")
    print(f"  - Rejetés (non connectés): {stats['disconnected']:,}")
    print(f"- Fichier généré: {output_file}")

if __name__ == '__main__':
    # Configuration des chemins
    DATA_DIR = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data")
    
    # Fichier de complexes source
    CORUM_FILE = DATA_DIR / "complexes" / "merged_unique_complexes.txt"
    
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