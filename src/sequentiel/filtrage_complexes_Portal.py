import os
import zipfile
import shutil
from xml.etree import ElementTree as ET

# Chemins
zip_path = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\yeast.zip"
extract_folder = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Complexes\yeast_extracted"
output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\Portal_complexes_levure.txt"

# 1. Nettoyage et extraction
if os.path.exists(extract_folder):
    shutil.rmtree(extract_folder)
os.makedirs(extract_folder, exist_ok=True)

with zipfile.ZipFile(zip_path, 'r') as zip_ref:
    zip_ref.extractall(extract_folder)

# 2. Traitement spécifique pour votre format XML
complexes = []
yeast_folder = os.path.join(extract_folder, "yeast")

for filename in os.listdir(yeast_folder):
    if filename.endswith(".xml"):
        filepath = os.path.join(yeast_folder, filename)
        
        try:
            tree = ET.parse(filepath)
            root = tree.getroot()
            
            # Namespace spécifique à vos fichiers
            ns = {'mif': 'http://psi.hupo.org/mi/mif300'}
            
            # Recherche des interactions complexes
            for interaction in root.findall(".//mif:abstractInteraction", ns):
                proteins = set()
                
                # Recherche des participants
                for participant in interaction.findall(".//mif:participant", ns):
                    # Référence à l'interacteur
                    interactor_ref = participant.find(".//mif:interactorRef", ns)
                    if interactor_ref is not None:
                        # Trouver l'interacteur correspondant
                        interactor = root.find(f".//mif:interactor[@id='{interactor_ref.text}']", ns)
                        if interactor is not None:
                            # Vérifier si c'est une protéine
                            interactor_type = interactor.find(".//mif:interactorType/mif:names/mif:shortLabel", ns)
                            if interactor_type is not None and interactor_type.text == "protein":
                                # Récupérer l'identifiant UniProt
                                uniprot = interactor.find(".//mif:xref/mif:primaryRef[@db='uniprotkb']", ns)
                                if uniprot is not None:
                                    proteins.add(uniprot.get("id"))
                
                if proteins:
                    complexes.append(sorted(proteins))
        
        except Exception as e:
            print(f"Erreur avec {filename}: {str(e)[:200]}")

# 3. Écriture du fichier final
with open(output_file, 'w', encoding='utf-8') as f_out:
    # Format: ID [tab] Liste_de_protéines (séparées par des espaces)
    for idx, proteins in enumerate(complexes, 1):
        f_out.write(f"{idx}\t{' '.join(proteins)}\n")

# 4. Rapport
print(f"Fichier généré: {output_file}")
print(f"Nombre total de complexes trouvés: {len(complexes)}")


################################################################################################################################


from collections import defaultdict, deque
from pathlib import Path

def load_complexes(file_path):
    """Charge les complexes depuis le fichier"""
    complexes = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                complex_id = parts[0]
                proteins = set(p.strip() for p in parts[1].split() if p.strip())
                if len(proteins) >= 3:  # Ne garder que les complexes avec ≥3 protéines
                    complexes.append((complex_id, proteins))
    return complexes

def load_ppi_network(file_path):
    """Charge le réseau PPI"""
    proteins = set()
    graph = defaultdict(set)
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("protein1"):
                continue
            parts = line.split('\t') if '\t' in line else line.split()
            if len(parts) >= 2:
                p1, p2 = parts[0].strip(), parts[1].strip()
                if p1 and p2:
                    proteins.add(p1)
                    proteins.add(p2)
                    graph[p1].add(p2)
                    graph[p2].add(p1)
    return proteins, graph

def is_single_connected_component(proteins, ppi_graph):
    """Vérifie la connectivité du complexe"""
    if len(proteins) < 2:
        return True
    
    visited = set()
    queue = deque([next(iter(proteins))])
    visited.add(queue[0])
    
    while queue:
        current = queue.popleft()
        for neighbor in ppi_graph[current]:
            if neighbor in proteins and neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    
    return visited == proteins

def filter_and_save_complexes(complexes, ppi_proteins, ppi_graph, output_file):
    """Filtre et sauvegarde les complexes valides"""
    stats = {
        'total': len(complexes),
        'kept': 0,
        'missing_proteins': 0,
        'disconnected': 0
    }
    
    with open(output_file, 'w', encoding='utf-8') as f_out:
        for complex_id, proteins in complexes:
            # Vérifier présence dans PPI
            if not all(p in ppi_proteins for p in proteins):
                stats['missing_proteins'] += 1
                continue
            
            # Vérifier connectivité
            if not is_single_connected_component(proteins, ppi_graph):
                stats['disconnected'] += 1
                continue
            
            # Écrire le complexe valide
            f_out.write(f"{complex_id}\t{' '.join(proteins)}\n")
            stats['kept'] += 1
    
    return stats

def main():
    # Configuration des chemins
    base_dir = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data")
    
    # Fichiers d'entrée/sortie
    complexes_file = base_dir / "complexes" / "complexes_levure.txt"
    reseau_files = [
        base_dir / "weighted_networks" / "weighted_STRING_levure.txt",
        base_dir / "weighted_networks" / "weighted_DIP_levure.txt",
        base_dir / "weighted_networks" / "weighted_BIOGRID_levure.txt"
    ]
    output_files = [
        base_dir / "complexes" / "complexes_STRING_levure.txt",
        base_dir / "complexes" / "complexes_DIP_levure.txt",
        base_dir / "complexes" / "complexes_BIOGRID_levure.txt"
    ]

    # Charger les complexes (déjà filtrés pour taille ≥3)
    complexes = load_complexes(complexes_file)
    print(f"Complexes chargés (taille ≥3): {len(complexes)}")

    # Traiter chaque réseau PPI
    for ppi_file, out_file in zip(reseau_files, output_files):
        print(f"\nTraitement de {ppi_file.name}...")
        
        try:
            # Charger le réseau PPI
            ppi_proteins, ppi_graph = load_ppi_network(ppi_file)
            print(f"- Protéines uniques: {len(ppi_proteins):,}")
            print(f"- Interactions: {sum(len(v) for v in ppi_graph.values())//2:,}")

            # Filtrer et sauvegarder
            stats = filter_and_save_complexes(complexes, ppi_proteins, ppi_graph, out_file)
            
            # Statistiques
            print(f"- Complexes analysés: {stats['total']:,}")
            print(f"- Complexes conservés: {stats['kept']:,} ({stats['kept']/stats['total']*100:.1f}%)")
            print(f"  - Rejetés (protéines manquantes): {stats['missing_proteins']:,}")
            print(f"  - Rejetés (non connectés): {stats['disconnected']:,}")
            print(f"- Fichier généré: {out_file}")

        except Exception as e:
            print(f"Erreur: {str(e)}")

    print("\nTerminé!")

if __name__ == "__main__":
    main()