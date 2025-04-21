import os
import zipfile
import shutil
from xml.etree import ElementTree as ET
from collections import defaultdict
import panda as pd
from unidecode import unidecode

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

# 2. Traitement des fichiers XML
complexes = []
all_proteins = set()
protein_to_complexes = defaultdict(list)

yeast_folder = os.path.join(extract_folder, "yeast")

for filename in os.listdir(yeast_folder):
    if filename.endswith(".xml"):
        filepath = os.path.join(yeast_folder, filename)
        
        try:
            tree = ET.parse(filepath)
            root = tree.getroot()
            ns = {'mif': 'http://psi.hupo.org/mi/mif300'}
            
            for interaction in root.findall(".//mif:abstractInteraction", ns):
                proteins = set()
                
                for participant in interaction.findall(".//mif:participant", ns):
                    interactor_ref = participant.find(".//mif:interactorRef", ns)
                    if interactor_ref is not None:
                        interactor = root.find(f".//mif:interactor[@id='{interactor_ref.text}']", ns)
                        if interactor is not None:
                            interactor_type = interactor.find(".//mif:interactorType/mif:names/mif:shortLabel", ns)
                            if interactor_type is not None and interactor_type.text == "protein":
                                uniprot = interactor.find(".//mif:xref/mif:primaryRef[@db='uniprotkb']", ns)
                                if uniprot is not None:
                                    protein_id = uniprot.get("id")
                                    proteins.add(protein_id)
                                    all_proteins.add(protein_id)
                
                if len(proteins)>2:
                    complex_id = f"{len(complexes)+1}"
                    complexes.append((complex_id, proteins))
                    for protein in proteins:
                        protein_to_complexes[protein].append(complex_id)
        
        except Exception as e:
            print(f"Erreur avec {filename}: {str(e)[:200]}")

# Flatten the list of complexes into a DataFrame
results = pd.DataFrame([
    {"id": complex_id, "Proteines": ';'.join(sorted(proteins))}
    for complex_id, proteins in complexes
])

# Remove duplicates based on the 'Proteines' column
results = results.drop_duplicates(subset=['Proteines'])

# Save the results to a file
results.to_csv(output_file, sep="\t", index=False, header=False)
"""# 3. Écriture du fichier final
with open(output_file, 'w', encoding='utf-8') as f_out:
    for complex_id, proteins in complexes:
        f_out.write(f"{complex_id}\t{' '.join(sorted(proteins))}\n")
"""
# 4. Calcul des statistiques
num_complexes = len(complexes)
num_unique_proteins = len(all_proteins)

# Distribution des tailles de complexes
size_distribution = defaultdict(int)
for _, proteins in complexes:
    size_distribution[len(proteins)] += 1

# Protéines les plus fréquentes
protein_frequency = {protein: len(complexes) for protein, complexes in protein_to_complexes.items()}
top_proteins = sorted(protein_frequency.items(), key=lambda x: x[1], reverse=True)[:10]

# 5. Affichage des résultats
print("\n" + "="*50)
print(f"Fichier généré: {output_file}")
print("="*50)
print(f"\nSTATISTIQUES GLOBALES:")
print(f"- Nombre total de complexes: {num_complexes}")
print(f"- Nombre de protéines uniques: {num_unique_proteins}")
print(f"- Nombre moyen de protéines par complexe: {sum(len(p) for _, p in complexes)/num_complexes:.2f}")

print("\nDISTRIBUTION DES TAILLES DE COMPLEXES:")
for size, count in sorted(size_distribution.items()):
    print(f"- Complexes de taille {size}: {count} ({(count/num_complexes)*100:.1f}%)")

print("\nTOP 10 DES PROTÉINES LES PLUS FRÉQUENTES:")
for protein, freq in top_proteins:
    print(f"- {protein}: présente dans {freq} complexes ({(freq/num_complexes)*100:.1f}%)")
print("="*50)

################################################################################################################################

from collections import defaultdict, deque
from pathlib import Path

def load_complexes(file_path):
    """Charge les complexes depuis le fichier et retourne une liste de tuples (id_complexe, set de protéines)"""
    complexes = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                complex_id = parts[0]
                proteins = set(p.strip() for p in parts[1].split() if p.strip())
                complexes.append((complex_id, proteins))
    return complexes

def load_ppi_network(file_path):
    """Charge le réseau PPI et retourne un set de protéines et le graphe PPI"""
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
    """Vérifie si les protéines forment un seul composant connecté dans le réseau PPI"""
    if not proteins:
        return False
    
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

def filter_complexes(complexes, ppi_proteins, ppi_graph, output_file):
    """Filtre les complexes et sauvegarde ceux valides"""
    stats = {
        'total': 0,
        'kept': 0,
        'missing_proteins': 0,
        'disconnected': 0,
        'too_small': 0  # Ajout d'une nouvelle statistique pour les complexes trop petits
    }
    
    with open(output_file, 'w', encoding='utf-8') as f_out:
        for complex_id, proteins in complexes:
            stats['total'] += 1
            
            # Étape 1: Filtrer les complexes avec moins de 3 protéines
            if len(proteins) < 3:
                stats['too_small'] += 1
                continue
            
            # Vérifier que toutes les protéines sont dans le PPI
            if not all(p in ppi_proteins for p in proteins):
                stats['missing_proteins'] += 1
                continue
            
            # Vérifier la connectivité
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

    # Charger les complexes
    complexes = load_complexes(complexes_file)
    print(f"Nombre total de complexes chargés: {len(complexes)}")

    # Traiter chaque réseau PPI
    for ppi_file, out_file in zip(reseau_files, output_files):
        print(f"\nTraitement de {ppi_file.name}...")
        
        try:
            # Charger le réseau PPI
            ppi_proteins, ppi_graph = load_ppi_network(ppi_file)
            print(f"- Protéines uniques dans PPI: {len(ppi_proteins):,}")
            print(f"- Interactions dans PPI: {sum(len(v) for v in ppi_graph.values())//2:,}")

            # Filtrer les complexes
            stats = filter_complexes(complexes, ppi_proteins, ppi_graph, out_file)
            
            # Afficher les statistiques
            print(f"- Complexes analysés: {stats['total']:,}")
            print(f"- Complexes conservés: {stats['kept']:,} ({stats['kept']/stats['total']*100:.1f}%)")
            print(f"  - Rejetés (protéines manquantes): {stats['missing_proteins']:,}")
            print(f"  - Rejetés (non connectés): {stats['disconnected']:,}")
            print(f"  - Rejetés (taille < 3): {stats['too_small']:,}")
            print(f"- Fichier généré: {out_file}")

        except Exception as e:
            print(f"Erreur avec {ppi_file}: {str(e)}")

    print("\nTerminé avec succès !")

if __name__ == "__main__":
    main()