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


from collections import defaultdict

# Chemins des fichiers
complexes_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_levure.txt"

reseau_files = [
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\STRING_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\DIP_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\BIOGRID_levure.txt"
]
output_files = [
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_STRING_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_DIP_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_BIOGRID_levure.txt"
]

# 1. Charger les complexes depuis complexes_levure.txt
def load_complexes(file_path):
    complexes = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                proteins = set(parts[1].split())
                complexes.append(proteins)
    return complexes

complexes = load_complexes(complexes_file)
print(f"Nombre total de complexes chargés: {len(complexes)}")

# 2. Charger un réseau PPI depuis un fichier (gère les TSV/CSV et ignore les en-têtes)
def load_ppi_network(file_path):
    network = defaultdict(set)
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("protein1"):  # Ignorer lignes vides/en-têtes
                continue
            parts = line.split('\t') if file_path.endswith('.tsv') else line.split()
            if len(parts) >= 2:
                prot1, prot2 = parts[0], parts[1]  # Prendre les 2 premières colonnes
                network[prot1].add(prot2)
                network[prot2].add(prot1)
    return network

# 3. Vérifier si un complexe est entièrement couvert par le réseau
def is_complex_covered(complexe, network):
    proteins = list(complexe)
    for i in range(len(proteins)):
        for j in range(i + 1, len(proteins)):
            prot1, prot2 = proteins[i], proteins[j]
            if prot2 not in network.get(prot1, set()):
                return False
    return True

# 4. Traiter chaque réseau PPI et sauvegarder les résultats
for reseau_file, output_file in zip(reseau_files, output_files):
    print(f"\nTraitement de {reseau_file}...")
    try:
        network = load_ppi_network(reseau_file)
        
        # Trouver les complexes couverts
        covered_complexes = []
        for complexe in complexes:
            if is_complex_covered(complexe, network):
                covered_complexes.append(complexe)
        
        # Sauvegarder dans le fichier de sortie
        with open(output_file, 'w', encoding='utf-8') as f_out:
            for idx, complexe in enumerate(covered_complexes, 1):
                f_out.write(f"{idx}\t{' '.join(complexe)}\n")
        
        print(f"→ {len(covered_complexes)} complexes trouvés dans ce réseau.")
        print(f"Fichier généré: {output_file}")
    except Exception as e:
        print(f"Erreur avec {reseau_file}: {str(e)}")

print("\nTerminé avec succès !")

#############################################################################################################################

from collections import defaultdict

# Chemins des fichiers
complexes_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_levure.txt"

reseau_files = [
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\weighted_networks\weighted_STRING_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\weighted_networks\weighted_DIP_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\weighted_networks\weighted_BIOGRID_levure.txt"
]
output_files = [
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_STRING_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_DIP_levure.txt",
    r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\complexes_BIOGRID_levure.txt"
]

# 1. Charger les complexes depuis complexes_levure.txt
def load_complexes(file_path):
    complexes = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                proteins = set(parts[1].split())
                complexes.append(proteins)
    return complexes

complexes = load_complexes(complexes_file)
print(f"Nombre total de complexes chargés: {len(complexes)}")

# 2. Charger un réseau PPI depuis un fichier (gère les TSV/CSV et ignore les en-têtes)
def load_ppi_network(file_path):
    network = defaultdict(set)
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("protein1"):  # Ignorer lignes vides/en-têtes
                continue
            parts = line.split('\t') if file_path.endswith('.tsv') else line.split()
            if len(parts) >= 2:
                prot1, prot2 = parts[0], parts[1]  # Prendre les 2 premières colonnes
                network[prot1].add(prot2)
                network[prot2].add(prot1)
    return network

# 3. Vérifier si un complexe est entièrement couvert par le réseau
def is_complex_covered(complexe, network):
    proteins = list(complexe)
    for i in range(len(proteins)):
        for j in range(i + 1, len(proteins)):
            prot1, prot2 = proteins[i], proteins[j]
            if prot2 not in network.get(prot1, set()):
                return False
    return True

# 4. Traiter chaque réseau PPI et sauvegarder les résultats
for reseau_file, output_file in zip(reseau_files, output_files):
    print(f"\nTraitement de {reseau_file}...")
    try:
        network = load_ppi_network(reseau_file)
        
        # Trouver les complexes couverts
        covered_complexes = []
        for complexe in complexes:
            if is_complex_covered(complexe, network):
                covered_complexes.append(complexe)
        
        # Sauvegarder dans le fichier de sortie
        with open(output_file, 'w', encoding='utf-8') as f_out:
            for idx, complexe in enumerate(covered_complexes, 1):
                f_out.write(f"{idx}\t{' '.join(complexe)}\n")
        
        print(f"→ {len(covered_complexes)} complexes trouvés dans ce réseau.")
        print(f"Fichier généré: {output_file}")
    except Exception as e:
        print(f"Erreur avec {reseau_file}: {str(e)}")

print("\nTerminé avec succès !")