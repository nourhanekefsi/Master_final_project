import xml.etree.ElementTree as ET
from pathlib import Path

def process_dip_interactions():
    # Configuration des chemins
    input_file = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\Protein_Interactions\DIP_Interactions.mif25")
    output_file = Path(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\DIP_filtered_interactions.txt")
    
    NS = {'mif': 'http://psi.hupo.org/mi/mif'}
    tree = ET.parse(input_file)
    root = tree.getroot()

    # 1. Extraction des protéines avec vérification
    id_to_protein = {}
    for interactor in root.findall(".//mif:interactor", NS):
        interactor_id = interactor.get("id")
        if not interactor_id:
            continue
            
        uniprot_ref = interactor.find(".//mif:xref/mif:secondaryRef[@db='uniprot knowledge base']", NS)
        refseq_ref = interactor.find(".//mif:xref/mif:secondaryRef[@db='refseq']", NS)
        
        protein_id = None
        if uniprot_ref is not None:
            protein_id = uniprot_ref.get("id")
        elif refseq_ref is not None:
            protein_id = refseq_ref.get("id")
        else:
            short_label = interactor.find(".//mif:names/mif:shortLabel", NS)
            protein_id = short_label.text if short_label is not None else None
        
        if protein_id:
            id_to_protein[interactor_id] = protein_id

    # 2. Extraction des interactions avec contrôle qualité
    unique_interactions = set()
    protein_set = set()  # Pour stocker les protéines uniques
    
    for interaction in root.findall(".//mif:interaction", NS):
        participants = interaction.findall(".//mif:participant/mif:interactorRef", NS)
        
        if len(participants) != 2:
            continue
            
        id1, id2 = participants[0].text, participants[1].text
        
        # Vérification que les deux protéines existent et sont différentes
        if (id1 not in id_to_protein or 
            id2 not in id_to_protein or 
            id_to_protein[id1] == id_to_protein[id2]):
            continue
            
        prot1, prot2 = id_to_protein[id1], id_to_protein[id2]
        
        # Ajout aux protéines uniques
        protein_set.add(prot1)
        protein_set.add(prot2)
        
        # Vérification du score
        score_element = interaction.find(".//mif:confidence/mif:value", NS)
        if score_element is not None:
            try:
                score = float(score_element.text)
                if score <= 0.8:
                    continue
            except (ValueError, TypeError):
                continue
                
        # Ajout sous forme triée pour éviter les doublons A-B vs B-A
        sorted_interaction = tuple(sorted((prot1, prot2)))
        unique_interactions.add(sorted_interaction)

    # 3. Sauvegarde
    with open(output_file, "w") as f:
        f.write("Protein1\tProtein2\n")
        for prot1, prot2 in unique_interactions:
            f.write(f"{prot1}\t{prot2}\n")

    # 4. Calcul et affichage des statistiques
    num_unique_proteins = len(protein_set)
    num_unique_interactions = len(unique_interactions)
    
    print("\nRésultats du traitement:")
    print(f"- Protéines uniques: {num_unique_proteins}")
    print(f"- Interactions uniques: {num_unique_interactions}")
    print(f"Fichier sauvegardé: {output_file}")

if __name__ == "__main__":
    process_dip_interactions()