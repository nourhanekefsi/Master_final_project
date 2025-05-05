import sys
from collections import defaultdict

def parse_mapping_file(mapping_file):
    """
    Parse the mapping file to create a dictionary from ENSP to Uniprot IDs.
    The mapping file is assumed to have lines with format like:
    P31946 STRING 9606.ENSP00000361930
    """
    ensembl_to_uniprot = defaultdict(list)
    
    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            uniprot_id, db, db_id = parts[:3]
            if db == "STRING":
                # Extract ENSP ID from STRING format (e.g., 9606.ENSP00000361930)
                if '.' in db_id:
                    ensembl_id = db_id.split('.')[1]
                    ensembl_to_uniprot[ensembl_id].append(uniprot_id)
    
    # Remove duplicates and convert to single value (take first Uniprot ID if multiple)
    return {ensp: uniprots[0] for ensp, uniprots in ensembl_to_uniprot.items() if uniprots}

def convert_interactions(input_file, output_file, mapping_dict):
    """
    Convert interaction file from ENSP IDs to Uniprot IDs using the mapping dictionary.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue
            prot1, prot2 = parts
            uniprot1 = mapping_dict.get(prot1, prot1)  # Keep original if not found
            uniprot2 = mapping_dict.get(prot2, prot2)
            outfile.write(f"{uniprot1}\t{uniprot2}\n")

def main():
    
    interactions_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\STRING_humain_filtered_interactions.txt"
    mapping_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\raw data\autres\HUMAN_9606_idmapping.dat"
    output_file = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\interactions\STRING_humain.txt"
    
    print("Parsing mapping file...")
    mapping_dict = parse_mapping_file(mapping_file)
    print(f"Found {len(mapping_dict)} ENSP to Uniprot mappings")
    
    print("Converting interactions...")
    convert_interactions(interactions_file, output_file, mapping_dict)
    print(f"Converted interactions written to {output_file}")

if __name__ == "__main__":
    main()