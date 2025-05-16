from collections import defaultdict
from typing import Dict, Set, Tuple, List

def load_ppi_network(file_path: str) -> Tuple[Dict[str, List[str]], Set[str]]:
    """
    Load PPI network from file.
    Returns:
        - relations: dictionary mapping proteins to their interaction partners
        - proteins_in_network: set of all protein identifiers
    """
    relations = defaultdict(list)
    proteins_in_network = set()
    
    try:
        with open(file_path, 'r') as f:
            # Check and skip header
            first_line = f.readline().strip()
            if first_line.startswith('protein1'):
                pass  # Header found and skipped
            else:
                f.seek(0)  # No header, reset to start
                
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split()
                if len(parts) < 2:
                    continue  # Skip malformed lines
                    
                prot1, prot2 = parts[0], parts[1]
                
                # Add undirected relationship
                relations[prot1].append(prot2)
                relations[prot2].append(prot1)
                
                proteins_in_network.update([prot1, prot2])
                
    except FileNotFoundError:
        raise FileNotFoundError(f"PPI network file not found: {file_path}")
    except Exception as e:
        raise RuntimeError(f"Error loading PPI network: {str(e)}")
        
    return dict(relations), proteins_in_network

def load_go_annotations(file_path: str, proteins_in_network: Set[str]) -> Dict[str, Set[str]]:
    """
    Load GO annotations for proteins in the network.
    Returns dictionary mapping protein IDs to their GO terms.
    """
    go_annotations = defaultdict(set)
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            # Skip header
            header = f.readline()
            
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 4:
                    continue  # Skip malformed lines
                    
                protein = parts[0]
                if protein not in proteins_in_network:
                    continue
                    
                # Extract GO terms more robustly
                go_terms = parts[3]
                for term in go_terms.split(';'):
                    term = term.strip()
                    if term.startswith('GO:'):
                        go_annotations[protein].add(term.split()[0])  # Take first part if multiple
                    elif '[GO:' in term:
                        # Extract GO:XXXXX from format "[GO:XXXXX]"
                        go_start = term.find('[GO:') + 1
                        go_end = term.find(']', go_start)
                        if go_end > go_start:
                            go_annotations[protein].add(term[go_start:go_end])
                            
    except FileNotFoundError:
        raise FileNotFoundError(f"GO annotations file not found: {file_path}")
    except Exception as e:
        raise RuntimeError(f"Error loading GO annotations: {str(e)}")
        
    return dict(go_annotations)

def load_detected_complexes(file_path: str) -> Dict[int, Set[str]]:
    complexes = {}
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):  # Ignore empty lines and comments
                continue
            
            # Split the line into parts
            parts = line.split()
            if len(parts) < 2:  # Skip lines without proteins
                print(f"Warning: Ligne sans protéines ignorée - {line}")
                continue
            
            try:
                complex_id = int(parts[0])
                # Split protein strings by semicolon and create a set
                proteins = set()
                for protein_str in parts[1:]:
                    proteins.update(protein_str.split(';'))
                complexes[complex_id] = proteins
            except ValueError as e:
                print(f"Warning: ID de complexe mal formaté ignoré - {line}")
                continue
                
    return complexes

def load_reference_complexes(file_path: str) -> Dict[int, Set[str]]:
    """Load reference complexes from either:
    1) Tab-separated format: complex_id\tprotein1;protein2;protein3
    2) Space-separated format: complex_id\tprotein1 protein2 protein3
    
    Handles both with and without header row.
    
    Args:
        file_path: Path to the reference complexes file
        
    Returns:
        Dictionary mapping complex IDs to sets of protein identifiers
    """
    complexes = {}
    with open(file_path) as f:
        # First try to read the first line to check for header
        first_line = next(f).strip()
        
        # Check if first line looks like a header (contains 'complex_id' or 'protein')
        has_header = ('complex_id' in first_line.lower()) or ('protein' in first_line.lower())
        
        # If it's a header, we've already consumed it with next(f)
        # If not, we need to process this first line as data
        if not has_header:
            # Process the first line as data
            try:
                parts = first_line.split('\t', 1)
                if len(parts) >= 2:
                    complex_id = int(parts[0])
                    protein_part = parts[1].strip()
                    if ';' in protein_part:
                        proteins = set(protein_part.split(';'))
                    else:
                        proteins = set(protein_part.split())
                    complexes[complex_id] = proteins
            except ValueError as e:
                print(f"Erreur ligne 1: {str(e)} - {first_line}")
        
        # Now process the rest of the file
        for line_num, line in enumerate(f, 2):  # Start counting from line 2
            line = line.strip()
            if not line:
                continue
                
            try:
                parts = line.split('\t', 1)
                if len(parts) < 2:
                    print(f"Warning: Ligne {line_num} mal formatée - {line}")
                    continue
                    
                complex_id = int(parts[0])
                protein_part = parts[1].strip()
                
                if ';' in protein_part:
                    proteins = set(protein_part.split(';'))
                else:
                    proteins = set(protein_part.split())
                    
                complexes[complex_id] = proteins
                
            except ValueError as e:
                print(f"Erreur ligne {line_num}: {str(e)} - {line}")
                continue
                
    return complexes