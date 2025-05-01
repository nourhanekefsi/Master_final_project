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
            if not line or line.startswith('#'):  # Ignore les lignes vides et commentaires
                continue
            
            # Séparation de l'ID et des protéines
            parts = line.split()
            try:
                complex_id = int(parts[0])
                proteins = set(parts[1:])
                complexes[complex_id] = proteins
            except (ValueError, IndexError) as e:
                print(f"Warning: Ligne mal formatée ignorée - {line}")
                continue
                
    return complexes

def load_reference_complexes(file_path: str) -> Dict[int, Set[str]]:
    
    complexes = {}
    with open(file_path) as f:
        header = next(f)  # Ignorer l'en-tête
        
        for line_num, line in enumerate(f, 2):  # Commencer à compter à partir de la ligne 2
            line = line.strip()
            if not line:
                continue
                
            try:
                complex_id, proteins = line.split('\t')
                complex_id = int(complex_id)
                protein_set = set(proteins.split(';'))
                complexes[complex_id] = protein_set
            except ValueError as e:
                print(f"Erreur ligne {line_num}: Format incorrect - {line}")
                continue
                
    return complexes

