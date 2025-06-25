from collections import deque
import random
from typing import Dict, List, Set
from population import MyPopulation
import numpy as np
import pandas as pd
from chromosome import Chromosome
from multiprocessing.pool import ThreadPool
import logging
logger = logging.getLogger(__name__)

def load_ppi(ppi_file: str) -> Dict[str, Dict[str, float]]:
    """Load protein-protein interaction network"""
    try:
       graph_df = pd.read_csv(ppi_file, sep="\t", header=0, names=["node1", "node2", "weight"])
    except Exception as e:
        logger.error(f"Erreur lors du chargement du fichier PPI: {e}")
        
    graph_df["weight"] = pd.to_numeric(graph_df["weight"], errors='coerce')
    graph_df = graph_df.dropna(subset=["weight"])
    
    ppi_graph = {}
    for _, row in graph_df.iterrows():
        node1, node2, weight = str(row['node1']), str(row['node2']), float(row['weight'])
        
        if node1 not in ppi_graph:
            ppi_graph[node1] = {}
        ppi_graph[node1][node2] = weight
        
        if node2 not in ppi_graph:
            ppi_graph[node2] = {}
        ppi_graph[node2][node1] = weight
        
    return ppi_graph

def load_complexes(file_path: str) -> List[List[str]]:
    """Charge les complexes de référence depuis le fichier."""
    try:
        complexes = pd.read_csv(file_path, sep='\t', header=0, names=['complex_id', 'proteins'])
        return [row.split() for row in complexes['proteins'].tolist()]
    except Exception as e:
        logger.error(f"Erreur lors du chargement des complexes de references: {e}")
        return []
    
import random
from collections import deque
from typing import List, Dict, Set
from multiprocessing import Pool, cpu_count
import numpy as np

_shared_graph = None
_shared_protein_nodes = None
_shared_min_size = 3
_shared_max_size = 25
_shared_max_complexes = 198

def init_globals(graph, min_size, max_size):
    global _shared_graph, _shared_protein_nodes, _shared_min_size, _shared_max_size
    _shared_graph = graph
    _shared_protein_nodes = list(graph.keys())
    _shared_min_size = min_size
    _shared_max_size = max_size

def create_connected_complex() -> Set[str]:
    if not _shared_protein_nodes:
        return set()
    start_node = random.choice(_shared_protein_nodes)
    while len(_shared_graph.get(start_node, {})) == 0:
        start_node = random.choice(_shared_protein_nodes)
    complex_nodes = set()
    queue = deque([start_node])
    target_size = random.randint(_shared_min_size, _shared_max_size)
    while queue and len(complex_nodes) < target_size:
        current = queue.popleft()
        if current in complex_nodes:
            continue
        complex_nodes.add(current)
        neighbors = sorted(
            [(n, w) for n, w in _shared_graph.get(current, {}).items() if n not in complex_nodes],
            key=lambda x: -x[1]
        )
        for neighbor, _ in neighbors[:3]:
            if len(complex_nodes) >= target_size:
                break
            queue.append(neighbor)
    return complex_nodes

def generate_individual(_):
    """Fonction utilisée par multiprocessing pour générer un individu"""
    try:
        ref_mean = _shared_max_complexes
        num_complexes = random.randint(int(0.75 * ref_mean), int(1.25 * ref_mean))
        complexes = []
        used_proteins = set()
        for _ in range(num_complexes):
            new_complex = create_connected_complex()
            if new_complex and len(new_complex & used_proteins) / len(new_complex) < 0.2:
                complexes.append(new_complex)
                used_proteins.update(new_complex)
        if complexes:
            return Chromosome(_shared_graph, complexes)
    except Exception as e:
        logger.error(f"Erreur dans generate_individual: {e}")
    return None

def random_initialization(reference_complexes: List[List[str]], 
                         graph: Dict[str, Dict[str, float]], 
                         population_size: int = 26,
                         max_complexes: int = 198, fact: int = 6, min_size=3, max_size=25) -> List:
    """
    Génère une population avec des complexes protéiques aléatoires et connectés
    Parallélisée avec multiprocessing
    """
    logger.info("Début de la génération de population initiale")
    init_globals(graph, min_size, max_size)
    global _shared_max_complexes
    _shared_max_complexes = max_complexes
    nb_solutions = fact * population_size

    with ThreadPool(processes=cpu_count()) as pool:
        results = pool.map(generate_individual, range(nb_solutions))

    individuals = [ind for ind in results if ind is not None]
    if not individuals:
        raise ValueError("Aucun individu valide n'a été généré.")

    fitnesses = [ind._fitness for ind in individuals]
    indices = np.argsort(fitnesses)[::-1]
    best_individuals = [individuals[i] for i in indices[:population_size]]
    
    return best_individuals


def load_initialization_complexes(file_path: str, graph : Dict[str, Dict[str,float]]) -> List[List[str]]:
    """Charge les complexes de référence depuis le fichier."""
    try:
        complexes = pd.read_csv(file_path, sep='\t', header=0, names=['solutionID','complexID','proteins'])
        id_s = 1
        solution = []
        population = []
        for _, row  in complexes.iterrows():
            s_id, complex = row['solutionID'], row['proteins']
            if s_id == id_s:
                solution.append(set(str(complex).split(' ')))
            else:
                id_s += 1 
                population.append(Chromosome(graph, solution))
                solution = []
        population.append(Chromosome(graph, solution))
        return population
    except Exception as e:
        logger.error(f"Erreur lors du chargement des complexes d'initialisation: {e}")
        return []

def save_complexes(complexes: List[Set[str]], filename: str) -> pd.DataFrame:
    data = []
    for i, complex in enumerate(complexes, 1):
        proteins = sorted(complex)
        data.append({'complex_id': i, 'proteins': ' '.join(proteins)})
    
    df = pd.DataFrame(data)
    df.to_csv(filename, sep='\t', index=False)
    return df
