import random
import os
import sys
from typing import Any, Callable, Dict, List, Optional, Set
from HPC_GA.core import GeneticAlgorithm
from population import MyPopulation 
from chromosome import Chromosome
from HPC_GA.core.operators import Crossover, Mutator
from HPC_Tabu.parallel import ParallelTabuSearch
import numpy as np
import pandas as pd
import ray 
from metrics import FF
from parallel_metrics import overlap_score
from chromosome import Chromosome
from solution import Solution
from HPC_Tabu.common.neighborhood import NeighborhoodGenerator
from HPC_Tabu.sequential import TabuSearch
from parallel_metrics import MasterSlaveMetrics
import networkx as nx
from sklearn.cluster import SpectralClustering
from multiprocessing.pool import ThreadPool
import networkx as nx
import logging
from data_loader import save_complexes

def setup_worker_logger(type, name=None):
    """Initialise le logger pour chaque worker avec un fichier de log unique"""
    log_dir=f"kefsi_mekhazni_workspace/results/{name}"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{type}_{os.getpid()}.log")

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Supprimer les anciens handlers pour éviter les doublons
    if logger.hasHandlers():
        logger.handlers.clear()

    handler = logging.FileHandler(log_file, mode='a')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Option 1 : forcer le flush à chaque log (très utile pour Ray ou multiprocessing)
    class FlushOnWrite:
        def __init__(self, stream):
            self.stream = stream
        def write(self, msg):
            self.stream.write(msg)
            self.stream.flush()
        def flush(self):
            self.stream.flush()

    sys.stdout = FlushOnWrite(open(log_file, "a", buffering=1))
    sys.stderr = FlushOnWrite(open(log_file, "a", buffering=1))

    logger.info(f"Logger initialisé pour le worker PID {os.getpid()}")


class Neighborhood(NeighborhoodGenerator):
    """Générateur de voisinage parallélisé avec des opérations intelligentes"""
    
    def __init__(self, min_complex_size: int = 3, max_complex_size: int = 15, n_workers: int = None):
        self.id = random.randint(0, 1000000000000000000)
        self.min_complex_size = min_complex_size
        self.max_complex_size = max_complex_size
        self.n_workers = 8
        
        
    def generate(self, solution: Solution) -> List[Solution]:
        
        #print(f"  Taille du complexe initial: {len(solution.representation)}")
         
        
        # Créer une liste de toutes les fonctions de génération de voisins
        neighbor_functions = [
            self._generate_add_protein_neighbors,
            self._generate_remove_protein_neighbors,
            self._generate_swap_neighbors,
            self._generate_merge_neighbors,
            self._generate_split_neighbors,
            self._generate_density_based_neighbors
        ]
        
        with ThreadPool(self.n_workers) as pool:
            results = pool.starmap(
                self._execute_neighbor_function,
                [(func, solution) for func in neighbor_functions]
            )
        
        
        # Fusionner tous les résultats
        neighbors = [neighbor for sublist in results for neighbor in sublist]
        
        # Filtrer les doublons et les complexes de taille invalide
        unique_neighbors = []
        seen = set()
        for neighbor in neighbors:
            neighbor_complex = frozenset(neighbor.representation)
            size = len(neighbor_complex)
            if (neighbor_complex not in seen and 
                self.min_complex_size <= size <= self.max_complex_size):
                seen.add(neighbor_complex)
                unique_neighbors.append(neighbor)
    
        #print(f"Généré {len(unique_neighbors)} voisins uniques")
         
        return unique_neighbors

    def _execute_neighbor_function(self, func: Callable, solution: Solution) -> List[Solution]:
        """Wrapper pour exécuter une fonction de génération de voisins"""
        try:
            return func(solution)
        except Exception as e:
            #print(f"Erreur dans {func.__name__}: {e}")
             
            return []

    def _generate_add_protein_neighbors(self, solution: Solution) -> List[Solution]:
        """Ajoute des protéines fortement connectées au complexe"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        # Calculer la connectivité moyenne avec les protéines externes
        external_proteins = list(set(graph.keys()) - set(current_complex))
        chunks = np.array_split(external_proteins, self.n_workers)
        
        with ThreadPool(self.n_workers) as pool:
            results = pool.starmap(
                self._process_add_protein_chunk,
                [(current_complex, graph, chunk) for chunk in chunks if len(chunk) > 0]
            )
        
        for chunk_neighbors in results:
            neighbors.extend(chunk_neighbors)
            
        return neighbors

    def _process_add_protein_chunk(self, current_complex: Set[str], graph: Dict, proteins: List[str]) -> List[Solution]:
        """Traite un chunk de protéines pour l'ajout"""
        chunk_neighbors = []
        protein_scores = []
        
        for protein in proteins:
            interactions = []
            for p in current_complex:
                if p in graph[protein]:
                    interactions.append(graph[protein][p])
            if interactions:
                avg_score = sum(interactions) / len(interactions)
                protein_scores.append((protein, avg_score))
    
        # Trier par score décroissant et prendre le top 20%
        protein_scores.sort(key=lambda x: x[1], reverse=True)
        top_proteins = [p[0] for p in protein_scores[:max(1, len(protein_scores)//3)]]
    
        for protein in top_proteins:
            new_complex = set(current_complex)
            new_complex.add(protein)
            chunk_neighbors.append(Solution(new_complex, graph))
    
        return chunk_neighbors

    def _generate_remove_protein_neighbors(self, solution: Solution) -> List[Solution]:
        """Supprime les protéines faiblement connectées du complexe"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        if len(current_complex) <= self.min_complex_size:
            return neighbors
        
        # Diviser les protéines du complexe en chunks
        proteins = list(current_complex)
        chunks = np.array_split(proteins, self.n_workers)
        
        with ThreadPool(self.n_workers) as pool:
            results = pool.starmap(
                self._process_remove_protein_chunk,
                [(current_complex, graph, chunk) for chunk in chunks if len(chunk) > 0]
            )
        
        for chunk_neighbors in results:
            neighbors.extend(chunk_neighbors)
            
        return neighbors

    def _process_remove_protein_chunk(self, current_complex: Set[str], graph: Dict, proteins: List[str]) -> List[Solution]:
        """Traite un chunk de protéines pour la suppression"""
        chunk_neighbors = []
        protein_scores = []
        
        for protein in proteins:
            interactions = []
            for p in current_complex:
                if p != protein and p in graph[protein]:
                    interactions.append(graph[protein][p])
            avg_score = sum(interactions)/len(interactions) if interactions else 0
            protein_scores.append((protein, avg_score))
    
        # Trier par score croissant et prendre le bottom 20%
        protein_scores.sort(key=lambda x: x[1])
        bottom_proteins = [p[0] for p in protein_scores[:max(1, len(protein_scores)//3)]]
    
        for protein in bottom_proteins:
            new_complex = set(current_complex)
            new_complex.remove(protein)
            chunk_neighbors.append(Solution(new_complex, graph))
    
        return chunk_neighbors

    def _generate_swap_neighbors(self, solution: Solution) -> List[Solution]:
        """Remplace des protéines faiblement connectées par des protéines externes fortement connectées"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        if len(current_complex) < 2:
            return neighbors
    
        # Identifier les protéines faiblement connectées dans le complexe (en parallèle)
        with ThreadPool(self.n_workers) as pool:
            weak_results = pool.starmap(
                self._find_weak_proteins,
                [(current_complex, graph, chunk) 
                 for chunk in np.array_split(list(current_complex), self.n_workers)])
        weak_proteins = [p for sublist in weak_results for p in sublist]
    
        # Identifier les protéines externes fortement connectées (en parallèle)
        external_proteins = list(set(graph.keys()) - set(current_complex))
        with ThreadPool(self.n_workers) as pool:
            strong_results = pool.starmap(
                self._find_strong_externals,
                [(current_complex, graph, chunk) 
                 for chunk in np.array_split(external_proteins, self.n_workers)])
        strong_externals = [p for sublist in strong_results for p in sublist]
    
        # Générer des voisins par échange
        swap_pairs = [(weak, strong) 
                     for weak in weak_proteins[:3] 
                     for strong in strong_externals[:3]]
        
        with ThreadPool(self.n_workers) as pool:
            neighbors = pool.starmap(
                self._create_swap_solution,
                [(current_complex, graph, weak, strong) 
                 for weak, strong in swap_pairs]
            )
    
        return [n for n in neighbors if n is not None]

    def _find_weak_proteins(self, current_complex: Set[str], graph: Dict, proteins: List[str]) -> List[str]:
        """Trouve les protéines faiblement connectées dans un chunk"""
        weak_proteins = []
        for protein in proteins:
            interactions = [graph[protein][p] for p in current_complex 
                          if p != protein and p in graph[protein]]
            if interactions:
                avg_score = sum(interactions)/len(interactions)
                if avg_score < 0.5:  # Seuil arbitraire
                    weak_proteins.append(protein)
        return weak_proteins

    def _find_strong_externals(self, current_complex: Set[str], graph: Dict, proteins: List[str]) -> List[str]:
        """Trouve les protéines externes fortement connectées dans un chunk"""
        strong_externals = []
        for protein in proteins:
            interactions = [graph[protein][p] for p in current_complex if p in graph[protein]]
            if interactions:
                avg_score = sum(interactions)/len(interactions)
                if avg_score > 0.7:  # Seuil arbitraire
                    strong_externals.append(protein)
        return strong_externals

    def _create_swap_solution(self, current_complex: Set[str], graph: Dict, 
                            weak: str, strong: str) -> Optional[Solution]:
        """Crée une nouvelle solution avec un échange de protéines"""
        new_complex = set(current_complex)
        new_complex.remove(weak)
        new_complex.add(strong)
        if self.min_complex_size <= len(new_complex) <= self.max_complex_size:
            return Solution(new_complex, graph)
        return None

    def _generate_merge_neighbors(self, solution: Solution) -> List[Solution]:
        """Fusionne avec des complexes voisins"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        # Trouver les protéines qui interagissent fortement avec le complexe actuel
        external_proteins = list(set(graph.keys()) - set(current_complex))
        chunks = np.array_split(external_proteins, self.n_workers)
        
        with ThreadPool(self.n_workers) as pool:
            results = pool.starmap(
                self._find_candidate_proteins,
                [(current_complex, graph, chunk) for chunk in chunks if len(chunk) > 0]
            )
        candidate_proteins = [p for sublist in results for p in sublist]
    
        # Créer des voisins en fusionnant avec ces protéines et leurs voisins
        merge_candidates = candidate_proteins[:5]  # Limiter le nombre de fusions
        with ThreadPool(self.n_workers) as pool:
            neighbors = pool.starmap(
                self._create_merge_solution,
                [(current_complex, graph, protein) for protein in merge_candidates]
            )
    
        return [n for n in neighbors if n is not None]

    def _find_candidate_proteins(self, current_complex: Set[str], graph: Dict, proteins: List[str]) -> List[str]:
        """Trouve les protéines candidates pour la fusion dans un chunk"""
        candidates = []
        for protein in proteins:
            interactions = [graph[protein][p] for p in current_complex if p in graph[protein]]
            if interactions and sum(interactions)/len(interactions) > 0.6:
                candidates.append(protein)
        return candidates

    def _create_merge_solution(self, current_complex: Set[str], graph: Dict, protein: str) -> Optional[Solution]:
        """Crée une nouvelle solution fusionnée"""
        new_complex = set(current_complex)
        new_complex.add(protein)
        for neighbor, weight in graph[protein].items():
            if weight > 0.7 and len(new_complex) < self.max_complex_size:
                new_complex.add(neighbor)
        
        if len(new_complex) > len(current_complex) and len(new_complex) <= self.max_complex_size:
            return Solution(new_complex, graph)
        return None

    def _generate_split_neighbors(self, solution: Solution) -> List[Solution]:
        """Divise le complexe en sous-complexes cohésifs"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph

        if len(current_complex) <= self.min_complex_size + 1:
            return neighbors

        try:
            # Construire le graphe induit par le complexe
            G_sub = nx.Graph()
            G_sub.add_nodes_from(current_complex)
            for u in current_complex:
                for v in graph[u]:
                    if v in current_complex:
                        G_sub.add_edge(u, v, weight=graph[u][v])

            # Traiter chaque composante connexe en parallèle
            components = list(nx.connected_components(G_sub))
            with ThreadPool(self.n_workers) as pool:
                results = pool.map(
                    self._process_component,
                    [(comp, graph, self.min_complex_size) for comp in components 
                     if len(comp) >= self.min_complex_size + 1]
                )
            
            for component_neighbors in results:
                neighbors.extend(component_neighbors)

        except Exception as e:
            print(f"  Erreur lors du clustering: {e}")
             
            
        return neighbors

    def _process_component(self, args: tuple) -> List[Solution]:
        """Traite une composante connexe pour la division"""
        component, graph, min_size = args
        neighbors = []
        proteins = list(component)
        n = len(proteins)
        similarity_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                if proteins[j] in graph[proteins[i]]:
                    weight = graph[proteins[i]][proteins[j]]
                    similarity_matrix[i][j] = weight
                    similarity_matrix[j][i] = weight

        n_clusters = min(3, len(proteins) // min_size)
        if n_clusters < 2:
            return neighbors

        clustering = SpectralClustering(
            n_clusters=n_clusters,
            affinity='precomputed',
            random_state=42
        ).fit(similarity_matrix)

        labels = clustering.labels_
        clusters = {}
        for i, label in enumerate(labels):
            clusters.setdefault(label, set()).add(proteins[i])

        for cluster in clusters.values():
            if len(cluster) >= min_size:
                neighbors.append(Solution(cluster, graph))

        return neighbors

    def _generate_density_based_neighbors(self, solution: Solution) -> List[Solution]:
        """Optimise le complexe en fonction de la densité"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        # Calculer la densité actuelle
        edges = 0
        total_possible = len(current_complex) * (len(current_complex) - 1) / 2
        for i, p1 in enumerate(current_complex):
            for p2 in list(current_complex)[i+1:]:
                if p2 in graph[p1]:
                    edges += 1
        current_density = edges / total_possible if total_possible > 0 else 0
    
        # Stratégie 1: Ajouter des protéines qui augmentent la densité (en parallèle)
        external_proteins = list(set(graph.keys()) - set(current_complex))
        add_args = [(current_complex, graph, edges, current_density, protein, self.max_complex_size)
                   for protein in external_proteins]
        
        with ThreadPool(self.n_workers) as pool:
            add_results = pool.starmap(self._process_add_for_density, add_args)
            neighbors.extend([sol for sol in add_results if sol is not None])
    
        # Stratégie 2: Supprimer des protéines qui augmentent la densité (en parallèle)
        if len(current_complex) > self.min_complex_size:
            remove_args = [(current_complex, graph, edges, current_density, protein, self.min_complex_size)
                         for protein in current_complex]
            
            with ThreadPool(self.n_workers) as pool:
                remove_results = pool.starmap(self._process_remove_for_density, remove_args)
                neighbors.extend([sol for sol in remove_results if sol is not None])
    
        return neighbors

    def _process_add_for_density(self, current_complex: Set[str], graph: Dict, 
                               edges: int, current_density: float, 
                               protein: str, max_size: int) -> Optional[Solution]:
        """Traite l'ajout d'une protéine pour améliorer la densité"""
        new_edges = sum(1 for p in current_complex if p in graph[protein])
        new_density = (edges + new_edges) / ((len(current_complex)+1)*len(current_complex)/2)
        if new_density > current_density and len(current_complex) < max_size:
            new_complex = set(current_complex)
            new_complex.add(protein)
            return Solution(new_complex, graph)
        return None

    def _process_remove_for_density(self, current_complex: Set[str], graph: Dict, 
                                  edges: int, current_density: float, 
                                  protein: str, min_size: int) -> Optional[Solution]:
        """Traite la suppression d'une protéine pour améliorer la densité"""
        edges_without = edges - sum(1 for p in current_complex if p != protein and p in graph[protein])
        new_size = len(current_complex) - 1
        denominator = (new_size)*(new_size-1)/2 if new_size > 1 else 1
        new_density = edges_without / denominator if denominator > 0 else 0
        if new_density > current_density and new_size >= min_size:
            new_complex = set(current_complex)
            new_complex.remove(protein)
            return Solution(new_complex, graph)
        return None

    @property
    def name(self) -> str:
        return "parallel_extended_complex_neighborhood"
    
class TabuMutator(Mutator):
    def __init__(self, 
                 graph: Dict[str, Dict[str, float]], 
                 tabu_tenure: int = 10, 
                 tabu_iterations: int = 10,
                 rate: float = 1.0,
                 min_complex_size: int = 3,
                 max_complex_size: int = 15,
                 threshold_ratio: int = 100):
        super().__init__(rate)
        self.graph = graph
        self.tabu_tenure = tabu_tenure
        self.tabu_iterations = tabu_iterations
        self.neighborhood = Neighborhood(min_complex_size, max_complex_size)
        self.threshold_ratio = threshold_ratio
        self.min_complex_size = min_complex_size
        self.max_complex_size = max_complex_size
    
    def __call__(self, chrom: Chromosome) -> Chromosome:
        if random.random() > self.rate:
            return chrom
        
        graph = chrom.graph
        complexes = chrom.genes
        fitnesses = chrom._fits
        indexes = np.argsort(fitnesses)
        complex = [complexes[i] for i in indexes]
        improved_complexes = complexes[:self.threshold_ratio] 
        complexes = [complex for complex in complexes if complex not in improved_complexes]
        new_complexes = []
        
        #self.logger.info(f"Nombre de complexes à améliorer: {len(improved_complexes)}")
         
        #self.logger.info("type de neighborhood", type(self.neighborhood))
         
        for c in improved_complexes:
            local_tabu = TabuSearch(
                initial_solution=Solution(c, chrom.graph),
                neighborhood_generator=self.neighborhood,
                get_move_hash=self._get_move_hash,
                tabu_tenure=self.tabu_tenure,
                update_history=None,
                apply_intensification=None,
                update_intensification_memory=None,
                apply_diversification=self.apply_diversification,
                update_diversification_memory=self.update_diversification_memory,
                max_iterations=self.tabu_iterations,
                diversification_necessity=5,
                diversification_frequency=20,
                intensification_threshold=20,
                patience=15
            )
            
            solution = local_tabu.run()
            new_complexes.append(solution.representation)
            #self.logger.info(f"Fitness ancien complexe: {FF(c, graph)} - Nouveau fitness: {FF(solution.representation, graph)}")
             
        complexes.extend(new_complexes)
        return Chromosome(self.graph, complexes)
        
    def update_diversification_memory(
        self, 
        solution: Solution, 
        memory: List[Dict[str, Any]],
        *,
        max_memory_size: int = 50,
        similarity_threshold: float = 0.6
    ) -> List[Dict[str, Any]]:
        """
        Met à jour la mémoire de diversification avec les caractéristiques des solutions explorées.
        
        Args:
            solution: Solution actuelle à ajouter à la mémoire
            memory: Mémoire actuelle de diversification
            max_memory_size: Nombre maximum de solutions à mémoriser
            similarity_threshold: Seuil de similarité pour éviter les doublons
            
        Returns:
            Mémoire mise à jour
        """
        # Caractéristiques de la solution actuelle
        current_features = {
            'proteins': frozenset(solution.representation),
            'structure': self._get_structure_hash(solution),
            'performance': solution.evaluate(),
            'size': len(solution.representation)
        }
        
        # Vérifier si la solution est similaire à celles en mémoire
        is_similar = False
        for entry in memory:
            similarity = self._calculate_solution_similarity(current_features, entry['features'])
            if similarity > similarity_threshold:
                entry['count'] += 1
                is_similar = True
                break
        
        # Ajouter une nouvelle entrée si la solution est suffisamment différente
        if not is_similar:
            memory.append({
                'features': current_features,
                'count': 1
            })
        
        # Trier et limiter la taille de la mémoire
        memory.sort(key=lambda x: -x['count'])
        if len(memory) > max_memory_size:
            memory = memory[:max_memory_size]
            
        return memory

    def apply_diversification(
        self, 
        current_solution: Solution,
        memory: List[Dict[str, Any]],
        *,
        intensity: float = 0.7,
        perturbation_size: int = 2
    ) -> Solution:
        """
        Applique une diversification basée sur la mémoire des solutions précédentes.
        
        Args:
            current_solution: Solution actuelle à diversifier
            memory: Mémoire de diversification
            intensity: Niveau d'intensité de la diversification [0,1]
            perturbation_size: Nombre de protéines à modifier
            
        Returns:
            Nouvelle solution diversifiée
        """
        if not memory or random.random() > intensity:
            return current_solution.copy()
        
        # Sélectionner une solution de référence sous-représentée
        underrepresented = [m for m in memory if m['count'] < np.mean([e['count'] for e in memory])]
        if not underrepresented:
            underrepresented = memory
            
        ref_entry = random.choice(underrepresented)
        ref_features = ref_entry['features']
        
        # Créer une nouvelle solution hybride
        current_set = set(current_solution.representation)
        ref_set = set(ref_features['proteins'])
        
        # Opérations de diversification
        if random.random() < 0.5:
            # Combinaison avec la solution de référence
            new_set = current_set.union(ref_set)
            if len(new_set) > self.max_complex_size:
                new_set = set(random.sample(list(new_set), self.max_complex_size))
        else:
            # Perturbation aléatoire guidée
            new_set = current_set.copy()
            to_remove = random.sample(list(current_set - ref_set), 
                                    min(perturbation_size, len(current_set - ref_set)))
            to_add = random.sample(list(ref_set - current_set), 
                                 min(perturbation_size, len(ref_set - current_set)))
            new_set.difference_update(to_remove)
            new_set.update(to_add)
        
        # Vérifier la taille minimale
        if len(new_set) < self.min_complex_size:
            needed = self.min_complex_size - len(new_set)
            candidates = set(self.graph.keys()) - new_set
            new_set.update(random.sample(list(candidates), needed))
        
        return Solution(new_set, self.graph)

    def _get_structure_hash(self, solution: Solution) -> int:
        """Calcule un hash basé sur la structure du complexe"""
        edges = []
        proteins = list(solution.representation)
        for i, p1 in enumerate(proteins):
            for p2 in proteins[i+1:]:
                if p2 in self.graph[p1]:
                    edges.append((p1, p2, self.graph[p1][p2]))
        return hash(frozenset(edges))

    def _calculate_solution_similarity(
        self, 
        features1: Dict[str, Any], 
        features2: Dict[str, Any]
    ) -> float:
        """Calcule la similarité entre deux solutions"""
        # Similarité basée sur le chevauchement des protéines
        overlap = len(features1['proteins'] & features2['proteins'])
        union = len(features1['proteins'] | features2['proteins']) or 1
        protein_sim = overlap / union
        
        # Similarité basée sur la performance (normalisée)
        perf_sim = 1 - abs(features1['performance'] - features2['performance'])/max(
            abs(features1['performance']), abs(features2['performance'])) if max(
            abs(features1['performance']), abs(features2['performance'])) != 0 else 0
        
        # Similarité basée sur la taille
        size_sim = 1 - abs(features1['size'] - features2['size'])/max(
            features1['size'], features2['size'])
        
        # Combinaison pondérée
        return 0.6 * protein_sim + 0.3 * perf_sim + 0.1 * size_sim

    def _get_move_hash(self, solution: Solution) -> int:
        """Génère un hash pour identifier le mouvement."""
        return hash(frozenset(solution.representation))
      
 
    
class ComplexCrossover(Crossover):
    """Opérateur de croisement intelligent pour les complexes"""
    def __init__(self, graph: Dict[str, Dict[str, float]]):
        self.graph = graph
        
    def __call__(self, p1: Chromosome, p2: Chromosome, log, rate:float = 1) -> Chromosome:
        if random.random() > rate:
            return p1 if p1._fitness > p2._fitness else p2
        # Combinaison des complexes des deux parents
        graph = p1.graph
        all_complexes = p1.genes + p2.genes
        all_fitness = np.concatenate((p1._fits, p2._fits))
        
        indexes = np.argsort(all_fitness)[::-1]
        all_complexes = [all_complexes[i] for i in indexes]
        all_fitness = [all_fitness[i] for i in indexes]
        
        min_len = min([len(p1.genes)*1.4, len(p2.genes)*1.4])
        max_len = max([len(p1.genes)*1, len(p2.genes)*1])
        log.info(f"max {max_len} min {min_len}")
         
        if min_len == max_len:
            min_len -= 1
        if min_len > max_len:
             min_len, max_len = max_len, min_len 
        target_size1 = np.random.randint(min_len, max_len)
        size1 = target_size1 + (len(all_complexes) - target_size1)//2

        target_size2 = np.random.randint(min_len, max_len)
        size2 = target_size2 + (len(all_complexes) - target_size2)//2

        #size = target_size
        
        all_fitness1 = [all_fitness[i] for i in indexes[:size1]]
        all_complexes1 = [all_complexes[i] for i in indexes[:size1]]
        all_fitness2 = [all_fitness[i] for i in indexes[:size2]]
        all_complexes2 = [all_complexes[i] for i in indexes[:size2]]

        # Taille cible basée sur une combinaison des tailles parentales
        child_genes1 = []
        child_genes2 = []
        while all_complexes1 and len(child_genes1) < target_size1:
            gene1 = np.random.randint(0, len(all_complexes1))
            gene2 = np.random.randint(0, len(all_complexes1))
            gene_child1 = all_complexes1[gene1] if all_fitness1[gene1] > all_fitness1[gene2] else all_complexes1[gene2]
            ind_child1 = gene1 if all_fitness1[gene1] > all_fitness1[gene2] else gene2
            all_complexes1.remove(gene_child1)
            all_fitness1.pop(ind_child1)
            for g in child_genes1:
                if overlap_score(gene_child1, g) > 0.2:
                    break
            child_genes1.append(set(gene_child1))
        
        while all_complexes2 and len(child_genes2) < target_size2:
            gene1 = np.random.randint(0, len(all_complexes2))
            gene2 = np.random.randint(0, len(all_complexes2))
            gene_child2 = all_complexes2[gene1] if all_fitness2[gene1] > all_fitness2[gene2] else all_complexes2[gene2]
            ind_child2 = gene1 if all_fitness2[gene1] > all_fitness2[gene2] else gene2
            all_complexes2.remove(gene_child2)
            all_fitness2.pop(ind_child2)
            
            for g in child_genes2:
                if overlap_score(gene_child2, g) > 0.2:
                    break
            child_genes2.append(set(gene_child2))
        c1 = Chromosome(p1.graph, child_genes1)
        c2 = Chromosome(p1.graph, child_genes2)
          
        return c1 if c1._fitness>c2._fitness else c2

class MyParallelTabuSearch(ParallelTabuSearch):
    def _run_1_control(self) -> List[Solution]:
        """1-control implementation (master-slave)."""
        workers = [
            self.TabuSearchWorker.remote(sol, config)
            for sol, config in zip(self.initial_solutions, self.worker_configs)
        ]
        results = ray.get([worker.run.remote() for worker in workers])
        for worker in workers:
            ray.kill(worker)
        return results
     
class TabuMutator_MasterSlave(Mutator):
    def __init__(self, 
                 graph: Dict[str, Dict[str, float]], 
                 tabu_tenure: int = 10, 
                 tabu_iterations: int = 10,
                 rate: float = 1.0,
                 min_complex_size: int = 3,
                 max_complex_size: int = 15,
                 threshold_ratio: int = 100,
                 nb_processors: int = 20):
        super().__init__(rate)
        self.graph = graph
        self.tabu_tenure = tabu_tenure
        self.tabu_iterations = tabu_iterations
        self.neighborhood = Neighborhood(min_complex_size, max_complex_size)
        self.threshold_ratio = threshold_ratio
        self.min_complex_size = min_complex_size
        self.max_complex_size = max_complex_size
        self.nb_processors = nb_processors
    def __call__(self, chrom: Chromosome, log, id) -> Chromosome:
        if random.random() > self.rate:
            return chrom
        
        graph = chrom.graph
        complexes = chrom.genes
        fitnesses = chrom._fits
        indexes = np.argsort(fitnesses)
        improved_complexes = [Solution(complexes[i], graph) for i in indexes[:self.threshold_ratio]]
        complexes = [complex for complex in complexes if complex not in improved_complexes]
        new_complexes = []
        
        log.info(f"Nombre de complexes à améliorer: {len(improved_complexes)}")
         
        solutions = []
        while improved_complexes:
            tabu_complexes = improved_complexes[:min(self.nb_processors, len(improved_complexes))]
            myParallelTabuSearch = MyParallelTabuSearch(initial_solutions=tabu_complexes, 
                                                    neighborhood_generator = self.neighborhood,
                                                    get_move_hash=self._get_move_hash,
                                                    tabu_tenure=self.tabu_tenure,
                                                    update_history=None,
                                                    apply_intensification=None,
                                                    update_intensification_memory=None,
                                                    apply_diversification=self.apply_diversification,
                                                    update_diversification_memory=self.update_diversification_memory,
                                                    max_iterations=self.tabu_iterations,
                                                    diversification_necessity=5,
                                                    diversification_frequency=20,
                                                    intensification_threshold=20,
                                                    control_cardinality = "1-control",
                                                    communication_type = "rigid",
                                                    patience=15)
            log.info(f"Lancement de la recherche Tabu parallèle pour {len(tabu_complexes)} complexes")
            solutions.extend(myParallelTabuSearch._run_1_control())
            del myParallelTabuSearch
            improved_complexes = [c for c in improved_complexes if c not in tabu_complexes]
        
        best_solutions = [s.representation for s in solutions]
        
        complexes.extend(best_solutions)
        return Chromosome(graph, complexes)
       
    def update_diversification_memory(
        self, 
        solution: Solution, 
        memory: List[Dict[str, Any]],
        *,
        max_memory_size: int = 50,
        similarity_threshold: float = 0.6
    ) -> List[Dict[str, Any]]:
        """
        Met à jour la mémoire de diversification avec les caractéristiques des solutions explorées.
        
        Args:
            solution: Solution actuelle à ajouter à la mémoire
            memory: Mémoire actuelle de diversification
            max_memory_size: Nombre maximum de solutions à mémoriser
            similarity_threshold: Seuil de similarité pour éviter les doublons
            
        Returns:
            Mémoire mise à jour
        """
        # Caractéristiques de la solution actuelle
        current_features = {
            'proteins': frozenset(solution.representation),
            'structure': self._get_structure_hash(solution),
            'performance': solution.evaluate(),
            'size': len(solution.representation)
        }
        
        # Vérifier si la solution est similaire à celles en mémoire
        is_similar = False
        for entry in memory:
            similarity = self._calculate_solution_similarity(current_features, entry['features'])
            if similarity > similarity_threshold:
                entry['count'] += 1
                is_similar = True
                break
        
        # Ajouter une nouvelle entrée si la solution est suffisamment différente
        if not is_similar:
            memory.append({
                'features': current_features,
                'count': 1
            })
        
        # Trier et limiter la taille de la mémoire
        memory.sort(key=lambda x: -x['count'])
        if len(memory) > max_memory_size:
            memory = memory[:max_memory_size]
            
        return memory

    def apply_diversification(
        self, 
        current_solution: Solution,
        memory: List[Dict[str, Any]],
        *,
        intensity: float = 0.7,
        perturbation_size: int = 2
    ) -> Solution:
        """
        Applique une diversification basée sur la mémoire des solutions précédentes.
        
        Args:
            current_solution: Solution actuelle à diversifier
            memory: Mémoire de diversification
            intensity: Niveau d'intensité de la diversification [0,1]
            perturbation_size: Nombre de protéines à modifier
            
        Returns:
            Nouvelle solution diversifiée
        """
        if not memory or random.random() > intensity:
            return current_solution.copy()
        
        # Sélectionner une solution de référence sous-représentée
        underrepresented = [m for m in memory if m['count'] < np.mean([e['count'] for e in memory])]
        if not underrepresented:
            underrepresented = memory
            
        ref_entry = random.choice(underrepresented)
        ref_features = ref_entry['features']
        
        # Créer une nouvelle solution hybride
        current_set = set(current_solution.representation)
        ref_set = set(ref_features['proteins'])
        
        # Opérations de diversification
        if random.random() < 0.5:
            # Combinaison avec la solution de référence
            new_set = current_set.union(ref_set)
            if len(new_set) > self.max_complex_size:
                new_set = set(random.sample(list(new_set), self.max_complex_size))
        else:
            # Perturbation aléatoire guidée
            new_set = current_set.copy()
            to_remove = random.sample(list(current_set - ref_set), 
                                    min(perturbation_size, len(current_set - ref_set)))
            to_add = random.sample(list(ref_set - current_set), 
                                 min(perturbation_size, len(ref_set - current_set)))
            new_set.difference_update(to_remove)
            new_set.update(to_add)
        
        # Vérifier la taille minimale
        if len(new_set) < self.min_complex_size:
            needed = self.min_complex_size - len(new_set)
            candidates = set(self.graph.keys()) - new_set
            new_set.update(random.sample(list(candidates), needed))
        
        return Solution(new_set, self.graph)

    def _get_structure_hash(self, solution: Solution) -> int:
        """Calcule un hash basé sur la structure du complexe"""
        edges = []
        proteins = list(solution.representation)
        for i, p1 in enumerate(proteins):
            for p2 in proteins[i+1:]:
                if p2 in self.graph[p1]:
                    edges.append((p1, p2, self.graph[p1][p2]))
        return hash(frozenset(edges))

    def _calculate_solution_similarity(
        self, 
        features1: Dict[str, Any], 
        features2: Dict[str, Any]
    ) -> float:
        """Calcule la similarité entre deux solutions"""
        # Similarité basée sur le chevauchement des protéines
        overlap = len(features1['proteins'] & features2['proteins'])
        union = len(features1['proteins'] | features2['proteins']) or 1
        protein_sim = overlap / union
        
        # Similarité basée sur la performance (normalisée)
        perf_sim = 1 - abs(features1['performance'] - features2['performance'])/max(
            abs(features1['performance']), abs(features2['performance'])) if max(
            abs(features1['performance']), abs(features2['performance'])) != 0 else 0
        
        # Similarité basée sur la taille
        size_sim = 1 - abs(features1['size'] - features2['size'])/max(
            features1['size'], features2['size'])
        
        # Combinaison pondérée
        return 0.6 * protein_sim + 0.3 * perf_sim + 0.1 * size_sim

    def _get_move_hash(self, solution: Solution) -> int:
        """Génère un hash pour identifier le mouvement."""
        return hash(frozenset(solution.representation))
      
    
    
class AlgorithmGA(GeneticAlgorithm):
    def __init__(
        self, 
        references: List[List[str]],
        ppi_name:str,
        graph: Dict[str, Dict[str, float]],
        population: MyPopulation,
        crossover: Crossover,
        mutator: Mutator,
        selection: Callable,
        update_population: Callable,
        max_generations: int,
        nb_parents : int=2,
        k_tournament: int = 3,
        update_type: str = "replace",
         **kwargs
    ):
        super().__init__(population, crossover, mutator, nb_parents, update_population, update_type, selection, max_generations, k_tournament)
        self.graph = graph
        self.references = references
        self.id = random.randint(0, 1000)
        setup_worker_logger(type="island", name=ppi_name)
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Island {self.id} initialisée.")
        best = self.population.best()
        masterSlaveMetrics = MasterSlaveMetrics(best.genes, self.references)
        metrics = masterSlaveMetrics.compute_all_metrics()
        self.iter = 0
        self.name = ppi_name
        
        self._statistiques = {
            'f_mesure': [metrics['F-mesure']],
            'PPV': [metrics['PPV']],
            'Recall(Sn)': [metrics["Recall (Sn)"]],
            'accuracy': [metrics['Accuracy']],
            'jaccard_index': [metrics['Jaccard']],
            'mmr': [metrics['MMR']],
            'cr': [metrics['Covered Rate']],
            'total_score': [metrics['Score Total']],
            'fitness': [best._fitness], 
            'avg': [population.average_fitness()]
        }
        
    def _evolve(self) -> MyPopulation:
        """Perform one generation of evolution."""
        self.iter +=1
        self.logger.info(f"  selection commence")
         
        parent1 = self.population.tournament_selection(2)
        parent2 = self.population.tournament_selection(2)
        
        worst_parent = parent1 if parent1._fitness < parent2._fitness else parent2
        self.logger.info(f"  parent1._fitness {parent1._fitness} parent2._fitness {parent2._fitness}")
         
        
        self.logger.info(f"  crossover commence")
         
        child = self.crossover(parent1, parent2, self.logger)
        self.logger.info(f" child fitness {child._fitness}")
         
        
        self.logger.info(f" crossover termine\n mutation commence")
         
        mutated_child = self.mutator(child, self.logger, self.id)
        self.logger.info(f"  mutated fitness {mutated_child._fitness}")
          
        if mutated_child.fitness > worst_parent._fitness:
            self.population.individuals.append(mutated_child)
            self.logger.info(f"  worst parent fitness {worst_parent._fitness} < child fitness {mutated_child._fitness}")
             
            if worst_parent in self.population.individuals:
                self.logger.info(f"  removing worst parent with fitness {worst_parent._fitness}") 
                 
                self.population.individuals.remove(worst_parent)
                
            self.logger.info(f"generation {self.iter} termine")
             
        self.iter+=1       
        return self.population
        
    def _update_history(self):
        best_chromosome = self.population.best()  # Get the best chromosome
        best_solution = best_chromosome.genes      # Get the genes (complexes)
        masterSlaveMetrics = MasterSlaveMetrics(best_solution, self.references)
        metrics = masterSlaveMetrics.compute_all_metrics()

        # Update history with best chromosome's fitness, not genes' fitness
        self.history["best"].append(best_chromosome._fitness)
        self.history["avg"].append(self.population.average_fitness())
    
        self._statistiques['accuracy'].append(metrics['Accuracy'])
        self._statistiques['PPV'].append(metrics['PPV'])
        self._statistiques['Recall(Sn)'].append(metrics['Recall (Sn)'])
        self._statistiques['f_mesure'].append(metrics['F-mesure'])
        self._statistiques['jaccard_index'].append(metrics['Jaccard'])
        self._statistiques['mmr'].append(metrics['MMR'])
        self._statistiques['cr'].append(metrics['Covered Rate'])
        self._statistiques['total_score'].append(metrics['Score Total'])
        self._statistiques['fitness'].append(best_chromosome._fitness)
        self._statistiques["avg"].append(self.population.average_fitness())
        if self.iter % 10 == 0:
            df = pd.DataFrame(self._statistiques)
            self.logger.info(f"  mise à jour des statistiques")
         
            df.to_csv(f"kefsi_mekhazni_workspace/results/{self.name}/metrics_{self.name}_{self.id}.csv", index=False)
        
            save_complexes(best_solution, f"kefsi_mekhazni_workspace/results/{self.name}/complexes_{self.name}_{self.id}.txt")
           