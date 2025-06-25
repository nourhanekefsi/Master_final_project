import numpy as np
from collections import defaultdict

# Fonctions pour évaluer un seul cluster
def cohesiveness(cluster, graph):
    """
    Calcule la cohésion d'un cluster.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI sous forme de dictionnaire {nœud: {voisins: poids}}
    :return: Score de cohésion
    """
    W_in = 0
    W_out = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor in cluster:
                W_in += weight
            else:
                W_out += weight
    return W_in / (W_in + W_out) if (W_in + W_out) > 0 else 0

def density(cluster, graph):
    """
    Calcule la densité d'un cluster.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score de densité
    """
    W_in = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor in cluster:
                W_in += weight
    n = len(cluster)
    return (2 * W_in) / (n * (n - 1)) if n > 1 else 0

def AIEW(cluster, graph):
    """
    Calcule le poids moyen des arêtes internes (Average Inner Edge Weight).
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score AIEW
    """
    W_in = 0
    E_C = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor in cluster:
                W_in += weight
                E_C += 1
    return W_in / E_C if E_C > 0 else 0

def ABEW(cluster, graph):
    """
    Calcule le poids moyen des arêtes frontalières (Average Border Edge Weight).
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score ABEW
    """
    W_out = 0
    BE_C = 0
    for node in cluster:
        for neighbor, weight in graph.get(node, {}).items():
            if neighbor not in cluster:
                W_out += weight
                BE_C += 1
    return W_out / BE_C if BE_C > 0 else 0

def AWM(cluster, graph):
    """
    Calcule la modularité pondérée moyenne (Average Weighted Modularity).
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score AWM
    """
    aiew = AIEW(cluster, graph)
    abew = ABEW(cluster, graph)
    return aiew / (aiew + abew) if (aiew + abew) > 0 else 0

def FF(cluster, graph):
    """
    Fonction de fitness pour un seul cluster.
    :param cluster: Ensemble des nœuds du cluster
    :param graph: Graphe PPI
    :return: Score FF combiné
    """
    return (
        cohesiveness(cluster, graph) +
        density(cluster, graph) +
        AIEW(cluster, graph) -
        ABEW(cluster, graph) +
        AWM(cluster, graph)
    )

# Fonction pour évaluer un individu (ensemble de clusters)
def FS_fitness(individual, graph):
    """
    Calcule le fitness total d'un individu (ensemble de clusters).
    :param individual: Liste de clusters
    :param graph: Graphe PPI
    :return: Score FS_fitness
    """
    fitness_values = [FF(cluster, graph) for cluster in individual]
    return np.sum(fitness_values) if fitness_values else 0


#========================================================================================================================#
                                        # Métriques de performance
#=========================================================================================================================
import numpy as np
from math import sqrt

def overlap_score(pred, real):
    return len(set(pred) & set(real)) / sqrt(len(pred) * len(real)) if pred and real else 0

def jaccard_index(pred, real):
    return len(set(pred) & set(real)) / len(set(pred) | set(real)) if set(pred) | set(real) else 0

def compute_metrics(predicted_complexes, real_complexes, threshold=0.2):
    m = len(predicted_complexes)
    n = len(real_complexes)

    # Intersection max for PPV
    ppv_numerator = sum(max(len(set(pred) & set(real)) for real in real_complexes) for pred in predicted_complexes)
    ppv_denominator = sum(len(pred) for pred in predicted_complexes)
    ppv = ppv_numerator / ppv_denominator if ppv_denominator > 0 else 0

    # Intersection max for Sn
    sn_numerator = sum(max(len(set(pred) & set(real)) for pred in predicted_complexes) for real in real_complexes)
    sn_denominator = sum(len(real) for real in real_complexes)
    sn = sn_numerator / sn_denominator if sn_denominator > 0 else 0

    # F-measure
    f_measure = (2 * ppv * sn) / (ppv + sn) if (ppv + sn) > 0 else 0

    # Accuracy (as sqrt(Sn * PPV))
    accuracy = sqrt(sn * ppv)

    # MMR
    mmr = np.mean([max(overlap_score(pred, real) for real in real_complexes) for pred in predicted_complexes]) if predicted_complexes else 0

    # Jaccard
    jaccard = np.mean([max(jaccard_index(pred, real) for real in real_complexes) for pred in predicted_complexes]) if predicted_complexes else 0

    # Covered Rate
    matched_reals = sum(
        any(overlap_score(pred, real) >= threshold for pred in predicted_complexes)
        for real in real_complexes
    )
    covered_rate = matched_reals / n if n > 0 else 0

    # Score total
    score_total = f_measure + covered_rate + mmr + jaccard + accuracy

    return {
        "PPV": ppv,
        "Recall (Sn)": sn,
        "F-mesure": f_measure,
        "Accuracy": accuracy,
        "MMR": mmr,
        "Jaccard": jaccard,
        "Covered Rate": covered_rate,
        "Score Total": score_total
    }



from typing import Any, List, Dict, Set, Tuple, Callable, Optional
from collections import deque
import numpy as np
import pandas as pd
import os
import random
from HPC_GA.common.population import Population
from HPC_GA.common.chromosome import Chromosome
from HPC_GA.core.genetic_algorithm import GeneticAlgorithm
from HPC_GA.core.operators import Crossover, Mutator
from HPC_Tabu.sequential.tabu_search import TabuSearch
from HPC_Tabu.common.neighborhood import NeighborhoodGenerator
from HPC_Tabu.common.solution import Solution

class Chromosome_complexes(Chromosome):
    """Chromosome représentant un ensemble de complexes protéiques avec standardisation des IDs"""
    def __init__(self, graph: Dict[str, Dict[str, float]], genes: Optional[List[Set[str]]] = None):
        
        super().__init__(genes)
        self.graph = graph
        self._fits = np.array([FF(gene, self.graph) for gene in self.genes])
        self._fitness = np.sum(self._fits)
        
    def evaluate(self) -> float:
        return FS_fitness(self.genes, self.graph)

    

class Solution_complexes(Solution):
    """Solution pour la recherche tabou avec gestion standardisée des protéines"""
    def __init__(self, representation: Set[str], graph: Dict[str, Dict[str, float]]):
        # Standardisation des IDs de protéines
        representation = set(sorted(set(representation)))
        super().__init__(representation)
        self.graph = graph
        self._value = FF(representation, graph)
        
    def _evaluate(self) -> float:
        return FF(self.representation, self.graph)
    
    def copy(self) -> 'Solution_complexes':
        return Solution_complexes(sorted(set([p for p in self.representation])), self.graph)
    
    
    def __hash__(self):
        # Crée un hash stable basé sur la représentation
        return hash(frozenset(self.representation))
    
    def __eq__(self, other):
        if not isinstance(other, Solution_complexes):
            return False
        return sorted(self.representation) == sorted(other.representation)

class ComplexNeighborhood(NeighborhoodGenerator):
    """Générateur de voisinage étendu avec des opérations plus intelligentes"""
    
    def __init__(self, min_complex_size: int = 3, max_complex_size: int = 15):
        self.min_complex_size = min_complex_size
        self.max_complex_size = max_complex_size
    
    def generate(self, solution: Solution_complexes) -> List[Solution_complexes]:
        print(f"Taille du complexe initial: {len(solution.representation)}")
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        # 1. Ajout de protéines hautement connectées
        neighbors.extend(self._generate_add_protein_neighbors(solution))
    
        # 2. Suppression de protéines peu connectées
        neighbors.extend(self._generate_remove_protein_neighbors(solution))
    
        # 3. Remplacement intelligent de protéines
        neighbors.extend(self._generate_swap_neighbors(solution))
    
        # 4. Fusion avec des complexes voisins
        neighbors.extend(self._generate_merge_neighbors(solution))
    
        # 5. Division du complexe en sous-complexes cohésifs
        neighbors.extend(self._generate_split_neighbors(solution))
    
        # 6. Optimisation locale basée sur la densité
        neighbors.extend(self._generate_density_based_neighbors(solution))
    
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
    
        print(f"Généré {len(unique_neighbors)} voisins uniques")
        return unique_neighbors

    def _generate_add_protein_neighbors(self, solution: Solution_complexes) -> List[Solution_complexes]:
        """Ajoute des protéines fortement connectées au complexe"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        # Calculer la connectivité moyenne avec les protéines externes
        external_proteins = set(graph.keys()) - set(current_complex)
        protein_scores = []
    
        for protein in external_proteins:
            # Calculer le score moyen d'interaction avec le complexe
            interactions = []
            for p in current_complex:
                if p in graph[protein]:
                    interactions.append(graph[protein][p])
            if interactions:
                avg_score = sum(interactions) / len(interactions)
                protein_scores.append((protein, avg_score))
    
        # Trier par score décroissant et prendre le top 20%
        protein_scores.sort(key=lambda x: x[1], reverse=True)
        top_proteins = [p[0] for p in protein_scores[:max(1, len(protein_scores)//5)]]
    
        for protein in top_proteins:
            new_complex = set(current_complex)
            new_complex.add(protein)
            neighbors.append(Solution_complexes(new_complex, graph))
    
        return neighbors

    def _generate_remove_protein_neighbors(self, solution: Solution_complexes) -> List[Solution_complexes]:
        """Supprime les protéines faiblement connectées du complexe"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        if len(current_complex) <= self.min_complex_size:
            return neighbors
    
        # Calculer le score de connectivité pour chaque protéine
        protein_scores = []
        for protein in current_complex:
            interactions = []
            for p in current_complex:
                if p != protein and p in graph[protein]:
                    interactions.append(graph[protein][p])
            avg_score = sum(interactions)/len(interactions) if interactions else 0
            protein_scores.append((protein, avg_score))
    
        # Trier par score croissant et prendre le bottom 20%
        protein_scores.sort(key=lambda x: x[1])
        bottom_proteins = [p[0] for p in protein_scores[:max(1, len(protein_scores)//5)]]
    
        for protein in bottom_proteins:
            new_complex = set(current_complex)
            new_complex.remove(protein)
            neighbors.append(Solution_complexes(new_complex, graph))
    
        return neighbors

    def _generate_swap_neighbors(self, solution: Solution_complexes) -> List[Solution_complexes]:
        """Remplace des protéines faiblement connectées par des protéines externes fortement connectées"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        if len(current_complex) < 2:
            return neighbors
    
        # Identifier les protéines faiblement connectées dans le complexe
        weak_proteins = []
        for protein in current_complex:
            interactions = [graph[protein][p] for p in current_complex if p != protein and p in graph[protein]]
            if interactions:
                avg_score = sum(interactions)/len(interactions)
                if avg_score < 0.5:  # Seuil arbitraire
                    weak_proteins.append(protein)
    
        # Identifier les protéines externes fortement connectées
        external_proteins = set(graph.keys()) - set(current_complex)
        strong_externals = []
        for protein in external_proteins:
            interactions = [graph[protein][p] for p in current_complex if p in graph[protein]]
            if interactions:
                avg_score = sum(interactions)/len(interactions)
                if avg_score > 0.7:  # Seuil arbitraire
                   strong_externals.append(protein)
    
        # Générer des voisins par échange
        for weak in weak_proteins[:3]:  # Limiter le nombre d'échanges
            for strong in strong_externals[:3]:
                new_complex = set(current_complex)
                new_complex.remove(weak)
                new_complex.add(strong)
                neighbors.append(Solution_complexes(new_complex, graph))
    
        return neighbors

    def _generate_merge_neighbors(self, solution: Solution_complexes) -> List[Solution_complexes]:
        """Fusionne avec des complexes voisins"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        # Trouver les protéines qui interagissent fortement avec le complexe actuel
        external_proteins = set(graph.keys()) - set(current_complex)
        candidate_proteins = []
    
        for protein in external_proteins:
            interactions = [graph[protein][p] for p in current_complex if p in graph[protein]]
            if interactions and sum(interactions)/len(interactions) > 0.6:
                candidate_proteins.append(protein)
    
        # Créer des voisins en fusionnant avec ces protéines et leurs voisins
        for protein in candidate_proteins[:5]:  # Limiter le nombre de fusions
            # Ajouter la protéine et ses voisins fortement connectés
            new_complex = set(current_complex)
            new_complex.add(protein)
            for neighbor, weight in graph[protein].items():
                if weight > 0.7 and len(new_complex) < self.max_complex_size:
                    new_complex.add(neighbor)
        
            if len(new_complex) > len(current_complex):
                neighbors.append(Solution_complexes(new_complex, graph))
    
        return neighbors

    def _generate_split_neighbors(self, solution: Solution_complexes) -> List[Solution_complexes]:
        """Divise le complexe en sous-complexes cohésifs"""
        neighbors = []
        current_complex = solution.representation
        graph = solution.graph
    
        if len(current_complex) <= self.min_complex_size + 1:
            return neighbors
    
        # Utiliser un algorithme de clustering pour diviser le complexe
        try:
            # Convertir en matrice de similarité
            proteins = list(current_complex)
            n = len(proteins)
            similarity_matrix = np.zeros((n, n))
        
            for i in range(n):
                for j in range(i+1, n):
                    if proteins[j] in graph[proteins[i]]:
                        similarity_matrix[i][j] = graph[proteins[i]][proteins[j]]
                        similarity_matrix[j][i] = similarity_matrix[i][j]
        
            # Utiliser un clustering spectral
            from sklearn.cluster import SpectralClustering
            n_clusters = min(3, len(current_complex)//self.min_complex_size)
            if n_clusters < 2:
                return neighbors
            
            clustering = SpectralClustering(n_clusters=n_clusters, 
                                      affinity='precomputed',
                                      random_state=42).fit(similarity_matrix)
        
            labels = clustering.labels_
            clusters = {}
            for i, label in enumerate(labels):
                if label not in clusters:
                    clusters[label] = set()
                clusters[label].add(proteins[i])
        
            # Créer des solutions pour chaque cluster assez grand
            for cluster in clusters.values():
                if len(cluster) >= self.min_complex_size:
                    neighbors.append(Solution_complexes(cluster, graph))
    
        except Exception as e:
            print(f"Erreur lors du clustering: {e}")
    
        return neighbors

    def _generate_density_based_neighbors(self, solution: Solution_complexes) -> List[Solution_complexes]:
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
    
        # Stratégie 1: Ajouter des protéines qui augmentent la densité
        external_proteins = set(graph.keys()) - set(current_complex)
        for protein in external_proteins:
            new_edges = sum(1 for p in current_complex if p in graph[protein])
            new_density = (edges + new_edges) / ((len(current_complex)+1)*len(current_complex)/2)
            if new_density > current_density and len(current_complex) < self.max_complex_size:
                new_complex = set(current_complex)
                new_complex.add(protein)
                neighbors.append(Solution_complexes(new_complex, graph))
    
        # Stratégie 2: Supprimer des protéines qui augmentent la densité
        if len(current_complex) > self.min_complex_size:
            for protein in current_complex:
                edges_without = edges - sum(1 for p in current_complex if p != protein and p in graph[protein])
                new_density = edges_without / ((len(current_complex)-1)*(len(current_complex)-2)/2) if len(current_complex) > 2 else 0
                if new_density > current_density:
                    new_complex = set(current_complex)
                    new_complex.remove(protein)
                    neighbors.append(Solution_complexes(new_complex, graph))
    
        return neighbors
          
    @property
    def name(self) -> str:
        return "extended_complex_neighborhood"

class TabuMutator(Mutator):
    def __init__(self, 
                 graph: Dict[str, Dict[str, float]], 
                 tabu_tenure: int = 10, 
                 tabu_iterations: int = 10,
                 rate: float = 1.0,
                 min_complex_size: int = 3,
                 max_complex_size: int = 30,
                 threshold_ratio: int = 100):
        super().__init__(rate)
        self.graph = graph
        self.tabu_tenure = tabu_tenure
        self.tabu_iterations = tabu_iterations
        self.neighborhood = ComplexNeighborhood(min_complex_size, max_complex_size)
        self.threshold_ratio = threshold_ratio
        self.min_complex_size = min_complex_size
        self.max_complex_size = max_complex_size
    
    def __call__(self, chrom: Chromosome_complexes) -> Chromosome_complexes:
        if random.random() > self.rate:
            return chrom
        
        graph = chrom.graph
        complexes = chrom.genes
        fitnesses = chrom._fits
        indexes = np.argsort(fitnesses)
        complexes = [complexes[i] for i in indexes]
        improved_complexes = complexes[:self.threshold_ratio] 
        complexes = [complex for complex in complexes if complex not in improved_complexes]
        new_complexes = []
        
        print(f"Nombre de complexes à améliorer: {len(improved_complexes)}")
        
        for c in improved_complexes:
            local_tabu = TabuSearch(
                initial_solution=Solution_complexes(c, chrom.graph),
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
            print(f"Fitness ancien complexe: {FF(c, graph)} - Nouveau fitness: {FF(solution.representation, graph)}")
        
        complexes.extend(new_complexes)
        return Chromosome_complexes(self.graph, complexes)
        
    def update_diversification_memory(
        self, 
        solution: Solution_complexes, 
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
        current_solution: Solution_complexes,
        memory: List[Dict[str, Any]],
        *,
        intensity: float = 0.7,
        perturbation_size: int = 2
    ) -> Solution_complexes:
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
        
        return Solution_complexes(new_set, self.graph)

    def _get_structure_hash(self, solution: Solution_complexes) -> int:
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
        
    def __call__(self, p1: Chromosome_complexes, p2: Chromosome_complexes, rate:float = 1) -> Chromosome_complexes:
        if random.random() > rate:
            return p1 if p1._fitness > p2._fitness else p2
        # Combinaison des complexes des deux parents
        graph = p1.graph
        all_complexes = p1.genes + p2.genes
        all_fitness = np.concatenate((p1._fits, p2._fits))
        
        indexes = np.argsort(all_fitness)[::-1]
        all_complexes = [all_complexes[i] for i in indexes]
        all_fitness = [all_fitness[i] for i in indexes]
        
        min_len = min([len(p1.genes)*1.2, len(p2.genes)*1.2])
        max_len = max([len(p1.genes)*1.6, len(p2.genes)*1.6])
        print(f"max {max_len} min {min_len}")
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
        c1 = Chromosome_complexes(p1.graph, child_genes1)
        c2 = Chromosome_complexes(p1.graph, child_genes2)
          
        return c1 if c1._fitness>c2._fitness else c2

class AlgorithmGA(GeneticAlgorithm):
    def __init__(
        self, 
        references: List[List[str]],
        graph: Dict[str, Dict[str, float]],
        population: Population,
        crossover: Crossover,
        mutator: Mutator,
        selection: Callable,
        max_generations: int
    ):
        super().__init__(population, crossover, mutator, selection, max_generations)
        self.graph = graph
        self.references = references
        best = self.population.best()
        metrics = compute_metrics(best.genes, self.references)
        
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
        
    def _evolve(self) -> Population:
        """Perform one generation of evolution."""
        print("selection commence")
        parent1 = self.selection(2)
        parent2 = self.selection(2)
        worst_parent = parent1 if parent1._fitness < parent2._fitness else parent2
        print(f"parent1._fitness {parent1._fitness} parent2._fitness {parent2._fitness}")
        
        print("crossover commence")
        child = self.crossover(parent1, parent2)
        print(f"child fitness {child._fitness}")
        
        print("crossover termine\n mutation commence")
        mutated_child = self.mutator(child)
        print(f"mutated fitness {mutated_child._fitness}")
        print(f"lenth of mutated child :{len(mutated_child.genes)}")
        
         
        fits = [ind._fitness for ind in self.population.individuals]
        print("fits avant = ", sorted(fits))
        if mutated_child.fitness > worst_parent._fitness:
            self.population.individuals.remove(worst_parent)
            self.population.individuals.append(mutated_child)
            print("new population", len(self.population.individuals))
            fits = [ind._fitness for ind in self.population.individuals]
            print("fits apres = ", sorted(fits))
        return self.population
        
    def _update_history(self):
        best_chromosome = self.population.best()  # Get the best chromosome
        best_solution = best_chromosome.genes      # Get the genes (complexes)
        metrics = compute_metrics(best_solution, self.references)

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
        
        df = pd.DataFrame(self._statistiques)
        df.to_csv(r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\results\mainMethod\OptClust-H\metrics\metrics_STRING_levure.csv", index=False)

        save_complexes(best_solution, r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\results\mainMethod\OptClust-H\complexes\complexes_STRING_levure.txt")
        solution_id = 0
        complexe_id = 0
        population_safe = {'solution_id': [], 'complexe_id': [], 'complexe': []}

        for ind in self.population.individuals:
            for c in ind.genes:
               population_safe['solution_id'].append(solution_id)
               population_safe['complexe_id'].append(complexe_id)
               population_safe['complexe'].append(c)
               complexe_id += 1
            solution_id += 1

        population_pd = pd.DataFrame(population_safe)
        population_pd.to_csv("/kaggle/working/population_STRING_levure.txt", index=False)
     
def load_ppi(ppi_file: str) -> Dict[str, Dict[str, float]]:
    """Load protein-protein interaction network"""
    graph_df = pd.read_csv(ppi_file, sep="\t", header=0, names=["node1", "node2", "weight"])
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
        print(f"Erreur lors du chargement des complexes: {e}")
        return []

def random_initialization(reference_complexes: List[List[str]], 
                         graph: Dict[str, Dict[str, float]], 
                         population_size: int = 26,
                         max_complexes: int = 198, fact : int = 6, min_size=3, max_size=25) -> Population:
    """
    Génère une population avec des complexes protéiques aléatoires et connectés
    Amélioré pour créer des complexes plus réalistes
    """
    print("Debut de generation de population initiale")
    protein_nodes = list(graph.keys())
    individuals = []
    nb_solutions = fact*population_size
    
    def create_connected_complex(min_size, max_size) -> Set[str]:
        """Crée un seul complexe protéique connecté"""
        if not protein_nodes:
            return set()
            
        # Choisir une protéine de départ avec des connexions
        start_node = random.choice(protein_nodes)
        while len(graph.get(start_node, {})) == 0:
            start_node = random.choice(protein_nodes)
        
        complex_nodes = set()
        queue = deque([start_node])
        target_size = random.randint(min_size, max_size)
        
        # Breadth-first search pour faire croître le complexe
        while queue and len(complex_nodes) < target_size:
            current = queue.popleft()
            if current in complex_nodes:
                continue
                
            complex_nodes.add(current)
            
            # Trier les voisins par force d'interaction
            neighbors = sorted(
                [(n, w) for n, w in graph.get(current, {}).items() 
                 if n not in complex_nodes],
                key=lambda x: -x[1]  # Plus fortes interactions en premier
            )
            
            # Ajouter les voisins les plus connectés
            for neighbor, _ in neighbors[:3]:  # Limiter à 3 meilleurs voisins
                if len(complex_nodes) >= target_size:
                    break
                queue.append(neighbor)
        
        return complex_nodes
    
    while len(individuals) < nb_solutions:
        try:
            # Créer entre 50% et 150% du nombre moyen de complexes de référence
            ref_mean = max_complexes
            num_complexes = random.randint(int(0.75 * ref_mean), int(1.25 * ref_mean))
            #num_complexes = max_complexes
            complexes = []
            used_proteins = set()
            
            for _ in range(num_complexes):
                # Essayer de créer un complexe connecté
                new_complex = create_connected_complex(min_size, max_size)
                if new_complex:
                    # Vérifier qu'on n'a pas trop de chevauchement
                    if len(new_complex & used_proteins) / len(new_complex) < 0.2:
                        complexes.append(new_complex)
                        used_proteins.update(new_complex)
            
            if complexes:
                new_individual = Chromosome_complexes(graph, complexes)
                individuals.append(new_individual)
        except Exception as e:
            print(f"Error generating individual: {e}")
            continue
    fitnesses = [individual._fitness for individual in individuals]
    indices = np.argsort(fitnesses)[::-1]
    individuals = [individuals[i] for i in indices[:population_size]]    
    return individuals
  
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
                population.append(Chromosome_complexes(graph, solution))
                solution = []
        population.append(Chromosome_complexes(graph, solution))
        return population
    except Exception as e:
        print(f"Erreur lors du chargement des complexes: {e}")
        return []

def save_complexes(complexes: List[Set[str]], filename: str) -> pd.DataFrame:
    data = []
    for i, complex in enumerate(complexes, 1):
        proteins = sorted(complex)
        data.append({'complex_id': i, 'proteins': ' '.join(proteins)})
    
    df = pd.DataFrame(data)
    df.to_csv(filename, sep='\t', index=False)
    return df

if __name__ == "__main__":
    # Chargement des données
    ppi_path = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\weighted_networks\weighted_STRING_levure.txt"
    complexes_path = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\Data\clean data\complexes\STRING_levure.txt"
    
    ppi_network = load_ppi(ppi_path)
    print(f"Réseau PPI chargé avec {len(ppi_network)} protéines")
    
    reference_complexes = load_complexes(complexes_path)
    print(f"{len(reference_complexes)} complexes de référence chargés")
    
    population = []
    complex_path = r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\results\initialization\complexes\detected_complexes_STRING_levure.txt"
    population1 = load_initialization_complexes(complex_path, ppi_network)
    population2 = random_initialization(reference_complexes, ppi_network, population_size=26, max_complexes=len(reference_complexes)*6, max_size=40)
    population.extend(population1)
    population.extend(population2)
    population = Population(population)
    print(f"la population initial est charge avec {len(population.individuals)} individus")
    
    
    # Configuration de l'algorithme
    ga = AlgorithmGA(
        reference_complexes,
        ppi_network,
        population=population,
        crossover=ComplexCrossover(ppi_network),
        mutator=TabuMutator(
            graph=ppi_network,
            tabu_tenure=20,
            tabu_iterations=100,
            rate=1,
            threshold_ratio = 400
        ),
        selection=population.tournament_selection,
        max_generations=1
    )
    
    # Exécution
    best_solution = ga.run()
    
    # Résultats
    save_complexes(best_solution.genes, r"C:\Users\PC\Documents\M2 HPC\PFE\PFE_CODE\results\mainMethod\OptClust-H\complexes\STRING_levure.txt")
    
    print("\nRésultats finaux:")
    print(f"Fitness: {best_solution._fitness}")
    print(f"Nombre de complexes: {len(best_solution.genes)}")
    print(f"metrics {compute_metrics(best_solution.genes, reference_complexes)}")
    