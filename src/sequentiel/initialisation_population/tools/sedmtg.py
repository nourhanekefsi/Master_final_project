import random
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import math
from dataclasses import dataclass
from tqdm import tqdm  # Pour la barre de progression

@dataclass
class ProteinNetwork:
    """Represents a protein-protein interaction network"""
    relations: Dict[str, List[str]]  # Protein ID to list of neighbors
    weights: Dict[Tuple[str, str], float]  # Edge weights
    
    @classmethod
    def from_weighted_network(cls, network_file: str, has_header: bool = True):
        """Create network from weighted interaction file"""
        relations = defaultdict(list)
        weights = {}
        
        # Ajout de la barre de progression pour le chargement
        print("Loading network data...")
        with open(network_file) as f:
            # Skip header if present
            if has_header:
                next(f)
            
            # Compter le nombre de lignes pour la barre de progression
            lines = list(f)
            for line in tqdm(lines, desc="Processing interactions"):
                if line.strip():
                    parts = line.strip().split()
                    # Ensure we have exactly 3 columns
                    if len(parts) != 3:
                        continue
                    
                    try:
                        protein1, protein2, weight = parts[0], parts[1], float(parts[2])
                    except (ValueError, IndexError):
                        continue
                    
                    # Add to relations
                    relations[protein1].append(protein2)
                    relations[protein2].append(protein1)
                    
                    # Add weights (undirected)
                    weights[(protein1, protein2)] = weight
                    weights[(protein2, protein1)] = weight
        
        return cls(relations=dict(relations), weights=weights)


class SEDMTG:
    """Implementation of the SEDMTG algorithm for protein complex detection"""
    
    def __init__(self, network: ProteinNetwork, iterations: int = 5):
        self.network = network
        self.iterations = iterations
    
    def find_seeds(self) -> Dict[str, float]:
        """Step 1: Generate seed queue based on node scores"""
        node_scores = {}
        proteins = list(self.network.relations.keys())
        
        # Ajout de la barre de progression
        for protein in tqdm(proteins, desc="Finding seeds"):
            neighbors = self.network.relations.get(protein, [])
            if len(neighbors) >= 2:
                degree_weight = sum(self.network.weights.get((protein, neighbor), 0) 
                                  for neighbor in neighbors)
                
                subgraph = neighbors + [protein]
                density, _ = self._calculate_density(subgraph)
                score = density * degree_weight
                node_scores[protein] = score
        
        # Shuffle seeds to introduce randomness
        seeds = list(node_scores.keys())
        random.shuffle(seeds)
        
        return {protein: node_scores[protein] for protein in seeds}
    
    def detect_complexes(self) -> Dict[int, List[str]]:
        """Main method to detect protein complexes"""
        all_complexes = defaultdict(list)
        
        # Barre de progression pour les itérations principales
        for _ in tqdm(range(self.iterations), desc="Detecting complexes", total=self.iterations):
            seeds = self.find_seeds()
            seed_proteins = list(seeds.keys())
            
            # Randomly sample seeds
            sampled_seeds = random.sample(seed_proteins, len(seed_proteins))
            seed_dict = {i: protein for i, protein in enumerate(sampled_seeds)}
            
            # Detect complexes from these seeds
            complexes, _, _, _ = self._detect_complexes_from_seeds(seed_dict)
            
            # Filter redundant complexes
            filtered_complexes = self._filter_redundant_complexes(complexes)
            
            # Store results
            for complex_id, proteins in filtered_complexes.items():
                all_complexes[len(all_complexes)] = proteins
        
        return dict(all_complexes)
    
    def save_complexes_to_file(self, complexes: Dict[int, List[str]], output_file: str):
        """Save detected complexes to a file with the required format"""
        with open(output_file, 'w') as f:
            # Tri des complexes par taille (décroissant) comme dans votre exemple
            sorted_complexes = sorted(complexes.items(), key=lambda x: len(x[1]), reverse=True)
            
            # Réindexation pour avoir des IDs séquentiels
            for new_id, (old_id, proteins) in enumerate(sorted_complexes, 1):
                # Format: ID suivi d'une tabulation puis les protéines séparées par des espaces
                f.write(f"{new_id}\t{' '.join(proteins)}\n")
    
    def _detect_complexes_from_seeds(self, seeds: Dict[int, str]) -> Tuple[Dict[int, List[str]], int, Dict[int, float], float]:
        """Step 2-3: Form initial clusters and extend/correct them"""
        complexes = defaultdict(list)
        complex_scores = {}
        total_iterations = 0
        visited = set()
        
        # Calculate average score threshold
        avg_score = self._calculate_average_seed_score(seeds)
        
        for seed_id, protein in seeds.items():
            if protein not in visited:
                neighbors = self.network.relations.get(protein, [])
                initial_cluster = neighbors + [protein]
                
                # Only proceed if cluster meets basic criteria
                if len(set(initial_cluster)) >= 3:
                    score, _, _ = self._composite_score(initial_cluster)
                    
                    if score >= avg_score:
                        # Form initial cluster
                        initial_graph = [protein]
                        
                        # Extend and refine cluster
                        final_cluster, iterations = self._refine_cluster(initial_graph)
                        total_iterations += iterations
                        
                        # Check final cluster quality
                        final_score, _, _ = self._composite_score(final_cluster)
                        if final_score > avg_score and len(final_cluster) >= 3:
                            complexes[len(complexes)] = sorted(list(set(final_cluster)))
                            complex_scores[len(complexes)] = final_score
                            visited.update(final_cluster)
        
        avg_iterations = total_iterations / len(complexes) if complexes else 0
        return complexes, len(complexes), complex_scores, avg_iterations
    
    def _refine_cluster(self, initial_cluster: List[str], max_iterations: int = 50) -> Tuple[List[str], int]:
        """Step 3: Extend and correct the cluster iteratively"""
        current_cluster = initial_cluster.copy()
        iterations = 0
        
        while iterations < max_iterations:
            iterations += 1
            old_cluster = current_cluster.copy()
            
            # Extension phase
            current_cluster = self._extend_cluster(current_cluster)
            
            # Correction phase
            current_cluster = self._correct_cluster(current_cluster)
            
            # Check for convergence
            if set(old_cluster) == set(current_cluster):
                break
        
        return list(set(current_cluster)), iterations
    
    def _extend_cluster(self, cluster: List[str]) -> List[str]:
        """Add relevant proteins to the cluster"""
        neighbors = self._get_cluster_neighbors(cluster)
        extended_cluster = cluster.copy()
        
        while neighbors:
            neighbors = list(set(neighbors))
            current_score, sum_weight, _ = self._composite_score(extended_cluster)
            
            # Find best node to add
            best_node, best_score = self._find_best_addition(extended_cluster, neighbors, sum_weight)
            
            # Calculate addition criteria
            common_neighbors = len(set(self.network.relations.get(best_node, [])) & set(extended_cluster))
            threshold = current_score * len(extended_cluster)
            
            # Check if addition improves the cluster
            if best_score > current_score and common_neighbors >= threshold:
                extended_cluster.append(best_node)
                neighbors.remove(best_node)
                neighbors.extend(self._get_new_neighbors(extended_cluster, best_node))
            else:
                break
        
        return extended_cluster
    
    def _correct_cluster(self, cluster: List[str]) -> List[str]:
        """Remove irrelevant proteins from the cluster"""
        if len(cluster) <= 2:
            return cluster.copy()
        
        corrected_cluster = cluster.copy()
        removable = self._get_removable_nodes(corrected_cluster)
        
        while removable:
            current_score, sum_weight, _ = self._composite_score(corrected_cluster)
            
            # Find worst node to remove
            worst_node, worst_score = self._find_worst_removal(corrected_cluster, removable, sum_weight)
            
            # Calculate removal criteria
            common_neighbors = len(set(self.network.relations.get(worst_node, [])) & set(corrected_cluster))
            threshold = current_score * len(corrected_cluster)
            
            # Check if removal improves the cluster
            if worst_score > current_score and common_neighbors <= threshold:
                corrected_cluster.remove(worst_node)
                removable.remove(worst_node)
            else:
                break
            
            if len(corrected_cluster) <= 2:
                break
        
        return corrected_cluster
    
    def _filter_redundant_complexes(self, complexes: Dict[int, List[str]], overlap_threshold: float = 0.8) -> Dict[int, List[str]]:
        """Step 4: Filter redundant protein complexes"""
        # Sort complexes by size (descending)
        sorted_complexes = sorted(complexes.items(), key=lambda x: len(x[1]), reverse=True)
        
        filtered = {}
        excluded = set()
        
        for i, (idx1, complex1) in enumerate(sorted_complexes):
            if idx1 not in excluded:
                filtered[idx1] = complex1
                
                # Compare with remaining complexes
                for j in range(i + 1, len(sorted_complexes)):
                    idx2, complex2 = sorted_complexes[j]
                    if idx2 not in excluded:
                        overlap = self._calculate_overlap(set(complex1), set(complex2))
                        if overlap >= overlap_threshold:
                            excluded.add(idx2)
        
        # Reindex filtered complexes
        return {i: proteins for i, (_, proteins) in enumerate(filtered.items())}
    
    # Helper methods for calculations
    def _calculate_density(self, subgraph: List[str]) -> Tuple[float, float]:
        """Calculate subgraph density and total weight"""
        if len(subgraph) <= 2:
            return 0.0, 0.0
        
        total_weight = 0.0
        nodes = list(set(subgraph))
        
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                protein1, protein2 = nodes[i], nodes[j]
                total_weight += self.network.weights.get((protein1, protein2), 0)
        
        density = 2 * total_weight / (len(nodes) * (len(nodes) - 1))
        return density, total_weight
    
    def _calculate_modularity(self, subgraph: List[str]) -> Tuple[float, float]:
        """Calculate subgraph modularity and external weight"""
        nodes = list(set(subgraph))
        internal_weight = 0.0
        external_weight = 0.0
        
        # Calculate internal weight
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                protein1, protein2 = nodes[i], nodes[j]
                internal_weight += self.network.weights.get((protein1, protein2), 0)
        
        # Calculate external weight
        for protein in nodes:
            for neighbor in self.network.relations.get(protein, []):
                if neighbor not in nodes:
                    external_weight += self.network.weights.get((protein, neighbor), 0)
        
        total_weight = internal_weight + external_weight
        modularity = internal_weight / total_weight if total_weight > 0 else 0.0
        return modularity, external_weight
    
    def _composite_score(self, subgraph: List[str]) -> Tuple[float, float, float]:
        """Combined score for subgraph evaluation"""
        density, sum_weight = self._calculate_density(subgraph)
        modularity, sum_out_weight = self._calculate_modularity(subgraph)
        
        if density == 0 or modularity == 0:
            score = (modularity + density) / 2
        else:
            score = (modularity + density + math.sqrt(density * modularity)) / 3
        
        return score, sum_weight, sum_out_weight
    
    def _calculate_average_seed_score(self, seeds: Dict[int, str]) -> float:
        """Calculate average score of seed proteins"""
        total_score = 0.0
        count = 0
        
        for seed_id, protein in seeds.items():
            neighbors = self.network.relations.get(protein, [])
            if len(neighbors) >= 2:
                score, _, _ = self._composite_score(neighbors + [protein])
                total_score += score
                count += 1
        
        return total_score / count if count > 0 else 0.0
    
    def _get_cluster_neighbors(self, cluster: List[str]) -> List[str]:
        """Get neighbors of a cluster not already in it"""
        neighbors = []
        for protein in cluster:
            for neighbor in self.network.relations.get(protein, []):
                if neighbor not in cluster:
                    neighbors.append(neighbor)
        return list(set(neighbors))
    
    def _get_removable_nodes(self, cluster: List[str]) -> List[str]:
        """Get nodes that are candidates for removal"""
        removable = []
        cluster_set = set(cluster)
        
        for protein in cluster:
            neighbors = set(self.network.relations.get(protein, []))
            if not neighbors.issubset(cluster_set):
                removable.append(protein)
        
        return list(set(removable))
    
    def _find_best_addition(self, cluster: List[str], candidates: List[str], current_weight: float) -> Tuple[str, float]:
        """Find the best node to add to the cluster"""
        best_node = None
        best_score = -1
        
        for node in candidates:
            # Calculate new weight if node were added
            new_weight = current_weight
            for protein in cluster:
                new_weight += self.network.weights.get((node, protein), 0)
            
            # Calculate score (using simplified formula)
            score = 2 * new_weight / (len(cluster) + 1)
            
            if score > best_score:
                best_score = score
                best_node = node
        
        return best_node, best_score
    
    def _find_worst_removal(self, cluster: List[str], candidates: List[str], current_weight: float) -> Tuple[str, float]:
        """Find the worst node to remove from the cluster"""
        worst_node = None
        worst_score = float('inf')
        
        for node in candidates:
            # Calculate new weight if node were removed
            new_weight = current_weight
            for protein in cluster:
                if protein != node:
                    new_weight -= self.network.weights.get((node, protein), 0)
            
            # Calculate score (using simplified formula)
            score = 2 * new_weight / (len(cluster))
            
            if score < worst_score:
                worst_score = score
                worst_node = node
        
        return worst_node, worst_score
    
    def _get_new_neighbors(self, cluster: List[str], new_node: str) -> List[str]:
        """Get new neighbors after adding a node to the cluster"""
        new_neighbors = []
        for neighbor in self.network.relations.get(new_node, []):
            if neighbor not in cluster:
                new_neighbors.append(neighbor)
        return new_neighbors
    
    def _calculate_overlap(self, set1: Set[str], set2: Set[str]) -> float:
        """Calculate overlap score between two protein sets"""
        intersection = len(set1 & set2)
        return (intersection * intersection) / (len(set1) * len(set2))


if __name__ == "__main__":
    # Load network data
    network_file = "/Users/ryham/Documents/mémoireM2/Master_final_project/Data/clean data/weighted_networks/tmp/GO_weighted_STRING_humain.txt"
    output_file = "/Users/ryham/Documents/mémoireM2/Master_final_project/Data/results/tmp/STRING_humain_complexes_sedmtg_v2.txt"
    
    # Create network (assuming the file has a header)
    print("Initializing network...")
    network = ProteinNetwork.from_weighted_network(network_file, has_header=True)
    
    # Run SEDMTG algorithm
    print("Running SEDMTG algorithm...")
    sedmtg = SEDMTG(network, iterations=10)
    protein_complexes = sedmtg.detect_complexes()
    
    # Save results to file
    print("Saving results...")
    sedmtg.save_complexes_to_file(protein_complexes, output_file)
    
    # Print summary
    print("\nProcessing complete!")
    print(f"Detected {len(protein_complexes)} protein complexes")
    print(f"Results saved to: {output_file}")