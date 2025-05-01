import numpy as np
from collections import defaultdict
from typing import Dict, List, Set, Tuple
import time
import math

class EWCA:
    def __init__(self, interaction_file: str, structural_similarity_threshold: float = 0.4):
        
        self.interaction_file = interaction_file
        self.structural_similarity_threshold = structural_similarity_threshold
        
        # Data structures
        self.relations: Dict[int, List[int]] = defaultdict(list)
        self.protein_to_id: Dict[str, int] = {}
        self.id_to_protein: Dict[int, str] = {}
        self.N = 0  # Number of proteins
        
        # Weights
        self.W1: np.ndarray = None
        self.Ws: np.ndarray = None
        
    def load_interactions(self) -> None:
        
        total_interactions = 0
        
        with open(self.interaction_file, 'r') as f:
            for line in f:
                if line.strip() == "":
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                    
                protein1, protein2 = parts[0], parts[1]
                
                # Assign IDs to proteins if not already assigned
                if protein1 not in self.protein_to_id:
                    self.protein_to_id[protein1] = self.N
                    self.id_to_protein[self.N] = protein1
                    self.N += 1
                    
                if protein2 not in self.protein_to_id:
                    self.protein_to_id[protein2] = self.N
                    self.id_to_protein[self.N] = protein2
                    self.N += 1
                    
                id1, id2 = self.protein_to_id[protein1], self.protein_to_id[protein2]
                
                # Add interaction (undirected)
                if id1 != id2 and id2 not in self.relations[id1]:
                    total_interactions += 1
                    self.relations[id1].append(id2)
                    self.relations[id2].append(id1)
        
        print(f"Total number of proteins: {self.N}")
        print(f"Total number of interactions: {total_interactions}")
        
        # Initialize weight matrices
        self.W1 = np.zeros((self.N, self.N))
        self.Ws = np.zeros((self.N, self.N))
    
    def calculate_jaccard_distance(self) -> None:
        
        for id1 in self.relations:
            neighbors1 = set(self.relations[id1])
            for id2 in self.relations[id1]:
                if id2 > id1:  # Process each pair only once
                    neighbors2 = set(self.relations[id2])
                    
                    if len(neighbors1) > 1 or len(neighbors2) > 1:
                        intersection = neighbors1 & neighbors2
                        union = neighbors1 | neighbors2
                        
                        if len(intersection) == 0:
                            weight = 0.0
                        else:
                            weight = len(intersection) / len(union)
                            
                        self.W1[id1, id2] = weight
                        self.W1[id2, id1] = weight
                    else:
                        self.W1[id1, id2] = 0.0
                        self.W1[id2, id1] = 0.0
    
    def calculate_ecv2_weights(self) -> None:
        
        for id1 in self.relations:
            neighbors1 = set(self.relations[id1])
            for id2 in self.relations[id1]:
                if id1 < id2:
                    neighbors2 = set(self.relations[id2])
                    common_neighbors = neighbors1 & neighbors2
                    
                    if len(common_neighbors) > 0:
                        sum_common_weight = sum(
                            self.W1[id1, iw] * self.W1[id2, iw] 
                            for iw in common_neighbors
                        )
                        h_second = sum_common_weight
                        h_second1 = (h_second + self.W1[id1, id2]) / (len(common_neighbors) + 1)
                        
                        self.Ws[id1, id2] = h_second1
                        self.Ws[id2, id1] = h_second1
        
        # Filter interactions based on weights
        self.filter_interactions()
    
    def filter_interactions(self) -> None:
        
        # Create a new relations dictionary with filtered interactions
        new_relations = defaultdict(list)
        
        for id1 in self.relations:
            for id2 in self.relations[id1]:
                if self.Ws[id1, id2] > 0.0:
                    new_relations[id1].append(id2)
        
        self.relations = {
            k: v for k, v in new_relations.items() 
            if len(v) > 1
        }
    
    def calculate_structural_similarity(self, id1: int, id2: int) -> float:
       
        neighbors1 = set(self.relations.get(id1, []))
        neighbors2 = set(self.relations.get(id2, []))
        
        # Include the proteins themselves in their neighborhoods
        neighbors1.add(id1)
        neighbors2.add(id2)
        
        common_neighbors = neighbors1 & neighbors2
        denominator = math.sqrt(len(neighbors1) * len(neighbors2))
        
        return len(common_neighbors) / denominator if denominator > 0 else 0.0
    
    def detect_core_complexes(self) -> Dict[int, List[int]]:
        
        core_complexes = {}
        
        for id1 in self.relations:
            complex_members = {id1}
            
            for id2 in self.relations[id1]:
                if id1 < id2:
                    similarity = self.calculate_structural_similarity(id1, id2)
                    if similarity > self.structural_similarity_threshold:
                        complex_members.add(id2)
            
            if len(complex_members) >= 2:
                core_complexes[id1] = list(complex_members)
        
        return core_complexes
    
    def find_attachments(self, core_complexes: Dict[int, List[int]]) -> Dict[int, List[int]]:
    
        complexes = {}
        
        for core_id, core_members in core_complexes.items():
            # Find all neighboring proteins not in the core
            attachment_candidates = set()
            for protein in core_members:
                attachment_candidates.update(
                    p for p in self.relations.get(protein, []) 
                    if p not in core_members
                )
            
            if not attachment_candidates:
                complexes[core_id] = core_members
                continue
            
            # Calculate average edge weight within the core
            in_edges, _ = self.clustering_coefficient(core_members)
            avg_core_edges = 2 * in_edges / len(core_members)
            
            # Classify attachment candidates
            local_attachments = set()
            overlapping_attachments = set()
            
            for candidate in attachment_candidates:
                in_weight = 0.0
                out_weight = 0.0
                count = 0
                
                for neighbor in self.relations.get(candidate, []):
                    if neighbor in core_members:
                        in_weight += self.Ws[candidate, neighbor]
                        count += 1
                    else:
                        out_weight += self.Ws[candidate, neighbor]
                
                if count >= 2:
                    if in_weight <= out_weight:
                        if in_weight >= 0.5 * avg_core_edges:
                            overlapping_attachments.add(candidate)
                    else:
                        local_attachments.add(candidate)
            
            # Combine all attachments
            all_attachments = local_attachments | overlapping_attachments
            complexes[core_id] = core_members + list(all_attachments)
        
        return complexes
    
    def clustering_coefficient(self, proteins: List[int]) -> Tuple[float, int]:
        """
        Calculate clustering coefficient for a set of proteins.
        Returns (sum_of_in_edges, count_of_edges)
        """
        protein_set = set(proteins)
        sum_edges = 0.0
        count = 0
        
        for i, id1 in enumerate(proteins):
            for id2 in proteins[i+1:]:
                if id2 in self.relations.get(id1, []):
                    sum_edges += self.Ws[id1, id2]
                    count += 1
        
        return sum_edges, count
    
    def overlap_score(self, set1: Set[int], set2: Set[int]) -> float:
        
        intersection = set1 & set2
        union = set1 | set2
        return len(intersection) / len(union) if union else 0.0
    
    def filter_redundant_complexes(self, complexes: Dict[int, List[int]], 
                                 overlap_threshold: float = 1.0) -> Dict[int, List[int]]:
    
        visited = set()
        result = {}
        complex_items = list(complexes.items())
        
        for i, (id1, members1) in enumerate(complex_items):
            if id1 in visited:
                continue
                
            current_members = set(members1)
            
            for j, (id2, members2) in enumerate(complex_items[i+1:], i+1):
                if id2 in visited:
                    continue
                    
                score = self.overlap_score(current_members, set(members2))
                if score >= overlap_threshold:
                    visited.add(id2)
                    current_members.update(members2)
            
            result[id1] = list(current_members)
        
        return result
    
    def save_results(self, complexes: Dict[int, List[int]], output_file: str) -> None:
        
        with open(output_file, 'w') as f:
            for complex_id, members in complexes.items():
                if len(members) >= 3:  # Only save complexes with at least 3 proteins
                    proteins = [self.id_to_protein[pid] for pid in members]
                    line = f"{complex_id}\t{' '.join(proteins)}\n"
                    f.write(line)
    
    def run(self, output_file: str) -> None:
        
        start_time = time.time()
        
        print("Loading interactions...")
        self.load_interactions()
        
        print("Calculating Jaccard distances...")
        self.calculate_jaccard_distance()
        
        print("Calculating ECV2 weights...")
        self.calculate_ecv2_weights()
        
        print("Detecting core complexes...")
        core_complexes = self.detect_core_complexes()
        
        print("Finding attachment proteins...")
        complexes = self.find_attachments(core_complexes)
        
        print("Filtering redundant complexes...")
        filtered_complexes = self.filter_redundant_complexes(complexes)
        
        print("Saving results...")
        self.save_results(filtered_complexes, output_file)
        
        elapsed = time.time() - start_time
        print(f"EWCA completed in {elapsed:.2f} seconds")
        print(f"Results saved to {output_file}")


# Example usage
if __name__ == "__main__":
    # Initialize with your interaction file and desired parameters
    ewca = EWCA(
        interaction_file="/Users/ryham/Documents/mémoireM2/Master_final_project/Data/clean data/weighted_networks/weighted_STRING_humain.txt",
        structural_similarity_threshold=0.4
    )
    
    # Run the algorithm
    ewca.run("/Users/ryham/Documents/mémoireM2/Master_final_project/Data/results/tmp/STRING_humain_complexes_ewca2.txt")