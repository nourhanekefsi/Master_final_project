import sys
from math import sqrt
import numpy as np
from collections import defaultdict
import copy
from numpy.matlib import zeros

class MPCC:
    def __init__(self):
        self.label_id = {}
        self.id_label = {}
        self.relations = defaultdict(list)
        self.weights = None
        
    def load_interactions(self, ppi_file):
        """Load PPI data with pre-calculated weights (handles header)"""
        index = 0
        protein_list = []
        relations = defaultdict(list)
        weights = {}
        has_header = True  # Modifier à False si le fichier n'a pas d'en-tête
        
        with open(ppi_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                
                # Skip header line if exists
                if has_header and index == 0:
                    has_header = False  # On a traité l'en-tête
                    continue  # Passe à la ligne suivante
                
                if len(parts) >= 3:  # Format attendu : protein1 protein2 weight
                    p1, p2, weight = parts[0], parts[1], float(parts[2])
                    
                    # Gestion des identifiants
                    if p1 not in self.label_id:
                        self.label_id[p1] = index
                        self.id_label[index] = p1
                        protein_list.append(p1)
                        index += 1
                    if p2 not in self.label_id:
                        self.label_id[p2] = index
                        self.id_label[index] = p2
                        protein_list.append(p2)
                        index += 1
                    
                    # Ajout des interactions
                    id1, id2 = self.label_id[p1], self.label_id[p2]
                    if id2 not in relations[id1] and id1 != id2:
                        relations[id1].append(id2)
                        relations[id2].append(id1)
                        weights[(id1, id2)] = weight
                        weights[(id2, id1)] = weight
        
        self.Protein_num = index
        self.relations = relations
        self.weights = weights
        return self.label_id, self.id_label, self.Protein_num, self.relations, protein_list
    
    def remove_false_positives(self):
        """Remove interactions with zero weights"""
        for id in self.relations:
            neighbors = self.relations[id]
            new_neighbors = [it for it in neighbors if self.weights.get((id, it), 0) > 0]
            self.relations[id] = new_neighbors
        return self.relations
    
    def calculate_topology_scores(self, network, N):
        """Calculate topology similarity scores"""
        Topology_weight = zeros((N, N))
        for id in network:
            neighbors = network[id]
            for it in neighbors:
                neighbors1 = network[it]
                score = self.TO_ij_score(neighbors, neighbors1)
                Topology_weight[id, it] = score
                Topology_weight[it, id] = score
        return Topology_weight
    
    def TO_ij_score(self, list1, list2):
        """Calculate topological overlap score"""
        common = set(list1) & set(list2)
        top_sum = len(common)
        top_down = min(len(list1), len(list2))
        
        if top_down > 0:
            score = (top_sum/sqrt(len(list1)*len(list2)) + 
                     top_sum/top_down + 
                     2*top_sum/(len(list1)+len(list2))) / 3
        else:
            score = 0.0
        return score
    
    def detect_seeds(self, id_list, network, weight_network):
        """Detect seed proteins for complex formation"""
        seeds = defaultdict(list)
        Graph_avgdensity = 0.0
        num = 0
        
        # Calculate average graph density
        for id in id_list:
            neighbors = network.get(id, [])
            if id not in neighbors:
                neighbors.append(id)
            if len(neighbors) >= 2:
                density, _ = self.graph_density(neighbors, network, weight_network)
                Graph_avgdensity += density
                num += 1
                
        if num > 0:
            Graph_avgdensity /= num
        
        # Identify seed proteins
        visit_node = []
        for id1 in id_list:
            if id1 not in visit_node:
                neighbors1 = network.get(id1, [])
                if id1 not in neighbors1:
                    neighbors1.append(id1)
                
                if len(neighbors1) >= 2:
                    # Find core proteins
                    protein_cores = [id1]
                    density1, weight_sum1 = self.graph_density(neighbors1, network, weight_network)
                    avg_deg = 2 * weight_sum1 / len(neighbors1)
                    
                    for ih in neighbors1:
                        neighbor_ih = network.get(ih, [])
                        sum_w = sum(weight_network[ih,it] for it in neighbor_ih 
                                   if it in neighbors1)
                        if sum_w >= avg_deg:
                            protein_cores.append(ih)
                    
                    # Find essential proteins connected to core
                    left_proteins = list(set(neighbors1) - set(protein_cores))
                    essential_cores = []
                    for ib in left_proteins:
                        neighbors_ib = network.get(ib, [])
                        common = set(neighbors_ib) & set(protein_cores)
                        if len(common) >= 2 and len(common) > 0.5*len(protein_cores):
                            essential_cores.append(ib)
                    
                    protein_cores = list(set(essential_cores) | set(protein_cores))
                    density2, _ = self.graph_density(protein_cores, network, weight_network)
                    
                    if density2 > Graph_avgdensity:
                        seeds[id1] = protein_cores
                        visit_node.extend(protein_cores)
        
        return seeds
    
    def graph_density(self, graph, network, weight):
        """Calculate graph density"""
        weight_sum = 0.0
        if len(graph) >= 2:
            for id in graph:
                neighbors = network.get(id, [])
                for it in neighbors:
                    if id < it and it in graph:
                        weight_sum += weight[id, it]
            density = 2 * weight_sum / (len(graph) * (len(graph)-1))
        else:
            density = 0.0
        return density, weight_sum
    
    def identify_complexes(self, seeds, network, weight_network):
        """Identify protein complexes from seeds"""
        complexes = defaultdict(list)
        k = 1
        
        for seed_id in seeds:
            initial_graph = seeds[seed_id]
            if len(initial_graph) >= 2:
                # Grow the complex iteratively
                protein_complex, _ = self.grow_complex(initial_graph, network, weight_network, 1)
                protein_complex = list(set(protein_complex))
                
                if len(protein_complex) >= 3:
                    complexes[k] = protein_complex
                    k += 1
                    
        return complexes
    
    def grow_complex(self, initial_graph, relations, weight, iterations):
        """Grow complex iteratively by adding/removing nodes"""
        subgraph_new = self.add_nodes(initial_graph, relations, weight)
        subgraph_new = self.remove_nodes(subgraph_new, relations, weight)
        
        if set(subgraph_new) != set(initial_graph) and iterations < 10:
            return self.grow_complex(subgraph_new, relations, weight, iterations+1)
        else:
            return subgraph_new, iterations
    
    def add_nodes(self, graph, relations, weight):
        """Add nodes to the complex that improve its score"""
        add_nodes = []
        for id in graph:
            neighbors = relations.get(id, [])
            for it in neighbors:
                if it not in graph and it not in add_nodes:
                    add_nodes.append(it)
        
        while add_nodes:
            start_score = self.graph_entropy(graph, relations, weight)
            max_node = self.find_max_node(add_nodes, graph, relations, weight)
            candidate = graph + [max_node]
            new_score = self.graph_entropy(candidate, relations, weight)
            
            if new_score > start_score:
                graph.append(max_node)
                add_nodes.remove(max_node)
            else:
                break
                
        return graph
    
    def remove_nodes(self, graph, relations, weight):
        """Remove nodes from complex that improve its score"""
        if len(graph) <= 2:
            return graph
            
        remove_list = [id for id in graph 
                      if not set(relations.get(id, [])).issubset(graph)]
        
        while remove_list:
            start_score = self.graph_entropy(graph, relations, weight)
            min_node = self.find_min_node(remove_list, graph, relations, weight)
            candidate = [x for x in graph if x != min_node]
            new_score = self.graph_entropy(candidate, relations, weight)
            
            if new_score > start_score:
                graph.remove(min_node)
                remove_list.remove(min_node)
                if len(graph) <= 2:
                    break
            else:
                break
                
        return graph
    
    def find_max_node(self, nodes, graph, relations, weight):
        """Find node with maximum connection weight to current graph"""
        max_node = nodes[0]
        max_weight = sum(weight[max_node, it] for it in relations.get(max_node, []) 
                        if it in graph)
        
        for id in nodes[1:]:
            sum_w = sum(weight[id, it] for it in relations.get(id, []) 
                    if it in graph)
            if sum_w > max_weight:
                max_node, max_weight = id, sum_w
                
        return max_node
    
    def find_min_node(self, nodes, graph, relations, weight):
        """Find node with minimum connection weight to current graph"""
        min_node = nodes[0]
        min_weight = sum(weight[min_node, it] for it in relations.get(min_node, []) 
                    if it in graph)
        
        for id in nodes[1:]:
            sum_w = sum(weight[id, it] for it in relations.get(id, []) 
                    if it in graph)
            if sum_w < min_weight:
                min_node, min_weight = id, sum_w
                
        return min_node
    
    def graph_entropy(self, graph, relations, weight):
        """Calculate graph clustering score"""
        if len(graph) < 2:
            return 0.0
            
        weight_in = 0.0
        weight_out = 0.0
        count_out = 0
        
        for id in graph:
            neighbors = relations.get(id, [])
            inner_sum = sum(weight[id, it] for it in neighbors if it in graph and it > id)
            outer_sum = sum(weight[id, it] for it in neighbors if it not in graph)
            
            weight_in += inner_sum
            if outer_sum > 0:
                weight_out += outer_sum
                count_out += 1
                
        weight_in = 2 * weight_in / len(graph)
        weight_out = weight_out / count_out if count_out > 0 else 0.0
        
        if weight_in + weight_out > 0:
            return weight_in / (weight_in + weight_out)
        return 0.0
    
    def read_essential_proteins(self, essential_file):
        """Read list of essential proteins"""
        essential_proteins = []
        with open(essential_file, "r") as f:
            for line in f:
                protein = line.strip().split()[0]
                if protein in self.label_id:
                    essential_proteins.append(self.label_id[protein])
        return list(set(essential_proteins))
    
    def overlap_score(self, clusterA, clusterB):
        """Calculate overlap score between two clusters"""
        intersect = len(set(clusterA) & set(clusterB))
        return (intersect * intersect) / (len(clusterA) * len(clusterB))
    
    def filter_redundant(self, complexes, threshold=0.8):
        """Filter redundant complexes based on overlap"""
        to_remove = set()
        complex_ids = list(complexes.keys())
        
        for i in range(len(complex_ids)):
            id1 = complex_ids[i]
            if id1 in to_remove:
                continue
                
            for j in range(i+1, len(complex_ids)):
                id2 = complex_ids[j]
                if id2 in to_remove:
                    continue
                    
                score = self.overlap_score(complexes[id1], complexes[id2])
                if score >= threshold:
                    to_remove.add(id2)
        
        return {k:v for k,v in complexes.items() if k not in to_remove}
    
    def sort_complexes(self, complexes, scores=None):
        """Sort complexes by size or score with robust key handling"""
        if scores:
            # Ensure all score keys exist in complexes
            valid_keys = [k for k in scores.keys() if k in complexes]
            sorted_ids = sorted(((k, scores[k]) for k in valid_keys), 
                              key=lambda x: x[1], reverse=True)
            return {i+1: complexes[k] for i, (k,v) in enumerate(sorted_ids)}
        else:
            sorted_ids = sorted(complexes.items(), 
                              key=lambda x: len(x[1]), 
                              reverse=True)
            return {i+1: v for i, (k,v) in enumerate(sorted_ids)}
    
    
    def run_static(self, ppi_file, output_file="resultats_complexes.txt"):
        """Main workflow with improved error handling"""
        try:
            self.load_interactions(ppi_file)
            self.remove_false_positives()
            print(f"Total proteins: {len(self.id_label)}")
            print(f"Total interactions: {sum(len(v) for v in self.relations.values())//2}")
            
            N = len(self.id_label)
            topo_weights = self.calculate_topology_scores(self.relations, N)
            combined_weights = zeros((N, N))
            
            # Create weight matrix
            weight_matrix = zeros((N, N))
            for (i,j), w in self.weights.items():
                weight_matrix[i,j] = w
            
            # Combine weights
            for i in range(N):
                for j in range(N):
                    if i < j and weight_matrix[i,j] > 0 and topo_weights[i,j] > 0:
                        combined_weights[i,j] = 2 * weight_matrix[i,j] * topo_weights[i,j] / (weight_matrix[i,j] + topo_weights[i,j])
                        combined_weights[j,i] = combined_weights[i,j]
            
            # Detect complexes
            seeds = self.detect_seeds(list(self.relations.keys()), self.relations, combined_weights)
            complexes = self.identify_complexes(seeds, self.relations, combined_weights)
            
            # Calculate scores with key consistency check
            complex_scores = {}
            final_complexes = defaultdict(list)
            count = 1
            
            for cid in list(complexes.keys()):  # Explicit conversion to list
                if len(complexes[cid]) >= 3:
                    final_complexes[count] = complexes[cid]
                    complex_scores[count] = self.graph_entropy(complexes[cid], self.relations, combined_weights)
                    count += 1
            print(f"Complexes keys: {set(final_complexes.keys())}")
            print(f"Scores keys: {set(complex_scores.keys())}")
            missing_keys = set(complex_scores.keys()) - set(final_complexes.keys())
            if missing_keys:
                print(f"WARNING: {len(missing_keys)} score keys missing from complexes")
            # Filter and sort
            filtered = self.filter_redundant(final_complexes)
            sorted_complexes = self.sort_complexes(filtered, complex_scores)
            
            # Save results
            with open(output_file, "w") as f:
                for cid in sorted_complexes:
                    proteins = [self.id_label.get(pid, "UNKNOWN") for pid in sorted_complexes[cid]]  # Safer protein name lookup
                    line = f"{cid} {' '.join(proteins)}\n"
                    f.write(line)
            
            return sorted_complexes
            
        except Exception as e:
            print(f"Error processing network: {str(e)}")
            raise

# Example usage with error handling
if __name__ == "__main__":
    try:
        mpcc = MPCC()
        complexes = mpcc.run_static(
            ppi_file="/Users/ryham/Documents/mémoireM2/Master_final_project/Data/clean data/weighted_networks/weighted_BIOGRID_humain.txt",
            output_file="/Users/ryham/Documents/mémoireM2/Master_final_project/Data/results/tmp/BIOGRID_humain_complexes_mpcc.txt"
        )
        print(f"Successfully identified {len(complexes)} protein complexes")
    except Exception as e:
        print(f"Failed to process: {str(e)}")