from typing import List, Optional, Set, Dict
import numpy as np
from HPC_GA.common.chromosome import Chromosome as BaseChromosome
from metrics import FF, FS_fitness

class Chromosome(BaseChromosome):
    """Chromosome représentant un ensemble de complexes protéiques avec standardisation des IDs"""
    def __init__(self, graph: Dict[str, Dict[str, float]], genes: Optional[List[Set[str]]] = None):
        
        super().__init__(genes)
        self.graph = graph
        self._fits = np.array([FF(gene, self.graph) for gene in self.genes])
        self._fitness = np.sum(self._fits)
        
    def evaluate(self) -> float:
        if self._fitness is None:
            self._fitness = np.sum(self._fits)
        return self._fitness