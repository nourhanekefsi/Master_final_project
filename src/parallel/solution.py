from typing import Callable, Dict, List, Optional, Set
from HPC_Tabu.common.solution import Solution as baseSolution
from metrics import FF
from HPC_Tabu.common.solution import Solution

class Solution(baseSolution):
    """Solution pour la recherche tabou avec gestion standardisée des protéines"""
    def __init__(self, representation: Set[str], graph: Dict[str, Dict[str, float]]):
        # Standardisation des IDs de protéines
        representation = set(sorted(set(representation)))
        super().__init__(representation)
        self.graph = graph
        self._value = FF(representation, graph)
        
    def _evaluate(self) -> float:
        return FF(self.representation, self.graph)
    
    def copy(self) -> 'Solution':
        return Solution(sorted(set([p for p in self.representation])), self.graph)
    
    
    def __hash__(self):
        # Crée un hash stable basé sur la représentation
        return hash(frozenset(self.representation))
    
    def __eq__(self, other):
        if not isinstance(other, Solution):
            return False
        return sorted(self.representation) == sorted(other.representation)