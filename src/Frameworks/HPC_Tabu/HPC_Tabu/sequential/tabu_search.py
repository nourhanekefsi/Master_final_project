from typing import Callable, List, Dict, Optional
from ..common.solution import Solution
from ..common.neighborhood import NeighborhoodGenerator
import numpy as np


class TabuSearch:
    def __init__(
        self,
        initial_solution: Solution,
        neighborhood_generator: NeighborhoodGenerator,
        tabu_tenure: int = 10,
        aspiration_criteria: Optional[List[Callable[[Solution], bool]]] = None,
        update_history: Optional[Callable[[Solution, Dict], Dict]] = None,
        intensification: bool = False,
        diversification: bool = False,
        max_iterations: int = 100,
        
    ):
        self.current_solution = initial_solution
        self.best_solution = initial_solution.copy()
        self.neighborhood = neighborhood_generator
        self.tabu_list = []
        self.tabu_tenure = tabu_tenure
        self.aspiration_criteria = aspiration_criteria or []
        self.intensification = intensification
        self.diversification = diversification
        self.max_iterations = max_iterations
        self.iterations = 0
        self.best_iteration = 0
        self._frequency: Dict[int, int] = {}  # Suivi des solutions visitées
        self._update_history = update_history
        self.history = {}
        
    def run(self) -> Solution:
        while not self._should_stop():
            neighbors = self._generate_neighbors()
            best_candidate = self._select_best_candidate(neighbors)
            
            if best_candidate:
                self._update_solution(best_candidate)
                self._update_tabu_list(best_candidate)
                self._update_frequency(best_candidate)
                if self._update_history is not None:
                    self.history = self._update_history(best_candidate, self.history)

            self.iterations += 1
        return self.best_solution

    def _generate_neighbors(self) -> List[Solution]:
        """Applique l'intensification si activée."""
        neighbors = self.neighborhood.generate(self.current_solution)
        if self.intensification:
            return [sol for sol in neighbors if self._frequency.get(hash(sol), 0) < 3]
        return neighbors

    def _should_stop(self) -> bool:
        """Gère les critères d'arrêt avec diversification."""
        if self.iterations >= self.max_iterations:
            return True
        if self.diversification and (self.iterations - self.best_iteration > 20):
            self._apply_diversification()
        return False

    def _apply_diversification(self):
        """Réinitialise avec une solution peu visitée."""
        if self._frequency:
            least_visited = min(self._frequency.items(), key=lambda x: x[1])[0]
            self.current_solution = least_visited
    
    
    def _update_tabu_list(self, best_candidate):
        """Met à jour la liste tabou avec le hash de la solution"""
        candidate_hash = hash(best_candidate)
        if len(self.tabu_list) >= self.tabu_tenure:
            self.tabu_list.pop(0)
        self.tabu_list.append(candidate_hash)
    
    def _select_best_candidate(self, neighbors):
        if not neighbors:
            return None
            
        # Évalue tous les voisins
        evaluated_neighbors = [(n, n.evaluate()) for n in neighbors]
        
        # Trie du meilleur (score le plus bas) au pire
        sorted_neighbors = sorted(evaluated_neighbors, key=lambda x: x[1])
        
        # Cherche le premier candidat non tabou ou qui satisfait les critères d'aspiration
        for candidate, score in sorted_neighbors:
            candidate_hash = hash(candidate)
            is_tabu = candidate_hash in self.tabu_list
            is_aspired = any(crit(candidate) for crit in self.aspiration_criteria)
            
            if not is_tabu or is_aspired:
                return candidate
                
        # Si tous sont tabous et aucun ne satisfait les critères, retourne le meilleur quand même
        return sorted_neighbors[0][0]
            
    
    def _update_solution(self, best_candidate):
        self.best_solution = best_candidate
    
    def _update_frequency(self, best_candidate):
        """Met à jour la fréquence de la solution visitée."""
        self._frequency[hash(best_candidate)] = self._frequency.get(hash(best_candidate), 0) + 1
        