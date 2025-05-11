from typing import Callable, List, Dict, Optional, Tuple
from ..common.solution import Solution
from ..common.neighborhood import NeighborhoodGenerator
import numpy as np
import random
from collections import defaultdict


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
        diversification_frequency: int = 20,
        intensification_threshold: int = 3,
        patience: int = 15
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
        self._frequency: Dict[int, int] = defaultdict(int)  # Suivi des solutions visitées
        self._update_history = update_history
        self.history = {}
        self.diversification_frequency = diversification_frequency
        self.intensification_threshold = intensification_threshold
        self.patience = patience
        self.no_improvement_count = 0
        
    def run(self) -> Solution:
        while not self._should_stop():
            neighbors = self._generate_neighbors()
            best_candidate = self._select_best_candidate(neighbors)
            
            if best_candidate:
                self._update_current_solution(best_candidate)
                self._update_best_solution(best_candidate)
                self._update_tabu_list(best_candidate)
                self._update_frequency(best_candidate)
                if self._update_history is not None:
                    self.history = self._update_history(best_candidate, self.history)

            self.iterations += 1
        return self.best_solution

    def _generate_neighbors(self) -> List[Solution]:
        """Génère les voisins avec stratégie d'intensification si activée."""
        neighbors = self.neighborhood.generate(self.current_solution)
        
        if self.intensification:
            # Intensification: favorise les solutions peu explorées
            neighbors = [
                sol for sol in neighbors 
                if self._frequency[hash(sol)] < self.intensification_threshold
            ]
            
        if not neighbors:  # Si intensification a filtré tous les voisins
            neighbors = self.neighborhood.generate(self.current_solution)
            
        return neighbors

    def _should_stop(self) -> bool:
        """Détermine si la recherche doit s'arrêter."""
        # Critère d'arrêt principal
        if self.iterations >= self.max_iterations:
            return True
            
        # Critère de stagnation
        if self.iterations - self.best_iteration > self.patience:
            return True
            
        # Diversification périodique
        if (self.diversification and 
            self.iterations > 0 and 
            self.iterations % self.diversification_frequency == 0 and
            len(self._frequency) > 0):
            self._apply_diversification()
            
        return False

    def _apply_diversification(self):
        """Applique une diversification en choisissant une solution peu visitée."""
        if not self._frequency:
            return
            
        # Sélectionne les solutions les moins fréquentées
        min_freq = min(self._frequency.values())
        candidates = [h for h, f in self._frequency.items() if f == min_freq]
        
        # Choisit aléatoirement parmi les moins visitées
        chosen_hash = random.choice(candidates)
        
        # Trouve la solution correspondante (nécessite un mécanisme de stockage)
        # Note: Cette partie nécessiterait un cache des solutions visitées
        # Pour l'instant, nous réinitialisons simplement avec la meilleure solution
        self.current_solution = self.best_solution.copy()
        self.no_improvement_count = 0
    
    def _update_tabu_list(self, best_candidate: Solution):
        """Met à jour la liste tabou avec une solution ou un attribut."""
        candidate_hash = self._get_move_hash(best_candidate)
        
        if len(self.tabu_list) >= self.tabu_tenure:
            self.tabu_list.pop(0)
            
        self.tabu_list.append(candidate_hash)
    
    def _get_move_hash(self, candidate: Solution) -> int:
        """Génère un hash pour le mouvement entre current et candidate."""
        # Pour les complexes protéiques, on peut hasher la différence entre les solutions
        diff_added = frozenset(candidate.representation - self.current_solution.representation)
        diff_removed = frozenset(self.current_solution.representation - candidate.representation)
        return hash((diff_added, diff_removed))
    
    def _select_best_candidate(self, neighbors: List[Solution]) -> Optional[Solution]:
        """Sélectionne le meilleur candidat selon les critères tabou."""
        if not neighbors:
            return None
            
        # Évalue tous les voisins
        evaluated_neighbors = [(n, n.evaluate()) for n in neighbors]
        
        # Trie du meilleur (score le plus élevé) au pire
        sorted_neighbors = sorted(evaluated_neighbors, key=lambda x: x[1], reverse=True)
        
        # Cherche le premier candidat non tabou ou qui satisfait les critères d'aspiration
        for candidate, score in sorted_neighbors:
            move_hash = self._get_move_hash(candidate)
            is_tabu = move_hash in self.tabu_list
            is_aspired = any(crit(candidate) for crit in self.aspiration_criteria)
            
            # Critère d'aspiration: permet de surpasser tabou si solution est meilleure que la meilleure globale
            is_best_aspired = score > self.best_solution.evaluate()
            
            if (not is_tabu) or is_aspired or is_best_aspired:
                if score > self.current_solution.evaluate():
                    return candidate
                    
        # Si tous sont tabous et aucun ne satisfait les critères, retourne None
        return None
            
    def _update_frequency(self, best_candidate: Solution):
        """Met à jour la fréquence des solutions visitées."""
        self._frequency[hash(best_candidate)] += 1
        
    def _update_best_solution(self, candidate: Solution):
        """Met à jour la meilleure solution globale."""
        if candidate.evaluate() > self.best_solution.evaluate():
            self.best_solution = candidate.copy()
            self.best_iteration = self.iterations
            self.no_improvement_count = 0
        else:
            self.no_improvement_count += 1
            
    def _update_current_solution(self, candidate: Solution):
        """Met à jour la solution courante."""
        self.current_solution = candidate.copy()