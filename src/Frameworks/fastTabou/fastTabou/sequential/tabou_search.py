from typing import Callable, List, Dict, Optional
from ..common.solution import Solution
from ..common.neighborhood import NeighborhoodGenerator

class TabouSearch:
    def __init__(
        self,
        initial_solution: Solution,
        neighborhood_generator: NeighborhoodGenerator,
        tabu_tenure: int = 10,
        aspiration_criteria: Optional[List[Callable]] = None,
        intensification: bool = False,
        diversification: bool = False,
        max_iterations: int = 100
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

    def run(self) -> Solution:
        while not self._should_stop():
            neighbors = self._generate_neighbors()
            best_candidate = self._select_best_candidate(neighbors)
            
            if best_candidate:
                self._update_solution(best_candidate)
                self._update_tabu_list(best_candidate)
                self._update_frequency(best_candidate)

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