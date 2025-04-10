from multiprocessing import Pool
from typing import List
from ..sequential.tabou_search import TabouSearch
from ..common import Solution, NeighborhoodGenerator

class ParallelTabouSearch:
    def __init__(
        self,
        initial_solutions: List[Solution],
        neighborhood_generator: NeighborhoodGenerator,
        tabu_tenure: int = 10,
        n_processes: int = 4,
        **kwargs  # Accepte intensification/diversification
    ):
        self.initial_solutions = initial_solutions
        self.neighborhood = neighborhood_generator
        self.tabu_tenure = tabu_tenure
        self.n_processes = n_processes
        self.kwargs = kwargs  # Passage des paramètres avancés

    def _worker(self, initial_solution: Solution) -> Solution:
        ts = TabouSearch(
            initial_solution=initial_solution,
            neighborhood_generator=self.neighborhood,
            tabu_tenure=self.tabu_tenure,
            **self.kwargs
        )
        return ts.run()

    def run(self) -> Solution:
        with Pool(self.n_processes) as pool:
            results = pool.map(self._worker, self.initial_solutions)
        return max(results, key=lambda x: x.evaluate())