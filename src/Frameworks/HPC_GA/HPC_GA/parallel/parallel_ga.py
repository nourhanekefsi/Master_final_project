from ..common.chromosome import Chromosome
from multiprocessing import Pool
from typing import List
from ..core.genetic_algorithm import GeneticAlgorithm
from ..common.population import Population

class ParallelGeneticAlgorithm:
    def __init__(
        self,
        populations: List[Population],
        crossover,
        mutator,
        migration_interval: int = 5,
        n_processes: int = 4,
        **kwargs
    ):
        self.populations = populations
        self.migration_interval = migration_interval
        self.n_processes = n_processes
        self.ga_kwargs = {
            "crossover": crossover,
            "mutator": mutator,
            **kwargs
        }

    def _run_island(self, pop: Population) -> Population:
        ga = GeneticAlgorithm(population=pop, **self.ga_kwargs)
        ga.run()
        return ga.population

    def run(self) -> 'Chromosome':
        with Pool(self.n_processes) as pool:
            for _ in range(self.ga_kwargs.get("max_generations", 100) // self.migration_interval):
                self.populations = pool.map(self._run_island, self.populations)
                self._migrate()
        return max([pop.best() for pop in self.populations])

    def _migrate(self):
        migrants = [pop.best() for pop in self.populations]
        for i, pop in enumerate(self.populations):
            pop.individuals[-1] = migrants[(i + 1) % len(self.populations)]