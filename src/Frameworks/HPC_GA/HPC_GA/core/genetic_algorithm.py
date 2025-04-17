from typing import List, Callable
from ..common.population import Population
from .operators import Crossover, Mutator
from ..common.chromosome import Chromosome

class GeneticAlgorithm:
    def __init__(
        self,
        population: Population,
        crossover: Crossover,
        mutator: Mutator,
        selection: Callable = Population.tournament_selection,
        max_generations: int = 100
    ):
        self.population = population
        self.crossover = crossover
        self.mutator = mutator
        self.selection = selection
        self.max_generations = max_generations
        self.history = {"best": [], "avg": []}

    def run(self) -> 'Chromosome':
        for gen in range(self.max_generations):
            self.population = self._evolve()
            self._update_history()
        return self.population.best()

    def _evolve(self) -> Population:
        offspring = [
            self.mutator(self.crossover(
                self.selection(self.population),
                self.selection(self.population)
            )) for _ in range(self.population.size)
        ]
        return Population(offspring)

    def _update_history(self):
        self.history["best"].append(self.population.best().fitness)
        self.history["avg"].append(self.population.average_fitness())