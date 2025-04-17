from typing import List
import numpy as np
from ..common.chromosome import Chromosome

class Population:
    def __init__(self, individuals: List['Chromosome']):
        self.individuals = individuals

    def best(self) -> 'Chromosome':
        return max(self.individuals)

    def average_fitness(self) -> float:
        return np.mean([ind.fitness for ind in self.individuals])

    def tournament_selection(self, k: int = 3) -> 'Chromosome':
        candidates = np.random.choice(self.individuals, size=k)
        return max(candidates)