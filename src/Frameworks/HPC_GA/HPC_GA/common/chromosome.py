from abc import ABC, abstractmethod
import numpy as np

class Chromosome(ABC):
    def __init__(self, genes):
        self.genes = genes
        self._fitness = None

    @abstractmethod
    def evaluate(self) -> float:
        pass

    @property
    def fitness(self) -> float:
        self._fitness = self.evaluate()
        return self._fitness

    def __lt__(self, other: 'Chromosome') -> bool:
        return self.fitness < other.fitness