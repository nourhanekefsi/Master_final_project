from abc import ABC, abstractmethod
import numpy as np
from ..common.chromosome import Chromosome

class Crossover(ABC):
    @abstractmethod
    def __call__(self, p1: 'Chromosome', p2: 'Chromosome') -> 'Chromosome':
        pass

class UniformCrossover(Crossover):
    def __call__(self, p1, p2):
        mask = np.random.randint(0, 2, size=len(p1.genes))
        new_genes = np.where(mask, p1.genes, p2.genes)
        return p1.__class__(new_genes)

class Mutator(ABC):
    def __init__(self, rate: float = 0.01):
        self.rate = rate

    @abstractmethod
    def __call__(self, chrom: 'Chromosome') -> 'Chromosome':
        pass

class GaussianMutator(Mutator):
    def __call__(self, chrom):
        noise = np.random.normal(0, 0.1, size=len(chrom.genes))
        mask = np.random.random(size=len(chrom.genes)) < self.rate
        new_genes = np.where(mask, chrom.genes + noise, chrom.genes)
        return chrom.__class__(new_genes)