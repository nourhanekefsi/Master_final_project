"""
fastGA - Framework for Genetic Algorithms

Expose les classes principales :
- GeneticAlgorithm : Algorithme séquentiel de base
- ParallelGeneticAlgorithm : Version parallèle
- Chromosome : Classe de base pour les solutions
"""

from .core.genetic_algorithm import GeneticAlgorithm
from .parallel.parallel_ga import ParallelGeneticAlgorithm
from .common.chromosome import Chromosome

__version__ = "1.0.0"
__all__ = [
    'GeneticAlgorithm',
    'ParallelGeneticAlgorithm',
    'Chromosome',
    'Population',
    'Crossover',
    'Mutator'
]