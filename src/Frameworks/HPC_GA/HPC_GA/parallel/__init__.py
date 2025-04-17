"""
Implémentations parallèles :
- ParallelGeneticAlgorithm : Modèle d'îlots parallèles
- Utilitaires pour MPI/multiprocessing
"""

from .parallel_ga import ParallelGeneticAlgorithm
from .utils import split_population

__all__ = [
    'ParallelGeneticAlgorithm',
    'split_population'
]