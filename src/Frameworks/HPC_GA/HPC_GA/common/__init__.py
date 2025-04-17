"""
Définitions communes :
- Chromosome : Représentation des solutions
- Population : Gestion des ensembles de solutions
"""

from .chromosome import Chromosome
from .population import Population

__all__ = [
    'Chromosome',
    'Population'
]