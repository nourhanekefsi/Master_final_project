"""
fastTabou - Framework for Tabu Search algorithms with sequential and parallel implementations
"""

from .sequential import TabuSearch
from .parallel import ParallelTabuSearch

__version__ = "0.1.0"
__all__ = ['TabuSearch', 'ParallelTabuSearch']