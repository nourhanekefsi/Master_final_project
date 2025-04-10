"""
fastTabou - Framework for Tabu Search algorithms with sequential and parallel implementations
"""

from .sequential import TabouSearch
from .parallel import ParallelTabouSearch

__version__ = "0.1.0"
__all__ = ['TabouSearch', 'ParallelTabouSearch']