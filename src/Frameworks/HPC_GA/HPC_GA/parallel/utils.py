import numpy as np
from typing import List
from ..common.population import Population

def split_population(pop: 'Population', n_islands: int) -> List['Population']:
    indices = np.array_split(np.arange(len(pop.individuals)), n_islands)
    return [
        Population([pop.individuals[i] for i in idx])
        for idx in indices
    ]