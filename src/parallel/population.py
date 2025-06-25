from typing import List
import random
import numpy as np
from chromosome import Chromosome
from HPC_GA.common.population import Population as BasePopulation

class MyPopulation(BasePopulation):
    def __init__(self, individuals: List[Chromosome]):
        super().__init__(individuals)
        
    def tournament_selection(self, k: int = 2) -> Chromosome:
        candidates = np.random.choice(self.individuals, size=k)
        return max(candidates, key=lambda ind: ind.evaluate())
    
    def integrate_migrants(self, migrants: List[Chromosome]):
        self.individuals = sorted(self.individuals, key=lambda x: x.fitness)[len(migrants):] + migrants
    
    def get_best_migrants(self, n: int) -> List[Chromosome]:
        return sorted(self.individuals, key=lambda x: x.fitness, reverse=True)[:n]
    
    def _my_update_population(self, worst_parent : Chromosome, child:Chromosome)->List['Chromosome']:
        if worst_parent.evaluate() < child.evaluate():
            self.individuals.remove(worst_parent)
            self.individuals.append(child)
            
            