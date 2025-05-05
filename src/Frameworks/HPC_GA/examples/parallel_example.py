from HPC_GA.common.population import Population
from HPC_GA.common.chromosome import Chromosome
from HPC_GA.core.genetic_algorithm import GeneticAlgorithm
from HPC_GA.core.operators import Crossover, Mutator
import numpy as np

from HPC_GA.core.operators import GaussianMutator

from HPC_GA.core.operators import UniformCrossover
from HPC_GA.parallel.utils import split_population
class SphereChromosome(Chromosome):
    def evaluate(self):
        return -np.sum(self.genes**2)

# Création d'une population divisée en îlots
pop = Population(SphereChromosome(np.random.rand(3)))
island = pop

# Configuration
pga = GeneticAlgorithm(
    population=island,
    crossover=UniformCrossover(),
    mutator=GaussianMutator(rate=0.1),
    max_generations=50
)

best = pga.run()
print(f"Best solution: {best.genes} (Fitness: {best.fitness})")