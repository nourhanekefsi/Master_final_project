from HPC_GA import Chromosome, Population, ParallelGeneticAlgorithm
from HPC_GA.core.operators import UniformCrossover, GaussianMutator
from HPC_GA.utils import split_population  # Import split_population if it exists in fastGA
import numpy as np

class SphereChromosome(Chromosome):
    def evaluate(self):
        return -np.sum(self.genes**2)

# Création d'une population divisée en îlots
pop = Population([SphereChromosome(np.random.rand(3)) for _ in range(100)])
islands = split_population(pop, n_islands=4)

# Configuration
pga = ParallelGeneticAlgorithm(
    populations=islands,
    crossover=UniformCrossover(),
    mutator=GaussianMutator(rate=0.1),
    migration_interval=10,
    n_processes=4,
    max_generations=50
)

best = pga.run()
print(f"Best solution: {best.genes} (Fitness: {best.fitness})")