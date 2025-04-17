# fastGA
Framework for Genetic Algorithms with parallel island model.

## Installation
```bash
pip install fastGA  # Version standard
pip install "fastGA[parallel]"  # Avec parallélisme
```

## Quickstart
```python
from fastGA import GeneticAlgorithm, Chromosome

class MySolution(Chromosome):
    def evaluate(self):
        return -sum(x**2 for x in self.genes)

ga = GeneticAlgorithm(...)
best = ga.run()
```