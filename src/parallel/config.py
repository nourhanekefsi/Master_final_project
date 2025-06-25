from dataclasses import dataclass

@dataclass
class GAConfig:
    max_generations: int = 100
    population_size: int = 50
    crossover_rate: float = 0.9
    mutation_rate: float = 0.1
    tournament_size: int = 3

@dataclass
class ParallelConfig:
    model_type: str = 'island'  # 'island', 'master_slave', or 'hybrid'
    n_islands: int = 4
    migration_interval: int = 5
    migration_size: int = 2
    metric_workers: int = None