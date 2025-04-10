from typing import Callable
from ..common.solution import Solution

# --- Critères d'aspiration ---
def default_aspiration(candidate: Solution, best: Solution) -> bool:
    return candidate.evaluate() > best.evaluate()

def frequency_aspiration(candidate: Solution, best: Solution, freq: dict) -> bool:
    """Accepte si la solution est rarement visitée."""
    return freq.get(hash(candidate), 0) < 2

# --- Conditions d'arrêt ---
def diversification_stop(algorithm, patience: int = 20) -> bool:
    return algorithm.iterations - algorithm.best_iteration > patience