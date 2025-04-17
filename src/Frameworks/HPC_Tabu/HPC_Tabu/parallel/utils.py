def default_stopping_condition(ts):
    """Condition d'arrêt par défaut : 50 itérations pour le parallèle"""
    return ts.iterations >= 50
from typing import Callable, List
from ..common.solution import Solution

# --- Aspiration Criteria ---
def basic_aspiration(candidate: Solution, best: Solution) -> bool:
    return candidate.evaluate() > best.evaluate()

def diversification_aspiration(candidate: Solution, best: Solution, freq_map: dict) -> bool:
    """Favorise les solutions rarement visitées."""
    return freq_map.get(hash(candidate), 0) < 3

# --- Intensification Helpers ---
def frequency_based_intensification(neighbors: List[Solution], freq_map: dict, threshold: int):
    return [sol for sol in neighbors if freq_map.get(hash(sol), 0) < threshold]

# --- Diversification Helpers ---
def restart_based_diversification(current_solution: Solution, freq_map: dict):
    """Réinitialise la recherche depuis une solution peu explorée."""
    least_visited = min(freq_map.items(), key=lambda x: x[1])[0]
    return least_visited