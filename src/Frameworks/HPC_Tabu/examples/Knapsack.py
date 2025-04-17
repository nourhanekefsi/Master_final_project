from typing import List
from fastTabou.common import Solution, NeighborhoodGenerator, CompositeNeighborhood
from fastTabou.sequential import TabouSearch
from fastTabou.sequential.utils import frequency_aspiration

# 0. Données du problème (à adapter)
items = [(60, 10), (100, 20), (120, 30)]  # Format: (valeur, poids)
max_weight = 50  # Poids maximum du sac

# 1. Définir votre Solution pour le problème du sac à dos
class KnapsackSolution(Solution):
    def __init__(self, representation, items, max_weight):
        """
        :param representation: Liste de booléens (True = item sélectionné)
        :param items: Liste de tuples (valeur, poids)
        :param max_weight: Poids maximum autorisé
        """
        super().__init__(representation)
        self.items = items
        self.max_weight = max_weight

    def _evaluate(self) -> float:
        total_value = 0
        total_weight = 0
        for (value, weight), selected in zip(self.items, self.representation):
            if selected:
                total_value += value
                total_weight += weight
        # Pénalité si le poids dépasse la limite
        return total_value if total_weight <= self.max_weight else -1

    def copy(self):
        return KnapsackSolution(
            self.representation.copy(),
            self.items,
            self.max_weight
        )

# 2. Définir les mouvements de voisinage
class SwapNeighborhood(NeighborhoodGenerator):
    def generate(self, solution: KnapsackSolution) -> List[KnapsackSolution]:
        neighbors = []
        n = len(solution.representation)
        for i in range(n):
            for j in range(i + 1, n):
                new_repr = solution.representation.copy()
                new_repr[i], new_repr[j] = new_repr[j], new_repr[i]  # Échange deux items
                neighbors.append(KnapsackSolution(
                    new_repr, 
                    solution.items, 
                    solution.max_weight
                ))
        return neighbors

    @property
    def name(self) -> str:
        return "swap"

class InsertNeighborhood(NeighborhoodGenerator):
    def generate(self, solution: KnapsackSolution) -> List[KnapsackSolution]:
        neighbors = []
        n = len(solution.representation)
        for i in range(n):
            for j in range(n):
                if i != j:
                    new_repr = solution.representation.copy()
                    item = new_repr.pop(i)
                    new_repr.insert(j, item)
                    neighbors.append(KnapsackSolution(
                        new_repr,
                        solution.items,
                        solution.max_weight
                    ))
        return neighbors

    @property
    def name(self) -> str:
        return "insert"

# 3. Configuration et exécution
if __name__ == "__main__":
    # Solution initiale aléatoire
    import random
    initial_repr = [random.random() < 0.5 for _ in items]
    initial_solution = KnapsackSolution(initial_repr, items, max_weight)

    # Combinaison de mouvements
    neighborhood = CompositeNeighborhood([
        SwapNeighborhood(),
        InsertNeighborhood()
    ])

    # Création et exécution de la recherche tabou
    ts = TabouSearch(
        initial_solution=initial_solution,
        neighborhood_generator=neighborhood,
        tabu_tenure=5,
        intensification=True,
        diversification=True,
        aspiration_criteria=[frequency_aspiration],
        max_iterations=100
    )
    
    best_solution = ts.run()
    
    # Affichage des résultats
    print(f"Meilleure solution trouvée: {best_solution.representation}")
    print(f"Valeur totale: {best_solution.evaluate()}")
    print(f"Poids total: {sum(weight for (_, weight), selected in zip(items, best_solution.representation) if selected)}")