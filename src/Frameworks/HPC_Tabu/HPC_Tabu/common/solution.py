from abc import ABC, abstractmethod
from typing import Any, Dict

class Solution(ABC):
    def __init__(self, representation: Any):
        self.representation = representation
        self._value = None
        self.metadata: Dict[str, Any] = {}  # Pour intensification/diversification

    @abstractmethod
    def _evaluate(self) -> float:
        pass

    def evaluate(self) -> float:
        if self._value is None:
            self._value = self._evaluate()
        return self._value

    @abstractmethod
    def copy(self) -> 'Solution':
        pass

    def __eq__(self, other) -> bool:
        return self.representation == other.representation

    def __hash__(self) -> int:
        return hash(str(self.representation))