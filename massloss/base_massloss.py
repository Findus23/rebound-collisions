from abc import ABC, abstractmethod
from typing import Tuple


class Massloss(ABC):
    name: str

    @abstractmethod
    def estimate(self, alpha, velocity, projectile_mass, gamma) -> Tuple[float, float, float]:
        pass
