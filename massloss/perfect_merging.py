from typing import Tuple

from massloss import Massloss


class PerfectMerging(Massloss):
    name = "perfectmerging"

    def estimate(self, alpha, velocity, projectile_mass, gamma) -> Tuple[float, float, float]:
        return 1, 1, 1
