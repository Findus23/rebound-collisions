from typing import Tuple

from numpy.random import default_rng

from . import Massloss


class LeiZhouMassloss(Massloss):
    """
    implement water loss estimation from
    "On the formation of terrestrial planets between two massive planets: The case of 55 Cancri"
    Lei Zhou, Rudolf Dvorak, Li-Yong Zhou
    https://arxiv.org/abs/2105.10105
    """
    name = "leizhou"

    def __init__(self):
        self.rng = default_rng()
        self.delta_w_low = 0.01
        self.delta_m_low = 0.01
        self.delta_w_up = 0.08
        self.delta_m_up = 0.10

    def estimate(self, alpha, velocity, projectile_mass, gamma) -> Tuple[float, float, float]:
        rand_water = self.rng.random()
        rand_mass = self.rng.random()
        water_loss = self.delta_w_low + rand_water * (self.delta_w_up - self.delta_w_low)
        mantle_loss = self.delta_m_low + rand_mass * (self.delta_m_up - self.delta_m_low)
        # the paper only considers objects with a core and water shell
        # so we assume the same mass retention for mantle and core
        core_loss = mantle_loss
        return 1 - water_loss, 1 - mantle_loss, 1 - core_loss
