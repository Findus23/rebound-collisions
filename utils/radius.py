from math import pi


class PlanetaryRadius:
    # all densities are till_rho_0
    ice_density = 917  # kg/m^3
    basalt_density = 2700  # kg/m^3
    iron_density = 7800  # kg/m^3

    def __init__(self, total_mass: float, water_fraction: float, core_fraction: float):
        self.total_mass = total_mass
        self.water_fraction = water_fraction
        self.core_fraction = core_fraction

    @property
    def mantle_fraction(self) -> float:
        return 1 - self.core_fraction - self.water_fraction

    @property
    def core_mass(self) -> float:
        return self.total_mass * self.core_fraction

    @property
    def mantle_mass(self) -> float:
        return self.total_mass * self.mantle_fraction

    @property
    def water_mass(self) -> float:
        return self.total_mass * self.water_fraction

    @property
    def core_radius(self) -> float:
        return (self.core_mass / self.iron_density * 3 / 4 / pi) ** (1 / 3)

    @property
    def mantle_radius(self) -> float:
        return (self.mantle_mass / self.basalt_density * 3 / 4 / pi + self.core_radius ** 3) ** (1 / 3)

    @property
    def total_radius(self) -> float:
        return (self.water_mass / self.ice_density * 3 / 4 / pi + self.mantle_radius ** 3) ** (1 / 3)


if __name__ == '__main__':
    testradius = PlanetaryRadius(1e21, 0.1, 0)
    assert testradius.mantle_radius == 430127.0069140495
    assert testradius.total_radius == 472683.51839024643
    testradius = PlanetaryRadius(1e21, 0.05, 0)
    assert testradius.total_radius == 459494.5247690266
