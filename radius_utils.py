from math import pi

ice_density = 0.917 / 1000 * 100 ** 3  # kg/m^3
basalt_density = 2.7 / 1000 * 100 ** 3  # kg/m^3


def core_radius(total_mass, water_fraction, density):
    core_mass = total_mass * (1 - water_fraction)
    return (core_mass / density * 3 / 4 / pi) ** (1 / 3)


def total_radius(total_mass, water_fraction, density, inner_radius):
    mantle_mass = total_mass * water_fraction
    return (mantle_mass / density * 3 / 4 / pi + inner_radius ** 3) ** (1 / 3)


def radius(mass: float, water_fraction: float) -> float:
    """
    :return: radius in Meter
    """
    return total_radius(mass, water_fraction, ice_density, core_radius(mass, water_fraction, basalt_density))
