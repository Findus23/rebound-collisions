from rebound import Particle

habitable_zone_inner=0.75
habitable_zone_outer=1.5


def is_potentially_habitable(planet: Particle):
    return habitable_zone_inner <= planet.a <= habitable_zone_outer
