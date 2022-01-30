from dataclasses import dataclass

from rebound import Particle
from scipy.constants import kilo, astronomical_unit

habitable_zone_inner = 0.75
habitable_zone_outer = 1.5


def is_potentially_habitable(planet: Particle):
    return habitable_zone_inner <= planet.a <= habitable_zone_outer


@dataclass
class ReferencePlanet:
    a: float  # km
    e: float
    mass: float

    @property
    def a_au(self):
        return self.a * kilo / astronomical_unit


# from NASA HORIZONS
# https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND=%27Mercury+Barycenter%27&OBJ_DATA=%27YES%27&MAKE_EPHEM=%27YES%27&EPHEM_TYPE=ELEMENTS&CENTER=%400&START_TIME=1961-08-06&STOP_TIME=1961-08-07&STEP_SIZE=%272+days%27
inner_solar_system_data = {
    "mercury": ReferencePlanet(mass=3.302e23, a=6.132591901350513E+07, e=2.330561345288357E-01),
    "venus": ReferencePlanet(mass=48.685e23, a=1.084671627283828E+08, e=3.297503114845414E-03),
    "earth": ReferencePlanet(mass=5.97219e24, a=1.473393857505501E+08, e=2.365126062626348E-02),
    "mars": ReferencePlanet(mass=6.4171e23, a=2.277188818140753E+08, e=9.685171612265037E-02),
}
