from ctypes import c_uint32
from random import randint
from typing import Dict

from numpy import linalg
from rebound import Simulation, Orbit, OrbitPlot, Particle
from scipy.constants import pi, gravitational_constant

from extradata import ExtraData
from utils import solar_mass


def random_hash() -> c_uint32:
    """
    returns a (hopefully) unique 32 bit integer to be used as a particle hash

    when a collision occurs and ruins the output data, please complain to the universe
    """
    return c_uint32(randint(0, 2 ** 32 - 1))


def unique_hash(ed: ExtraData) -> c_uint32:
    ed.meta.hash_counter += 1
    return c_uint32(ed.meta.hash_counter)


def innermost_period(sim: Simulation) -> float:
    """
    returns the orbital period in years of the innerpost object
    for comparison with symplectic time step
    """
    minp = None
    mina = float("inf")
    for p in sim.particles[1:]:
        if abs(p.a) < mina:
            mina = abs(p.a)
            minp = p
    orb: Orbit = minp.orbit
    return abs(orb.P)


def third_kepler_law(orbital_period: float):
    return (
                   gravitational_constant * solar_mass / (4 * pi ** 2)
                   * orbital_period ** 2
           ) ** (1 / 3)


def total_momentum(sim: Simulation) -> float:
    total = 0
    for p in sim.particles:
        total += linalg.norm(p.vxyz) * p.m
    return total


def total_mass(sim: Simulation) -> float:
    total = 0
    for p in sim.particles:
        total += p.m
    return total


def show_orbits(sim: Simulation):
    from matplotlib import pyplot as plt
    OrbitPlot(sim, slices=1, color=True)
    plt.show()

def reorder_particles(sim: Simulation, ed: ExtraData) -> None:
    """
    probably not needed anymore
    """
    exit()
    particles_by_hash: Dict[int, Particle] = {}
    hashes = []
    suns = []
    gas_giants = []
    embryos = []
    planetesimals = []
    original_N = sim.N
    p: Particle
    for p in sim.particles:
        hash_value = p.hash.value
        # save a copy of the particles
        particles_by_hash[hash_value] = p.copy()
        type = ed.pd(p).type
        if type == "sun":
            suns.append(hash_value)
            # keep sun in the simulation to avoid empty particles list
            continue
        elif type == "gas giant":
            gas_giants.append(hash_value)
        elif type == "embryo":
            embryos.append(hash_value)
        elif type == "planetesimal":
            planetesimals.append(hash_value)
        else:
            raise ValueError(f"unknown type: {type}")
        hashes.append(p.hash)

    for hash in hashes:
        sim.remove(hash=hash)
    ordered_particles = gas_giants + embryos + planetesimals
    assert len(suns) + len(ordered_particles) == original_N
    particle: Particle
    for h in ordered_particles:
        particle = particles_by_hash[h]
        sim.add(particle)
    # print(list(sim.particles))
    # exit()
    sim.N_active = len(suns) + len(ordered_particles) - len(planetesimals)
    # mark objects > N_active as testparticle_type = 1
    # this means they become semi-active bodies
    # TODO: double-check meaning
    sim.testparticle_type = 1
    assert sim.N == original_N
