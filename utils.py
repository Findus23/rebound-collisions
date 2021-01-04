from ctypes import c_uint32
from pathlib import Path
from random import randint
from sys import argv
from typing import Tuple, Dict

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy import linalg
from rebound import Simulation, Orbit, Particle, OrbitPlot

from extradata import ExtraData


def unique_hash() -> c_uint32:
    """
    returns a (hopefully) unique 32 bit integer to be used as a particle hash

    when a collision occurs and ruins the output data, please complain to the universe
    """
    return c_uint32(randint(0, 2 ** 32 - 1))


def innermost_period(sim: Simulation) -> float:
    """
    returns the orbital period in years of the innerpost object
    for comparison with symplectic time step
    """
    minp = None
    mina = float("inf")

    for p in sim.particles[1:]:
        if p.a < mina:
            mina = p.a
            minp = p
    orb: Orbit = minp.orbit
    return orb.P


def total_impulse(sim: Simulation) -> float:
    total = 0
    p: Particle
    for p in sim.particles:
        total += linalg.norm(p.vxyz) * p.m
    return total


def show_orbits(sim: Simulation):
    from matplotlib import pyplot as plt
    OrbitPlot(sim, slices=1, color=True)
    plt.show()


def clamp(n: float, smallest: float, largest: float) -> float:
    assert smallest < largest
    return max(smallest, min(n, largest))


def filename_from_argv(argument: str = None) -> Path:
    if len(argv) < 2:
        raise ValueError("specify filename")
    if argument:
        fn = argument
    else:
        fn = argv[1]
    fn = fn.replace(".bin", "").replace(".meta.json", "")
    if fn.endswith("."):
        fn = fn[:-1]
    return Path(fn.replace(".bin", "").replace(".meta.json", ""))


def create_figure() -> Tuple[Figure, Axes]:
    """
    helper function for matplotlib OOP interface with proper typing
    """
    fig: Figure = plt.figure()
    ax: Axes = fig.gca()
    return fig, ax


def reorder_particles(sim: Simulation, ed: ExtraData) -> None:
    particles_by_hash: Dict[str, Particle] = {}
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
