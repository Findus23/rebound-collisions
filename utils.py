from ctypes import c_uint32
from pathlib import Path
from random import randint
from sys import argv

from numpy import linalg
from rebound import Simulation, Orbit, Particle, OrbitPlot


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
