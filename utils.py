import os
import socket
import subprocess
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
from scipy.constants import pi, gravitational_constant
from setproctitle import setproctitle

from extradata import ExtraData

solar_radius = 6.957e8
solar_mass = 1.9885e+30


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
    """
    probably not needed anymore
    """
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


def git_hash() -> str:
    output = subprocess.run(["git", "rev-parse", "--verify", "HEAD"], capture_output=True)
    return output.stdout.decode()


def process_friendlyness(fn: Path) -> None:
    if socket.gethostname() == "standpc":
        # only handle other computers specially
        return
    setproctitle(f"[rebound-watersim] [{fn.stem}] read /home/winklerl23/sim-info.txt for more information")
    os.nice(5)
