import time
from math import radians

import numpy as np
from numpy import linalg
from rebound import Simulation, Particle, reb_simulation_integrator_mercurius, NoParticles, Collision, Orbit, Escape
from scipy.constants import astronomical_unit

from extradata import ExtraData, ParticleData
from merge import merge_particles, handle_escape
from radius_utils import radius
from utils import unique_hash, filename_from_argv


def innermost_period(sim: Simulation) -> float:
    minp = None
    mina = float("inf")

    for p in sim.particles[1:]:
        if p.a < mina:
            mina = p.a
            minp = p
    orb: Orbit = minp.orbit
    return orb.P


start = time.perf_counter()

sim = Simulation()

fn = filename_from_argv()

sim.units = ('yr', 'AU', 'kg')
# sim.boundary = "open"
# boxsize = 100
# sim.configure_box(boxsize)
sim.exit_max_distance = 30
sim.integrator = "mercurius"
# sim.collision = 'line'
# sim.ri_whfast.coordinates = "democraticheliocentric"
sim.dt = 1e-2
sim.ri_ias15.min_dt = 0.0001 / 365
sim.ri_whfast.safe_mode = 0
sim.collision = "direct"
sim.ri_mercurius.hillfac = 3.
sim.ri_whfast.corrector = 11
# sim.collision_resolve = "merge"

extradata = ExtraData()

i = 0
with open("initcon/conditions_many.input") as f:
    for line in f:
        if line.startswith("#") or line.startswith("ERROR") or line == "\n":
            continue
        columns = list(map(float, line.split()))
        hash = unique_hash()
        if len(columns) > 7:
            # print(columns[7:])
            cmf, mmf, wmf = columns[7:]
            total_fractions = cmf + mmf + wmf
            if total_fractions != 1:
                diff = 1 - total_fractions
                print(f"fractions don't add up by {diff}")
                print("adding rest to cmf")
                cmf += diff
            assert cmf + mmf + wmf - 1 <= 1e-10
            extradata.pdata[hash.value] = ParticleData(water_mass_fraction=wmf)
        else:
            wmf = 0
        if columns[1] == 0:  # that should not be needed, but nevertheless is
            part = Particle(m=columns[0], hash=hash)
        else:
            part = Particle(
                m=columns[0], a=columns[1], e=columns[2],
                inc=radians(columns[3]), omega=columns[4],
                Omega=columns[5], M=columns[6],
                simulation=sim,
                hash=hash,
                r=radius(columns[0], wmf) / astronomical_unit
            )
        sim.add(part)

        i += 1
        # if i == 29:  # only use 29 objects to make results comparable with lie-simulation
        #     break
sim.move_to_com()

print(sim.N)
assert sim.dt < innermost_period(sim) / 20
# exit()
# from matplotlib import pyplot as plt
# OrbitPlot(sim,slices=1,color=True)
# plt.show()

max_n = sim.N
print("start")
tmax = 300
savesteps = 300
times = np.linspace(0., tmax, savesteps)
extradata.energy.set_initial_energy(sim.calculate_energy())

extradata.meta.tmax = tmax
extradata.meta.savesteps = savesteps
extradata.meta.max_n = max_n
abort = False
try:
    fn.with_suffix(".bin").unlink()
except OSError:
    pass
# sim.automateSimulationArchive("water.bin", interval=1, deletefile=True)
for i, t in enumerate(times, start=1):
    try:
        sim.integrate(t, exact_finish_time=0)
        print("dt", sim.dt)
        merc: reb_simulation_integrator_mercurius = sim.ri_mercurius
        print("t", t)
        print(merc.mode)
    except Collision:
        merge_particles(sim, extradata)
    except Escape:
        handle_escape(sim, extradata)
    except NoParticles:
        print("No Particles left")
        abort = True
    # except Escape:
    #     print("something escaped")
    print(f"{i / savesteps * 100:.2f}% ({sim.N})")
    sim.simulationarchive_snapshot(str(fn.with_suffix(".bin")))
    extradata.meta.walltime = time.perf_counter() - start
    extradata.meta.cputime = time.process_time()
    extradata.meta.current_time = t
    extradata.meta.current_steps = i
    extradata.energy.add_energy_value(sim.calculate_energy())
    total = 0
    p: Particle
    for p in sim.particles:
        total += linalg.norm(p.vxyz) * p.m
    print("total", total)
    extradata.save(fn.with_suffix(".extra.json"))
    # assert sim.dt < innermost_period(sim) / 20
    print("fraction", innermost_period(sim) / 20)

    if abort:
        exit(1)

# OrbitPlot(sim,slices=1,color=True)
# plt.show()
