import time
from math import radians

import numpy as np
from rebound import Simulation, Particle, Collision, Escape, reb_simulation_integrator_mercurius
from scipy.constants import astronomical_unit

from extradata import ExtraData, ParticleData
from merge import merge_particles
from radius_utils import radius
from utils import unique_hash, filename_from_argv

start = time.perf_counter()

sim = Simulation()

fn = filename_from_argv()

sim.units = ('yr', 'AU', 'kg')
sim.boundary = "open"
boxsize = 20
sim.configure_box(boxsize)
# sim.integrator = "mercurius"
# sim.dt = 1e-3
# sim.ri_whfast.safe_mode = 0
sim.collision = "direct"
# sim.ri_mercurius.hillfac = 3.
# sim.collision_resolve = "merge"

extradata = ExtraData()

i = 0
with open("initcon/conditions.input") as f:
    for line in f:
        if line.startswith("#") or line.startswith("ERROR") or line == "\n":
            continue
        columns = list(map(float, line.split()))
        hash = unique_hash()
        if len(columns) > 7:
            print(columns[7:])
            cmf, mmf, wmf = columns[7:]
            assert cmf + mmf + wmf == 1
            extradata.pdata[hash.value] = ParticleData(water_mass_fraction=wmf)
        else:
            wmf = 0
        if columns[1] == 0:  # that should not be needed, but nevertheless is
            part = Particle(m=columns[0], hash=hash)
        else:
            part = Particle(m=columns[0], a=columns[1], e=columns[2],
                            inc=radians(columns[3]), omega=radians(columns[4]),
                            Omega=radians(columns[5]), M=radians(columns[6]),
                            simulation=sim,
                            hash=hash,
                            r=radius(columns[0], wmf) / astronomical_unit
                            )
        sim.add(part)

        i += 1
        if i == 29:  # only use 29 objects to make results comparable with lie-simulation
            break

print(sim.N)
# from matplotlib import pyplot as plt
# OrbitPlot(sim,slices=1,color=True)
# plt.show()

max_n = sim.N
print("start")
tmax = 5e6
savesteps = 6000
times = np.linspace(0., tmax, savesteps)
sim.move_to_com()
# sim.exit_max_distance = 15

extradata.meta.tmax = tmax
extradata.meta.savesteps = savesteps
extradata.meta.max_n = max_n

try:
    fn.with_suffix(".bin").unlink()
except OSError:
    pass
# sim.automateSimulationArchive("water.bin", interval=1, deletefile=True)
for i, t in enumerate(times, start=1):
    try:
        sim.integrate(t)
        print(sim.dt)
        merc: reb_simulation_integrator_mercurius = sim.ri_mercurius
        print(merc._encounterN)
        print(merc.mode)
    except Collision:
        merge_particles(sim, extradata)
    # except Escape:
    #     print("something escaped")
    print(f"{i / savesteps * 100:.2f}% ({sim.N})")
    sim.simulationarchive_snapshot(str(fn.with_suffix(".bin")))
    extradata.meta.walltime = time.perf_counter() - start
    extradata.meta.cputime = time.process_time()
    extradata.meta.current_time = t
    extradata.meta.current_steps = i
    extradata.save(fn.with_suffix(".extra.json"))

# OrbitPlot(sim,slices=1,color=True)
# plt.show()
