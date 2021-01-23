import re
import time
from math import radians
from pathlib import Path
from shutil import copy

from rebound import Simulation, Particle, NoParticles, Escape, \
    SimulationArchive
from rebound.simulation import POINTER_REB_SIM, reb_collision
from scipy.constants import astronomical_unit, mega

from extradata import ExtraData, ParticleData
from merge import merge_particles, handle_escape
from radius_utils import radius
from utils import unique_hash, filename_from_argv, innermost_period, total_impulse, process_friendlyness

MIN_TIMESTEP_PER_ORBIT = 20
PERFECT_MERGING = False
INITCON_FILE = Path("initcon/conditions_many.input")


def main(fn: Path):
    start = time.perf_counter()

    if not fn.with_suffix(".bin").exists():
        # set up a fresh simulation
        sim = Simulation()

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
        tmax = 200 * mega
        num_savesteps = 20000
        per_savestep = tmax / num_savesteps
        extradata = ExtraData()
        # times = np.linspace(0., tmax, savesteps)
        extradata.meta.tmax = tmax
        extradata.meta.per_savestep = per_savestep
        extradata.meta.num_savesteps = num_savesteps
        extradata.meta.perfect_merging = PERFECT_MERGING

        initcon = INITCON_FILE.read_text()
        num_embryos = int(re.search(r"Generated (\d+) minor bodies", initcon, re.MULTILINE).group(1))
        num_planetesimals = int(re.search(r"Generated (\d+) small bodies", initcon, re.MULTILINE).group(1))
        sim.N_active = num_embryos + 3
        i = 1
        for line in initcon.split("\n"):
            if line.startswith("#") or line.startswith("ERROR") or line == "\n" or not line:
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
                if i > num_embryos + 3:
                    object_type = "planetesimal"
                else:
                    object_type = "embryo"
            else:
                wmf = mmf = 0
                cmf = 1
                if columns[1] == 0:
                    object_type = "sun"
                else:
                    object_type = "gas giant"
            extradata.pdata[hash.value] = ParticleData(
                water_mass_fraction=wmf,
                type=object_type
            )

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
        assert sim.N == num_planetesimals + num_embryos + 3
        sim.move_to_com()
        extradata.meta.initial_N = sim.N
        extradata.meta.initial_N_planetesimal = num_planetesimals
        extradata.meta.initial_N_embryo = num_embryos
        extradata.energy.set_initial_energy(sim.calculate_energy())
        cputimeoffset = walltimeoffset = 0
        t = 0
    else:
        if fn.with_suffix(".lock").exists():
            raise FileExistsError("Lock file found, is the simulation currently running?")
        copy(fn.with_suffix(".bin"), fn.with_suffix(".bak.bin"))
        copy(fn.with_suffix(".extra.json"), fn.with_suffix(".extra.bak.json"))
        sa = SimulationArchive(str(fn.with_suffix(".bin")))
        extradata = ExtraData.load(fn.with_suffix(".extra.json"))
        tmax = extradata.meta.tmax
        per_savestep = extradata.meta.per_savestep
        t = extradata.meta.current_time - per_savestep
        sim = sa.getSimulation(t=t)
        sim.move_to_com()
        sim.ri_whfast.recalculate_coordinates_this_timestep = 1
        sim.integrator_synchronize()

        num_savesteps = extradata.meta.num_savesteps
        cputimeoffset = extradata.meta.cputime
        walltimeoffset = extradata.meta.walltime

    assert sim.dt < innermost_period(sim) / MIN_TIMESTEP_PER_ORBIT

    def collision_resolve_handler(sim_p: POINTER_REB_SIM, collision: reb_collision) -> int:
        return merge_particles(sim_p, collision, ed=extradata)

    sim.collision_resolve = collision_resolve_handler

    # show_orbits(sim)

    fn.with_suffix(".lock").touch()
    print("start")

    abort = False
    while t <= tmax:
        print()
        print(f"{t / tmax * 100:.2f}%")
        try:
            print(f"integrating until {t}")
            sim.integrate(t, exact_finish_time=0)
            print("dt", sim.dt)
            print("t", t)
            t += per_savestep
        except Escape:
            handle_escape(sim, extradata)
        except NoParticles:
            print("No Particles left")
            abort = True
        print("N", sim.N)
        print("N_active", sim.N_active)
        sim.simulationarchive_snapshot(str(fn.with_suffix(".bin")))
        extradata.meta.walltime = time.perf_counter() - start + walltimeoffset
        extradata.meta.cputime = time.process_time() + cputimeoffset
        extradata.meta.current_time = t
        # extradata.meta.current_steps = i
        extradata.energy.add_energy_value(sim.calculate_energy())
        print("total impulse", total_impulse(sim))
        extradata.save(fn.with_suffix(".extra.json"))
        assert sim.dt < innermost_period(sim) / MIN_TIMESTEP_PER_ORBIT
        print("fraction", innermost_period(sim) / MIN_TIMESTEP_PER_ORBIT)
        if abort:
            exit(1)
    print("finished")
    fn.with_suffix(".lock").unlink()


if __name__ == '__main__':
    fn = filename_from_argv()
    process_friendlyness(fn)
    try:
        main(fn)
    except KeyboardInterrupt:
        print("aborting")
        lockfile = fn.with_suffix(".lock")
        print(f"deleting {lockfile}")
        lockfile.unlink()
