import re
import time
from ctypes import Structure, c_uint32, c_double, c_uint, cdll, c_int, create_string_buffer, c_char_p
from dataclasses import dataclass
from pathlib import Path
from shutil import copy
from sys import argv

import rebound
import yaml
from rebound import Simulation, Particle, NoParticles, SimulationArchive
from rebound.simulation import POINTER_REB_SIM, reb_collision
from scipy.constants import astronomical_unit, mega, year

from extradata import ExtraData, ParticleData
from merge import merge_particles
from radius_utils import PlanetaryRadius
from utils import unique_hash, filename_from_argv, innermost_period, total_momentum, process_friendlyness, total_mass, \
    third_kepler_law, solar_radius, git_hash, check_heartbeat_needs_recompile

MIN_TIMESTEP_PER_ORBIT = 20

abort = False


class hb_event(Structure):
    _fields_ = [("hash", c_uint32),
                ("time", c_double),
                ("new", c_uint)]
    hash: int
    time: float
    new: int


hb_event_list = hb_event * 500


@dataclass
class Parameters:
    initcon_file: str
    perfect_merging: bool
    no_merging: bool = False


def main(fn: Path, testrun=False):
    global abort
    start = time.perf_counter()

    if not fn.with_suffix(".bin").exists():
        if not testrun:
            with open(fn.with_suffix(".yaml")) as f:
                parameters = Parameters(**yaml.safe_load(f))
        else:
            parameters = Parameters(perfect_merging=True, initcon_file="initcon/conditions_many.input")
        # set up a fresh simulation
        sim = Simulation()

        sim.units = ('yr', 'AU', 'kg')
        # sim.boundary = "open"
        # boxsize = 100
        # sim.configure_box(boxsize)
        sim.integrator = "mercurius"
        sim.dt = 1e-2
        sim.ri_ias15.min_dt = 0.0001 / 365
        if not parameters.no_merging:
            sim.collision = "direct"
        sim.ri_mercurius.hillfac = 3.
        sim.testparticle_type = 1
        tmax = 200 * mega
        num_savesteps = 20000
        if testrun:
            tmax /= 200000
            num_savesteps /= 1000
        per_savestep = tmax / num_savesteps
        extradata = ExtraData()
        # times = np.linspace(0., tmax, savesteps)
        extradata.meta.tmax = tmax
        extradata.meta.per_savestep = per_savestep
        extradata.meta.num_savesteps = num_savesteps
        extradata.meta.git_hash = git_hash()
        extradata.meta.git_hash = rebound.__githash__
        extradata.meta.perfect_merging = parameters.perfect_merging
        extradata.meta.initcon_file = parameters.initcon_file
        extradata.meta.no_merging = parameters.no_merging

        initcon = Path(parameters.initcon_file).read_text()
        num_embryos = int(re.search(r"Generated (\d+) minor bodies", initcon, re.MULTILINE).group(1))
        num_planetesimals = int(re.search(r"Generated (\d+) small bodies", initcon, re.MULTILINE).group(1))
        sim.N_active = num_embryos + 3
        i = 1
        for line in initcon.split("\n"):
            if line.startswith("#") or line.startswith("ERROR") or line == "\n" or not line:
                continue
            columns = list(map(float, line.split()))
            hash = unique_hash(extradata)
            if len(columns) > 7:
                cmf, mmf, wmf = columns[7:]
                total_fractions = cmf + mmf + wmf
                if total_fractions != 1:
                    diff = 1 - total_fractions
                    print(f"fractions don't add up by {diff}")
                    print("adding rest to mmf")
                    mmf += diff
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
                core_mass_fraction=cmf,
                type=object_type,
                total_mass=columns[0]
            )

            if columns[1] == 0:  # that should not be needed, but nevertheless is
                part = Particle(m=columns[0], hash=hash, r=solar_radius / astronomical_unit)
            else:
                part = Particle(
                    m=columns[0], a=columns[1], e=columns[2],
                    inc=columns[3], omega=columns[4],
                    Omega=columns[5], M=columns[6],
                    simulation=sim,
                    hash=hash,
                    r=PlanetaryRadius(columns[0], wmf, cmf).total_radius / astronomical_unit
                )
            sim.add(part)
            i += 1
        assert sim.N == num_planetesimals + num_embryos + 3
        sim.move_to_com()
        extradata.meta.initial_N = sim.N
        extradata.meta.initial_N_planetesimal = num_planetesimals
        extradata.meta.initial_N_embryo = num_embryos
        extradata.history.append(
            energy=sim.calculate_energy(),
            momentum=total_momentum(sim),
            total_mass=total_mass(sim),
            time=sim.t,
            N=sim.N,
            N_active=sim.N_active
        )
        cputimeoffset = walltimeoffset = 0
        t = 0
    else:
        if fn.with_suffix(".lock").exists():
            raise FileExistsError("Lock file found, is the simulation currently running?")
        copy(fn.with_suffix(".bin"), fn.with_suffix(".bak.bin"))
        copy(fn.with_suffix(".extra.json"), fn.with_suffix(".extra.bak.json"))
        sa = SimulationArchive(str(fn.with_suffix(".bin")))
        extradata = ExtraData.load(fn)
        tmax = extradata.meta.tmax
        per_savestep = extradata.meta.per_savestep
        sim = sa[-1]
        t = round(sim.t + per_savestep)
        print(f"continuing from {t}")
        sim.move_to_com()
        sim.ri_mercurius.recalculate_coordinates_this_timestep = 1
        sim.integrator_synchronize()

        if extradata.meta.git_hash != git_hash():
            print("warning: The saved output was originally run with another version of the code")
            print(f"original: {extradata.meta.git_hash}")
            print(f"current: {git_hash()}")
        num_savesteps = extradata.meta.num_savesteps
        cputimeoffset = extradata.meta.cputime
        walltimeoffset = extradata.meta.walltime

    check_heartbeat_needs_recompile()
    clibheartbeat = cdll.LoadLibrary("heartbeat/heartbeat.so")
    clibheartbeat.init_logfile.argtypes = [c_char_p]
    logfile = create_string_buffer(128)
    logfile.value = str(fn.with_suffix(".energylog.csv")).encode()
    clibheartbeat.init_logfile(logfile)
    sim.heartbeat = clibheartbeat.heartbeat
    innermost_semimajor_axis = third_kepler_law(
        orbital_period=sim.dt * year * MIN_TIMESTEP_PER_ORBIT
    ) / astronomical_unit
    print(f"innermost semimajor axis is {innermost_semimajor_axis}")

    c_double.in_dll(clibheartbeat, "min_distance_from_sun_squared").value = innermost_semimajor_axis ** 2
    c_double.in_dll(clibheartbeat, "max_distance_from_sun_squared").value = 30 ** 2

    assert sim.dt < innermost_period(sim) / MIN_TIMESTEP_PER_ORBIT

    def collision_resolve_handler(sim_p: POINTER_REB_SIM, collision: reb_collision) -> int:
        global abort  # needed as exceptions don't halt integration
        try:
            return merge_particles(sim_p, collision, ed=extradata)
        except BaseException as exception:
            print("exception during collision_resolve")
            print(exception)
            abort = True
            sim_p.contents._status = 1
            raise exception

    sim.collision_resolve = collision_resolve_handler

    # show_orbits(sim)

    fn.with_suffix(".lock").touch()
    print("start")

    while t <= tmax:
        print()
        print(f"{t / tmax * 100:.2f}%")
        try:
            print(f"integrating until {t}")
            sim.integrate(t, exact_finish_time=0)
            print("dt", sim.dt)
            print("t", t)
            t += per_savestep
        except NoParticles:
            print("No Particles left")
            abort = True
        print("N", sim.N)
        print("N_active", sim.N_active)

        print("fraction", innermost_period(sim) / MIN_TIMESTEP_PER_ORBIT)
        assert sim.dt < innermost_period(sim) / MIN_TIMESTEP_PER_ORBIT

        escape: hb_event
        sun_collision: hb_event
        for escape in hb_event_list.in_dll(clibheartbeat, "hb_escapes"):
            if not escape.new:
                continue
            print("escape:", escape.time, escape.hash)
            extradata.pdata[escape.hash].escaped = escape.time
            escape.new = 0  # make sure to not handle it again
        c_int.in_dll(clibheartbeat, "hb_escape_index").value = 0
        for sun_collision in hb_event_list.in_dll(clibheartbeat, "hb_sun_collisions"):
            if not sun_collision.new:
                continue
            print("sun collision:", sun_collision.time, sun_collision.hash)
            extradata.pdata[sun_collision.hash].collided_with_sun = sun_collision.time
            sun_collision.new = 0
        print(c_int.in_dll(clibheartbeat, "hb_sun_collision_index").value, "found")
        c_int.in_dll(clibheartbeat, "hb_sun_collision_index").value = 0
        sim.simulationarchive_snapshot(str(fn.with_suffix(".bin")))
        extradata.meta.walltime = time.perf_counter() - start + walltimeoffset
        extradata.meta.cputime = time.process_time() + cputimeoffset
        extradata.meta.current_time = t
        extradata.history.append(
            energy=sim.calculate_energy(),
            momentum=total_momentum(sim),
            total_mass=total_mass(sim),
            time=sim.t,
            N=sim.N,
            N_active=sim.N_active
        )
        extradata.save(fn)
        if abort:
            print("aborted")
            exit(1)
    print("finished")
    fn.with_suffix(".lock").unlink()


if __name__ == '__main__':
    fn = filename_from_argv()
    process_friendlyness(fn)
    testrun = False
    if len(argv) > 2 and argv[2] == "test":
        testrun = True
    try:
        main(fn, testrun)
    except KeyboardInterrupt:
        print("aborting")
        lockfile = fn.with_suffix(".lock")
        print(f"deleting {lockfile}")
        lockfile.unlink()
