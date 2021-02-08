import sys
from copy import copy
from pathlib import Path
from pprint import pprint
from typing import Tuple

import numpy as np
from numpy import linalg, sqrt
from rebound import Simulation, Particle, reb_simulation_integrator_mercurius
from rebound.simulation import POINTER_REB_SIM, reb_collision
from scipy.constants import astronomical_unit, G

from extradata import ExtraData, ParticleData, CollisionMeta, Input
from radius_utils import radius
from utils import unique_hash, clamp

sys.path.append("./bac")

from bac.simulation_list import SimulationList
from bac.CustomScaler import CustomScaler
from bac.interpolators.rbf import RbfInterpolator

simulations = SimulationList.jsonlines_load(Path("./save.jsonl"))

scaler = CustomScaler()
scaler.fit(simulations.X)

scaled_data = scaler.transform_data(simulations.X)
water_interpolator = RbfInterpolator(scaled_data, simulations.Y_water)
mass_interpolator = RbfInterpolator(scaled_data, simulations.Y_mass)


def interpolate(alpha, velocity, projectile_mass, gamma):
    hard_coded_water_mass_fraction = 0.15  # workaround to get proper results for water poor collisions
    testinput = [alpha, velocity, projectile_mass, gamma,
                 hard_coded_water_mass_fraction, hard_coded_water_mass_fraction]

    print("# alpha velocity projectile_mass gamma target_water_fraction projectile_water_fraction\n")
    print(" ".join(map(str, testinput)))

    scaled_input = list(scaler.transform_parameters(testinput))
    water_retention = water_interpolator.interpolate(*scaled_input)
    mass_retention = mass_interpolator.interpolate(*scaled_input)
    return float(water_retention), float(mass_retention)


def get_mass_fractions(input_data: Input) -> Tuple[float, float, CollisionMeta]:
    print("v_esc", input_data.escape_velocity)
    print("v_orig,v_si", input_data.velocity_original, input_data.velocity_si)
    print("v/v_esc", input_data.velocity_esc)
    data = copy(input_data)
    if data.gamma > 1:
        data.gamma = 1 / data.gamma
    data.alpha = clamp(data.alpha, 0, 60)
    data.velocity_esc = clamp(data.velocity_esc, 1, 5)

    m_ceres = 9.393e+20
    m_earth = 5.9722e+24
    data.projectile_mass = clamp(data.projectile_mass, 2 * m_ceres, 2 * m_earth)
    data.gamma = clamp(data.gamma, 1 / 10, 1)

    water_retention, mass_retention = interpolate(data.alpha, data.velocity_esc, data.projectile_mass, data.gamma)

    metadata = CollisionMeta()
    metadata.interpolation_input = [data.alpha, data.velocity_esc, data.projectile_mass, data.gamma]
    metadata.input = input_data
    metadata.adjusted_input = data
    metadata.raw_water_retention = water_retention
    metadata.raw_mass_retention = mass_retention

    water_retention = clamp(water_retention, 0, 1)
    mass_retention = clamp(mass_retention, 0, 1)

    metadata.water_retention = water_retention
    metadata.mass_retention = mass_retention

    return water_retention, mass_retention, metadata


def merge_particles(sim_p: POINTER_REB_SIM, collision: reb_collision, ed: ExtraData):
    print("--------------")
    print("colliding")
    sim: Simulation = sim_p.contents
    print("current time step", sim.dt, )
    print("mode", sim.ri_mercurius.mode)

    # the assignment to cp1 or cp2 is mostly random
    # (cp1 is the one with a lower index in sim.particles)
    # naming them projectile or target is therefore also arbitrary
    # also look at a copy instead of the original particles
    # to avoid issues after they have been modified
    cp1: Particle = sim.particles[collision.p1].copy()  # projectile
    cp2: Particle = sim.particles[collision.p2].copy()  # target

    # just calling the more massive one the main particle to keep its type/name
    # Sun<->Protoplanet -> Sun
    if cp1.m > cp2.m:
        main_particle_id = collision.p1
        main_particle = cp1
    else:
        main_particle_id = collision.p2
        main_particle = cp2

    print(f"colliding {ed.pd(cp1).type} with {ed.pd(cp2).type}")

    projectile_wmf = ed.pd(cp1).water_mass_fraction
    target_wmf = ed.pd(cp2).water_mass_fraction

    # get the velocities, velocity differences and unit vector as numpy arrays
    # all units are in sytem units (so AU/year)
    v1 = np.array(cp1.vxyz)
    v2 = np.array(cp2.vxyz)
    r1 = np.array(cp1.xyz)
    r2 = np.array(cp2.xyz)
    vdiff = v2 - v1
    rdiff = r2 - r1
    vdiff_n = linalg.norm(vdiff)
    rdiff_n = linalg.norm(rdiff)
    print("dt", sim.dt)
    merc: reb_simulation_integrator_mercurius = sim.ri_mercurius
    print("current mode", "ias15" if merc.mode else "whfast")
    # during a collision ias15 should always be used, otherwise something weird has happend
    assert merc.mode == 1

    print("rdiff", rdiff)
    print("vdiff", vdiff)
    print("sum_radii", cp1.r + cp2.r)
    print("rdiff_n", rdiff_n)
    print("vdiff_n", vdiff_n)
    ang = float(np.degrees(np.arccos(np.dot(rdiff, vdiff) / (rdiff_n * vdiff_n))))
    if ang > 90:
        ang = 180 - ang

    print("angle_deg", ang)
    print()
    # get mass fraction
    # if it is >1 it will be inverted during interpolation
    gamma = cp1.m / cp2.m

    # calculate mutual escape velocity (for norming the velocities in the interpolation) in SI units
    escape_velocity = sqrt(2 * G * (cp1.m + cp2.m) / ((cp1.r + cp2.r) * astronomical_unit))

    print("interpolating")

    # let interpolation calculate water and mass retention fraction
    # meta is just a bunch of intermediate results that will be logged to help
    # understand the collisions better
    input_data = Input(
        alpha=ang,
        velocity_original=vdiff_n,
        escape_velocity=escape_velocity,
        gamma=gamma,
        projectile_mass=cp1.m,
        target_water_fraction=target_wmf,
        projectile_water_fraction=projectile_wmf,
    )

    water_ret, stone_ret, meta = get_mass_fractions(input_data)
    print("interpolation finished")
    print(water_ret, stone_ret)

    meta.collision_velocities = (v1.tolist(), v2.tolist())
    meta.collision_positions = (cp1.xyz, cp2.xyz)
    meta.collision_radii = (cp1.r, cp2.r)

    hash = unique_hash()  # hash for newly created particle

    # handle loss of water and core mass
    water_mass = cp1.m * projectile_wmf + cp2.m * target_wmf
    stone_mass = cp1.m + cp2.m - water_mass

    water_mass *= water_ret
    stone_mass *= stone_ret

    total_mass = water_mass + stone_mass
    final_wmf = water_mass / total_mass
    print(final_wmf)
    # create new object preserving momentum
    merged_planet = (cp1 * cp1.m + cp2 * cp2.m) / total_mass
    merged_planet.m = total_mass
    merged_planet.hash = hash

    merged_planet.r = radius(merged_planet.m, final_wmf) / astronomical_unit
    ed.pdata[hash.value] = ParticleData(
        water_mass_fraction=final_wmf,
        type=ed.pd(main_particle).type,
        total_mass=total_mass
    )

    meta.final_wmf = final_wmf
    meta.final_radius = merged_planet.r
    meta.target_wmf = target_wmf
    meta.projectile_wmf = projectile_wmf
    meta.time = sim.t
    pprint(meta)

    ed.tree.add(cp1, cp2, merged_planet, meta)

    sim.particles[main_particle_id] = merged_planet

    sim.move_to_com()
    sim.integrator_synchronize()
    sim.ri_mercurius.recalculate_coordinates_this_timestep = 1
    sim.ri_mercurius.recalculate_dcrit_this_timestep = 1

    print("collision finished")
    print("--------------")
    # from rebound docs:
    # A return value of 0 indicates that both particles remain in the simulation.
    # A return value of 1 (2) indicates that particle 1 (2) should be removed from the simulation.
    # A return value of 3 indicates that both particles should be removed from the simulation.

    if main_particle_id == collision.p1:
        return 2
    else:
        return 1
