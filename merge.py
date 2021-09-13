from copy import copy
from pprint import pprint
from typing import Tuple, Optional

import numpy as np
from numpy import linalg, sqrt
from rebound import Simulation, Particle
from rebound.simulation import POINTER_REB_SIM, reb_collision
from scipy.constants import astronomical_unit, G

from extradata import ExtraData, ParticleData, CollisionMeta, Input
from massloss import RbfMassloss, Massloss, LeiZhouMassloss, SimpleNNMassloss
from massloss.perfect_merging import PerfectMerging
from utils import unique_hash, clamp, PlanetaryRadius

massloss_estimator: Optional[Massloss] = None  # global waterloss estimator cache


def get_mass_fractions(input_data: Input) -> Tuple[float, float, float, CollisionMeta]:
    global massloss_estimator
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

    water_retention, mantle_retention, core_retention = \
        massloss_estimator.estimate(data.alpha, data.velocity_esc, data.projectile_mass, data.gamma, )

    metadata = CollisionMeta()
    metadata.interpolation_input = [data.alpha, data.velocity_esc, data.projectile_mass, data.gamma]
    metadata.input = input_data
    metadata.adjusted_input = data
    metadata.raw_water_retention = water_retention
    metadata.raw_mantle_retention = mantle_retention
    metadata.raw_core_retention = core_retention

    water_retention = clamp(water_retention, 0, 1)
    mantle_retention = clamp(mantle_retention, 0, 1)
    core_retention = clamp(core_retention, 0, 1)

    metadata.water_retention = water_retention
    metadata.mantle_retention = mantle_retention
    metadata.core_retention = core_retention

    return water_retention, mantle_retention, core_retention, metadata


def merge_particles(sim_p: POINTER_REB_SIM, collision: reb_collision, ed: ExtraData):
    global massloss_estimator
    print("--------------")
    print("colliding")
    sim: Simulation = sim_p.contents
    print("current time step", sim.dt)
    print(f"p1 is {collision.p1}")
    print(f"p2 is {collision.p2}")
    # the assignment to cp1 or cp2 seems to be random
    # also look at a copy instead of the original particles
    # to avoid issues after they have been modified
    cp1: Particle = sim.particles[collision.p1].copy()
    cp2: Particle = sim.particles[collision.p2].copy()

    # just calling the more massive one the target to keep its type/name
    # Sun<->Protoplanet -> Sun
    # and to keep collsions mostly reproducable
    if cp1.m > cp2.m:
        target = cp1
        projectile = cp2
    else:  # also when masses are the same
        target = cp2
        projectile = cp1

    if collision.p1 > collision.p2:
        lower_index_particle_index = collision.p2
    else:
        lower_index_particle_index = collision.p1

    print(f"colliding {target.hash.value} ({ed.pd(target).type}) "
          f"with {projectile.hash.value} ({ed.pd(projectile).type})")

    projectile_wmf = ed.pd(projectile).water_mass_fraction
    projectile_cmf = ed.pd(projectile).core_mass_fraction
    target_wmf = ed.pd(target).water_mass_fraction
    target_cmf = ed.pd(target).core_mass_fraction

    # get the velocities, velocity differences and unit vector as numpy arrays
    # all units are in sytem units (so AU/year)
    v1 = np.array(target.vxyz)
    v2 = np.array(projectile.vxyz)
    r1 = np.array(target.xyz)
    r2 = np.array(projectile.xyz)
    vdiff = v2 - v1
    rdiff = r2 - r1
    vdiff_n = linalg.norm(vdiff)
    rdiff_n = linalg.norm(rdiff)
    print("dt", sim.dt)
    # during a collision ias15 should always be used, otherwise something weird has happend
    assert sim.ri_mercurius.mode == 1

    print("rdiff", rdiff)
    print("vdiff", vdiff)
    print("sum_radii", target.r + projectile.r)
    print("rdiff_n", rdiff_n)
    print("vdiff_n", vdiff_n)
    ang = float(np.degrees(np.arccos(np.dot(rdiff, vdiff) / (rdiff_n * vdiff_n))))
    if ang > 90:
        ang = 180 - ang

    print("angle_deg", ang)
    print()
    # get mass fraction
    gamma = projectile.m / target.m

    # calculate mutual escape velocity (for norming the velocities in the interpolation) in SI units
    escape_velocity = sqrt(2 * G * (target.m + projectile.m) / ((target.r + projectile.r) * astronomical_unit))

    print("interpolating")

    if not massloss_estimator:
        methods = [RbfMassloss, LeiZhouMassloss, PerfectMerging, SimpleNNMassloss]
        per_name = {}
        for method in methods:
            per_name[method.name] = method
        try:
            estimator_class = per_name[ed.meta.massloss_method]
        except KeyError:
            print("invalid mass loss estimation method")
            print("please use one of these:")
            print(per_name)
            raise
        massloss_estimator = estimator_class()

    # let interpolation calculate water and mass retention fraction
    # meta is just a bunch of intermediate results that will be logged to help
    # understand the collisions better
    input_data = Input(
        alpha=ang,
        velocity_original=vdiff_n,
        escape_velocity=escape_velocity,
        gamma=gamma,
        projectile_mass=projectile.m,
        target_water_fraction=target_wmf,
        projectile_water_fraction=projectile_wmf,
    )

    water_ret, mantle_ret, core_ret, meta = get_mass_fractions(input_data)
    print("mass retentions:", water_ret, mantle_ret, core_ret)

    meta.collision_velocities = (v1.tolist(), v2.tolist())
    meta.collision_positions = (target.xyz, projectile.xyz)
    meta.collision_radii = (target.r, projectile.r)

    hash = unique_hash(ed)  # hash for newly created particle

    # handle loss of water and core mass
    water_mass = target.m * target_wmf + projectile.m * projectile_wmf
    core_mass = target.m * target_cmf + projectile.m * projectile_cmf
    mantle_mass = target.m + projectile.m - water_mass - core_mass

    water_mass *= water_ret
    mantle_mass *= mantle_ret
    core_mass *= core_ret

    total_mass = water_mass + mantle_mass + core_mass
    final_wmf = water_mass / total_mass
    final_cmf = core_mass / total_mass
    print(final_wmf)
    # create new object preserving momentum
    merged_planet = (target * target.m + projectile * projectile.m) / total_mass
    merged_planet.m = total_mass
    merged_planet.hash = hash

    merged_planet.r = PlanetaryRadius(merged_planet.m, final_wmf, final_cmf).total_radius / astronomical_unit
    ed.pdata[hash.value] = ParticleData(
        water_mass_fraction=final_wmf,
        core_mass_fraction=final_cmf,
        type=ed.pd(target).type,
        total_mass=total_mass
    )

    meta.final_wmf = final_wmf
    meta.final_radius = merged_planet.r
    meta.target_wmf = target_wmf
    meta.projectile_wmf = projectile_wmf
    meta.time = sim.t
    pprint(meta)

    ed.tree.add(target, projectile, merged_planet, meta)

    sim.particles[lower_index_particle_index] = merged_planet

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
    # always keep lower index particle and delete other one
    # this keeps the N_active working
    if lower_index_particle_index == collision.p1:
        print("deleting p2")
        return 2
    elif lower_index_particle_index == collision.p2:
        print("deleting p1")
        return 1
    else:
        raise ValueError("invalid index")
