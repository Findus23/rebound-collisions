import sys
from typing import List

import numpy as np
from numpy import linalg, sqrt
from rebound import Simulation, Particle
from scipy.constants import astronomical_unit, G, year

from extradata import ExtraData, ParticleData
from radius_utils import radius
from utils import unique_hash, clamp

sys.path.append("./bac")

from bac.simulation_list import SimulationList
from bac.CustomScaler import CustomScaler
from bac.interpolators.rbf import RbfInterpolator

simulations = SimulationList.jsonlines_load()

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


def get_mass_fractions(alpha, velocity_original, escape_velocity, gamma, projectile_mass, target_water_fraction,
                       projectile_water_fraction):
    velocity_si = velocity_original * astronomical_unit / year
    print("v_esc", escape_velocity)
    velocity = velocity_si / escape_velocity
    print("v_orig,v_si", velocity_original, velocity_si)
    print("v", velocity)
    if alpha > 90:
        alpha = 180 - alpha
    if gamma > 1:
        gamma = 1 / gamma
    alpha = clamp(alpha, 0, 60)
    orig_velocity = velocity
    velocity = clamp(velocity, 1, 5)

    m_ceres = 9.393e+20
    m_earth = 5.9722e+24
    projectile_mass = clamp(projectile_mass, 2 * m_ceres, 2 * m_earth)
    gamma = clamp(gamma, 1 / 10, 1)

    water_retention, mass_retention = interpolate(alpha, velocity, projectile_mass, gamma)

    water_retention = clamp(water_retention, 0, 1)
    mass_retention = clamp(mass_retention, 0, 1)

    metadata = {"water_retention": water_retention, "mass_retention": mass_retention,
                "testinput": [alpha, velocity, projectile_mass, gamma],
                "velocity_si": velocity_si, "escape_velocity": escape_velocity, "orig_velocity": orig_velocity}

    return water_retention, mass_retention, metadata


def merge_particles(sim: Simulation, ed: ExtraData):
    print("colliding")
    collided: List[Particle] = []
    p: Particle
    for p in sim.particles:
        # print(p.lastcollision, sim.t)
        # if p.lastcollision == sim.t:
        if p.lastcollision >= sim.t - 1:
            collided.append(p)
    # if not collided:
    #     print("empty collision")
    #     return
    print(collided)
    assert len(collided) == 2, "More or fewer than 2 objects collided with each other"
    cp1: Particle  # projectile
    cp2: Particle  # target
    cp1, cp2 = collided

    projectile_wmf = ed.pdata[cp1.hash.value].water_mass_fraction
    target_wmf = ed.pdata[cp2.hash.value].water_mass_fraction

    v1 = np.array(cp1.vxyz)
    v2 = np.array(cp2.vxyz)
    vdiff = linalg.norm(v2 - v1)  # AU/year
    v1_u = v1 / linalg.norm(v1)
    v2_u = v2 / linalg.norm(v2)
    # https://stackoverflow.com/a/13849249/4398037
    ang = np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

    gamma = cp1.m / cp2.m

    escape_velocity = sqrt(2 * G * (cp1.m + cp2.m) / ((cp1.r + cp2.r) * astronomical_unit))

    print("interpolating")
    water_ret, stone_ret, meta = get_mass_fractions(
        alpha=ang, velocity_original=vdiff, escape_velocity=escape_velocity, gamma=gamma, projectile_mass=cp1.m,
        target_water_fraction=target_wmf, projectile_water_fraction=projectile_wmf)
    print(water_ret, stone_ret)
    hash = unique_hash()

    water_mass = cp1.m * projectile_wmf + cp2.m * target_wmf
    stone_mass = cp1.m + cp2.m - water_mass

    water_mass *= water_ret
    stone_ret *= stone_ret

    total_mass = water_mass + stone_mass
    final_wmf = water_mass / total_mass
    print(final_wmf)

    merged_planet = (cp1 * cp1.m + cp2 * cp2.m) / (cp1.m + cp2.m)
    merged_planet.m = total_mass
    merged_planet.hash = hash

    merged_planet.r = radius(merged_planet.m, final_wmf) / astronomical_unit
    ed.pdata[hash.value] = ParticleData(water_mass_fraction=final_wmf)

    meta["total_mass"] = total_mass
    meta["final_wmf"] = final_wmf
    meta["final_radius"] = merged_planet.r
    meta["target_wmf"] = target_wmf
    meta["projectile_wmf"] = projectile_wmf
    ed.tree.add(cp1, cp2, merged_planet, meta)

    cp1_hash = cp1.hash
    cp2_hash = cp2.hash

    # don't use cp1 and cp2 from now on as they will change

    print("removing", cp1_hash.value, cp2_hash.value)
    sim.remove(hash=cp1_hash)
    sim.remove(hash=cp2_hash)
    sim.add(merged_planet)
