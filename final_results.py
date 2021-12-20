import glob
from typing import List

import pandas as pd
from rebound import SimulationArchive, Simulation
from scipy.constants import mega

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, is_potentially_habitable, Particle, earth_mass, earth_water_mass

# pd.set_option('display.max_columns', None)
pd.options.display.max_columns = None
pd.options.display.width = 0

minmax=True

def chunks(lst, lst2):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), 2):
        yield lst[i:i + 2] + lst2[i:i + 2]


def print_row(mean, std, name: str):
    strings = [name]
    for value, pm in zip(list(mean), list(std)):
        strings.append(f"{value:.1f}\\pm {pm:.1f}")
    print(" & ".join(strings) + r" \\")


def get_masses(ed: ExtraData, planets: List[Particle]):
    total = 0
    water = 0
    for p in planets:
        data = ed.pd(p)
        # print(data.water_mass)
        total += data.total_mass / earth_mass
        water += data.water_mass / earth_water_mass
    return total, water


methods = {
    "rbf": "data/final_rbf_*.bin",
    "nn": "data/final_nn_*.bin",
    "lz": "data/final_lz_*.bin",
    "pm": "data/final_pm_*.bin",
}
columns = ["num_planets", "num_planets_pot", "M_planets", "M_planets_pot", "M_water", "M_water_pot", "escaped_mass",
           "escaped_water_mass", "sun_mass", "sun_water_mass", "gas_giant_mass", "gas_giant_water_mass",
           "collision_mass", "collision_water_mass","last_collision_time"]

maintable = []

for name, filepath in methods.items():
    table = []
    rows = []
    for file in sorted(glob.glob(filepath)):
        fn = filename_from_argv(file)
        # print(fn)

        ed = ExtraData.load(fn)
        if ed.meta.current_time < ed.meta.tmax:
            print("not yet finished")
            continue

        sa = SimulationArchive(str(fn.with_suffix(".bin")))
        last_sim: Simulation = sa[-1]
        planets = []
        for particle in last_sim.particles:
            if ed.pd(particle).type in ["sun", "gas giant"]:
                continue
            if ed.pd(particle).type == "planetesimal":
                continue
            # print(particle.r * astronomical_unit / earth_radius)
            planets.append(particle)

        pothab_planets = [p for p in planets if is_potentially_habitable(p)]
        num_planets = len(planets)
        num_planets_pot = len(pothab_planets)
        M_planets, M_water = get_masses(ed, planets)
        M_planets_pot, M_water_pot = get_masses(ed, pothab_planets)

        # mass of particles thrown into Sun/escaped

        escaped_mass = 0
        escaped_water_mass = 0
        sun_mass = 0
        sun_water_mass = 0
        for particle in ed.pdata.values():
            if particle.escaped:
                if particle.type in ["sun", "gas giant"]:
                    print(fn, particle, "escaped", particle.type)
                    continue
                escaped_mass += particle.total_mass / earth_mass
                escaped_water_mass += particle.water_mass / earth_water_mass
            elif particle.collided_with_sun:
                sun_mass += particle.total_mass / earth_mass
                sun_water_mass += particle.water_mass / earth_water_mass

        # count mass lost to gas giants
        gas_giant_mass = 0
        gas_giant_water_mass = 0
        for col_id, collision in ed.tree.get_tree().items():
            gas_parent = None
            other_parent = None
            gas_parent_counter = 0
            for parent in collision["parents"]:
                if ed.pdata[parent].type == "gas giant":
                    gas_parent_counter += 1
                    gas_parent = parent
                else:
                    other_parent = parent
            if gas_parent_counter > 1:
                print(f"it seems like {gas_parent_counter} gas giants collided with each other")
                continue
            if gas_parent:
                previous_body = ed.pdata[other_parent]
                gas_giant_mass += previous_body.total_mass / earth_mass
                gas_giant_water_mass += previous_body.water_mass / earth_water_mass

        # count mass lost in collisions
        collision_mass = 0
        collision_water_mass = 0
        for col_id, collision in ed.tree.get_tree().items():
            parent_mass = 0
            parent_water_mass = 0
            for parent in collision["parents"]:
                parent_body = ed.pdata[parent]
                parent_mass += parent_body.total_mass / earth_mass
                parent_water_mass += parent_body.water_mass / earth_water_mass
            child = ed.pdata[col_id]
            diff = parent_mass - child.total_mass / earth_mass
            diff_water = parent_water_mass - child.water_mass / earth_water_mass
            collision_mass += diff
            collision_water_mass += diff_water
        if collision_water_mass < 1e-10:
            collision_water_mass = 0
        if collision_mass < 1e-10:
            collision_mass = 0

        # last collision time

        last_collision_time = 0
        for col_id, collision in ed.tree.get_tree().items():
            meta: CollisionMeta = collision["meta"]
            if meta.time > last_collision_time:
                last_collision_time = meta.time
        last_collision_time /= mega

        values = [num_planets, num_planets_pot, M_planets, M_planets_pot, M_water, M_water_pot, escaped_mass,
                  escaped_water_mass, sun_mass, sun_water_mass, gas_giant_mass, gas_giant_water_mass,
                  collision_mass, collision_water_mass, last_collision_time]
        # print(values)
        table.append(values)
        rows.append(str(fn))

    pd.set_option('display.float_format', lambda x: '%.1f' % x)
    df = pd.DataFrame(table, index=rows, columns=columns)
    # print(df)
    # print("\n-----\n")
    # print_row(df.mean(), df.std(), name)
    if minmax:
        maintable.append(list(zip(df.min(), df.max())))
    else:
        maintable.append(list(zip(df.mean(), df.std())))

transposed = list(map(list, zip(*maintable)))

for row, name in zip(transposed, columns):
    strings = [name.replace("_", "\\_")]
    for value, pm in row:
        if minmax:
            strings.append(f"{value:.1f} -- {pm:.1f}")
        else:
            strings.append(f"{value:.1f}\\pm {pm:.1f}")
    print(" & ".join(strings) + r" \\")
