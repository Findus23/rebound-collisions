from sys import argv
from typing import List

import pandas as pd
from rebound import SimulationArchive, Simulation

from extradata import ExtraData
from utils import filename_from_argv, is_potentially_habitable, Particle, earth_mass, earth_water_mass

# pd.set_option('display.max_columns', None)
pd.options.display.max_columns = None
pd.options.display.width = 0


def get_masses(ed: ExtraData, planets: List[Particle]):
    total = 0
    water = 0
    for p in planets:
        data = ed.pd(p)
        print(data.water_mass)
        total += data.total_mass / earth_mass
        water += data.water_mass / earth_water_mass
    return total, water


table = []
rows = []
for file in argv[1:]:
    fn = filename_from_argv(file)
    print(fn)

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
    num_pothab_planets = len(pothab_planets)
    M_planets, M_water = get_masses(ed, planets)
    M_planets_pot, M_water_pot = get_masses(ed, pothab_planets)

    values = [num_planets, num_pothab_planets, M_planets, M_planets_pot, M_water, M_water_pot]
    print(values)
    table.append(values)
    rows.append(str(fn))

print("\n-----\n")
pd.set_option('display.float_format', lambda x: '%.1f' % x)
df = pd.DataFrame(table,
                  index=rows,
                  columns=["num_planets", "num_pothab_planets", "M_planets", "M_planets_pot", "M_water", "M_water_pot"])
print(df)
print("\n-----\n")


def chunks(lst):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), 2):
        yield lst[i:i + 2]


def print_row(row, name: str):
    strings = [name]
    for n, pn in chunks(row):
        if int(n) == n:
            strings.append(f"{n:.0f} ({pn:.0f})")
        else:
            strings.append(f"{n:.1f} ({pn:.1f})")
    print(" & ".join(strings) + r" \\")


print_row(df.mean(), "mean")
if len(df) > 1:
    print_row(df.std(), "std")
print_row(df.min(), "min")
print_row(df.max(), "max")
print(df.describe())
