from math import isclose
from os.path import expanduser
from statistics import mean
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rebound import SimulationArchive, Simulation
from scipy.constants import mega

from extradata import ExtraData, CollisionMeta
from utils import filename_from_argv, is_potentially_habitable, Particle, earth_mass, earth_water_mass, \
    habitable_zone_inner, habitable_zone_outer, get_water_cmap, create_figure, add_au_e_label, \
    inner_solar_system_data, is_ci, get_cb_data, initial_saturn_a, initial_jupiter_a

# pd.set_option('display.max_columns', None)
pd.options.display.max_columns = None
pd.options.display.width = 0

minmax = True
plot = True


def chunks(lst, lst2):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), 2):
        yield lst[i:i + 2] + lst2[i:i + 2]


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
    "rbf": "data/final_rbf_NUM.bin",
    "rbf_sm": "data/final_rbf_NUM.bin",
    "nn": "data/final_nn_NUM.bin",
    "lz": "data/final_lz_correct_NUM.bin",
    "pm": "data/final_pm_NUM.bin",
}
columns = [
    "$N_\\text{planets}$ [1]", "$N_\\text{planets,pot}$ [1]", "$M_\\text{planets}$ [$M_\\oplus$]",
    "$M_\\text{planets,pot}$ [$M_\\oplus$]", "$M_\\text{water}$ [$M_{w,\\oplus}$]",
    "$M_\\text{water,pot}$ [$M_{w,\\oplus}$]", "$M_\\text{esc}$ [$M_\\oplus$]",
    "$M_\\text{esc,water}$ [$M_{w,\\oplus}$]", "$M_\\text{sun}$ [$M_\\oplus$]",
    "$M_\\text{sun,water}$ [$M_{w,\\oplus}$]", "$M_\\text{gas-giant}$  [$M_\\oplus$]",
    "$M_\\text{gas-giant,water}$ [$M_{w,\\oplus}$]", "$M_\\text{col}$ [$M_\\oplus$]",
    "$M_\\text{col,water}$ [$M_{w,\\oplus}$]", "$t_\\text{last-col}$ [Myr]"
]

maintable = []

final_bodies = {}

for name, filepath in methods.items():
    table = []
    rows = []
    fin_as = []
    fin_es = []
    fin_mass = []
    fin_core_mass = []
    fin_wmf = []
    max_num = 41 if name == "rbf" else 21
    for i in range(1, max_num):
        fn = filename_from_argv(filepath.replace("NUM", str(i)))
        print(fn)
        try:
            ed = ExtraData.load(fn)
        except FileNotFoundError as e:
            print(e.filename)
            continue
        if ed.meta.current_time < ed.meta.tmax:
            print("not yet finished")
            continue

        sa = SimulationArchive(str(fn.with_suffix(".bin")))
        last_sim: Simulation = sa[-1]
        planets = []
        a_values_per_planet = {}
        e_values_per_planet = {}
        for sim in sa:
            if sim.t < ed.meta.tmax - 10 * mega:
                continue
            for particle in sim.particles[1:]:
                hash = particle.hash.value
                orbit = particle.calculate_orbit()
                if hash not in a_values_per_planet:
                    a_values_per_planet[hash] = []
                    e_values_per_planet[hash] = []
                a_values_per_planet[hash].append(orbit.a)
                e_values_per_planet[hash].append(orbit.e)
        gas_giants = []
        for particle in last_sim.particles:
            particle_data = ed.pd(particle)
            if particle_data.type in ["sun", "planetesimal"]:
                continue
            if particle_data.type == "gas giant":
                gas_giants.append(particle)
                continue
            # print(particle.r * astronomical_unit / earth_radius)
            planets.append(particle)
            fin_as.append(mean(a_values_per_planet[particle.hash.value]))
            fin_es.append(mean(e_values_per_planet[particle.hash.value]))
            fin_mass.append(particle_data.total_mass)
            fin_core_mass.append(particle_data.total_mass * particle_data.core_mass_fraction)
            fin_wmf.append(particle_data.water_mass_fraction)

        if len(gas_giants) != 2:
            print(f"it seems like gas giants got lost ({len(gas_giants)} left)")
        elif not (
                isclose(gas_giants[0].a, initial_saturn_a, rel_tol=0.05)
                and
                isclose(gas_giants[1].a, initial_jupiter_a, rel_tol=0.05)
        ):
            print("gas giants moved")

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
                print(f"it seems like {gas_parent_counter} gas giants collided with each other in {col_id}")
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

        values = [num_planets, num_planets_pot, M_planets, M_planets_pot, M_water, M_water_pot,
                  sun_mass, sun_water_mass, escaped_mass, escaped_water_mass,
                  gas_giant_mass, gas_giant_water_mass,
                  collision_mass, collision_water_mass, last_collision_time]
        # print(values)
        table.append(values)
        rows.append(str(fn))

    if plot:
        mean_mass = 5.208403167890638e+24  # constant between plots
        size_mult = 100
        fin_mass = np.array(fin_mass)
        fin_core_mass = np.array(fin_core_mass)
        print("mean", mean_mass)
        sizes = (np.array(fin_mass) / mean_mass) ** (2 / 3) * size_mult
        core_sizes = (np.array(fin_core_mass) / mean_mass) ** (2 / 3) * size_mult
        with np.errstate(divide='ignore'):  # allow 0 water (becomes -inf)
            color_val = (np.log10(fin_wmf) + 5) / 5
        cmap = get_water_cmap()
        print(color_val)
        print(np.log10(fin_wmf))
        colors = cmap(color_val)
        fig, ax = create_figure()
        ax.scatter(fin_as, fin_es, s=sizes, zorder=10, cmap=(), c=colors)
        ax.scatter(fin_as, fin_es, s=core_sizes, zorder=12, color="black")
        for pname, planet in inner_solar_system_data.items():
            if pname == "earth":
                earth_wmf = earth_water_mass / earth_mass
                earth_color_val = (np.log10(earth_wmf) + 5) / 5
                fill_color = cmap(earth_color_val)
            else:
                fill_color = "white"
            ax.scatter(planet.a_au, planet.e, s=(planet.mass / mean_mass) ** (2 / 3) * size_mult,
                       color=fill_color, edgecolors="black", zorder=5)
        ax.axvspan(habitable_zone_inner, habitable_zone_outer, color="#eee")
        # add_water_colormap(fig, ax, cmap=cmap)
        add_au_e_label(ax)
        ax.set_xlim(0.25, 3.6)
        ax.set_ylim(-0.05, 0.55)
        fig.tight_layout()
        if not is_ci():
            plt.savefig(expanduser(f"~/tmp/final_bodies_{name}.pdf"))
        # plt.show()
    pd.set_option('display.float_format', lambda x: '%.1f' % x)
    df = pd.DataFrame(table, index=rows, columns=columns)
    print([a[2] for a in table])
    # print(df)
    # print("\n-----\n")
    # print_row(df.mean(), df.std(), name)
    if minmax:
        maintable.append(list(zip(df.min(), df.max())))
    else:
        maintable.append(list(zip(df.mean(), df.std())))

maintable.append(list(zip(*get_cb_data(minmax))))
maintable.append(list(zip(*get_cb_data(minmax, pm=True))))

transposed = list(map(list, zip(*maintable)))

for row, name in zip(transposed, columns):
    n, unit = name.split()

    strings = [" ".join([n, unit])]
    for value, pm in row:
        if np.isnan(value):
            strings.append("{-}")
        elif minmax:  # in that case (value,pm) = (min,max)
            strings.append(f"{value:.1f} -- {pm:.1f}")
        else:
            strings.append(f"{value:.1f}\\pm {pm:.1f}")
    print(" & ".join(strings) + r" \\")
