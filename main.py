import json
from math import radians

import numpy as np
from rebound import Particle, Simulation
from reboundx import Extras

data = [
    [5.202887e0, 0.04838624e0, 1.30439695e0, 274.2545707e0, 100.4739091e0, 19.66796068e0, 0.000954792e0, "JUPITER"],
    [9.53667594e0, 0.05386179e0, 2.48599187e0, 338.9364538e0, 113.6624245e0, 317.3553659e0, 0.000285886e0, "SATURN"],
    [19.18916464e0, 0.04725744e0, 0.77263783e0, 96.93735127e0, 74.01692503e0, 142.2838282e0, 4.36624e-05, "URANUS"],
    [30.06992276e0, 0.00859048e0, 1.77004347e0, 273.1805365e0, 131.7842257e0, 259.915208e0, 5.15139e-05, "NEPTUNE"]
]


def main():
    sim = Simulation()
    rebx = Extras(sim)
    rebx.register_param("water", "REBX_TYPE_DOUBLE")
    sim.units = ('yr', 'AU', 'Msun')
    sim.add(m=1.0)  # Sun
    for planet in data:
        part = Particle(a=planet[0], e=planet[1], inc=radians(planet[2]), omega=radians(planet[3]),
                        Omega=radians(planet[4]), M=radians(planet[5]), m=planet[6], simulation=sim)
        sim.add(part)
    max_n = sim.N
    sim.particles[1].params["water"] = 0.3
    print("start")
    tmax = 1e5
    savesteps = 1000
    times = np.linspace(0., tmax, savesteps)
    sim.move_to_com()
    sim.automateSimulationArchive("sa.bin", interval=savesteps)
    for i, t in enumerate(times):
        sim.integrate(t)
        print(f"done: {i / 10}% ({sim.N}")
    with open("sa.meta.json", "w") as f:
        meta = {"tmax": tmax, "savesteps": savesteps, "max_n": max_n}
        json.dump(meta, f)


if __name__ == '__main__':
    main()
