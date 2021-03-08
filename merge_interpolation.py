import sys
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy.interpolate import Rbf

sys.path.append("./bac")

from bac.simulation_list import SimulationList
from bac.CustomScaler import CustomScaler

testrun = len(sys.argv) > 2 and sys.argv[2] == "test"

print("loading interpolation dataset")
simulations = SimulationList.jsonlines_load(Path("./rsmc_dataset.jsonl"))

scaler = CustomScaler()
scaler.fit(simulations.X)

scaled_data = scaler.transform_data(simulations.X)
output_data = np.array([simulations.Y_water, simulations.Y_mantle, simulations.Y_core])
if testrun:
    # keep memory usage low in tests
    output_data = output_data[::, :100]
    scaled_data = scaled_data[:100]
interpolator = Rbf(*scaled_data.T, output_data.T, function="linear", mode="N-D")
del simulations

print("finished loading interpolation dataset")


def interpolate(alpha, velocity, projectile_mass, gamma) -> Tuple[float, float, float]:
    hard_coded_water_mass_fraction = 1e-5  # workaround to get proper results for water poor collisions
    testinput = [alpha, velocity, projectile_mass, gamma,
                 hard_coded_water_mass_fraction, hard_coded_water_mass_fraction]

    print("# alpha velocity projectile_mass gamma target_water_fraction projectile_water_fraction\n")
    print(" ".join(map(str, testinput)))

    scaled_input = list(scaler.transform_parameters(testinput))
    water_retention, mantle_retention, core_retention = interpolator(*scaled_input)
    print(mantle_retention, core_retention)
    return float(water_retention), float(mantle_retention), float(core_retention)


if __name__ == '__main__':
    print(interpolate(32, 3, 7.6e22, 0.16))
