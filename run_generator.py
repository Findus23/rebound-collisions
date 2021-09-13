from pathlib import Path

import yaml

from water_sim import Parameters

basename = "final_rbf_"
outdir = Path("data")
for num in range(1, 11):
    param = Parameters(
        massloss_method="rbf",
        initcon_file=f"initcon/conditions_final{num}.input"
    )
    print(yaml.dump(param.__dict__))
    outfile = outdir / f"{basename}{num}.yaml"
    if outfile.exists():
        print("file exists", outfile)
        continue
    with outfile.open("w") as f:
        yaml.dump(param.__dict__, f)
