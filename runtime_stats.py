from pathlib import Path

import numpy as np
from scipy.constants import hour

from extradata import ExtraData
from utils import filename_from_argv, mode_from_fn

files = Path("data/").glob("final*.bin")
times = []
for file in files:
    fn = filename_from_argv(str(file))
    mode = mode_from_fn(fn)
    if mode == "lz":
        continue
    if mode != "lz_correct":
        continue
    print(fn)
    ed = ExtraData.load(fn)
    if str(fn) == "data/final_lz_correct_15":
        print(ed.meta.cputime / hour)
    times.append(ed.meta.cputime)

times = np.asarray(times)
print(times.mean() / hour)
print(times.std() / hour)
