"""
The .energylog.csv files logged the energy every 100 years which resulted in huge files.
This script converts the data to a compressed hdf5 archive.
"""
import h5py

from utils import filename_from_argv

times = []
values = []
fn = filename_from_argv()
with fn.with_suffix(".energylog.csv").open() as f:
    for line in f:
        time_str, val_str = line.split(",")
        times.append(int(float(time_str)))
        values.append(float(val_str))
print("creating hdf5 file")
with h5py.File(fn.with_suffix(".energylog.hdf5"), "w") as f:
    times_dset = f.create_dataset("times", data=times, compression="gzip", shuffle=True, fletcher32=True)
    print(times_dset.dtype)
    values_dset = f.create_dataset("values", data=values, compression="gzip", shuffle=True, fletcher32=True)
    print(values_dset.dtype)
