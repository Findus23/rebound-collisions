import numpy as np


def reorder(data, water_lost_median: float, std=False):
    return [
        *data[:9], water_lost_median,
        np.nan, np.nan,
        0 if std else data[9], np.nan,
        data[13]
    ]


def get_cb_data(pm=False):
    ncols = 14
    data = np.loadtxt("cb_data2.txt" if pm else "cb_data1.txt", usecols=range(1, ncols + 1))
    means = data.mean(axis=0).tolist()
    stds = data.std(axis=0).tolist()

    print(means)
    water_lost_median = 210 if pm else 257
    return reorder(means, water_lost_median=water_lost_median), \
           reorder(stds, water_lost_median=water_lost_median, std=True)
