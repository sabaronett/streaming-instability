import numpy as np

runs = ['AB/0.05', 'AB/0.05-256', 'BA/0.05', 'BA/0.05-256']

for run in runs:
    path = '../nonlinear/Pi/%s/output/'%run
    data = np.loadtxt(path+'growth.txt')
    times = data[:, 0]
    rhopmax = data[:, 1]
    data = np.loadtxt(path+'cpdd_min.txt')
    cdf = data[:, 1]
    mins = data[:, 0]
    data = np.loadtxt(path+'cpdd_max.txt')
    maxs = data[:, 0]
    data = np.loadtxt(path+'cpdd_avg.txt')
    avgs = data[:, 0]
    data = np.loadtxt(path+'cpdd_std.txt')
    stds = data[:, 0]
    np.savez_compressed(path+'growth', times=times, rhopmax=rhopmax)
    np.savez_compressed(path+'cpdd', cdf=cdf, mins=mins, maxs=maxs, avgs=avgs,
                        stds=stds)
