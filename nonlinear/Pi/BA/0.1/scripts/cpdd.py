"""Cumulative Particle-Density Distribution.

Computes and outputs the time-averaged cumulative particle-density
distribution of the SI run during the saturated, turbulent state, including
time-varying minima, maxima, and standard deviations in logarithmic space.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

# Set earliest simulation time of saturated state
t_sat = 10                                   # / T

# Collect .athdf inputs, outputs, sim consts.
athinput = athena_read.athinput('../athinput.si')
Omega = athinput['problem']['omega']         # local Keplerian ang. freq.
epsilon = athinput['problem']['epsilon']     # avg. BG dust/gas ρ-ratio
dt = athinput['output2']['dt']               # time between rhop outputs
i_sat = int(t_sat / dt)                      # output index of sat. state
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                            '.out2.*.athdf')))
rhops = []                                   # sorted dust density snapshots

# Load & process saturated-state data into memory
sat_outputs = outputs[i_sat:]
for output in sat_outputs:
    data = athena_read.athdf(output)
    temp = data['rhop'].flatten() / epsilon  # flatten & convert
    rhops.append(np.sort(temp))              # sort

# Find min., max., avg. of each ordered rhop over saturated state
mins = np.amin(rhops, axis=0)
maxs = np.amax(rhops, axis=0)
avgs = np.average(rhops, axis=0)

# Find std. dev. computed in log space
lnrhops = np.log(rhops)
stds = np.exp(np.std(lnrhops, axis=0))
cdf = np.linspace(1, 0 , mins.size, endpoint=False)

np.savez_compressed('../output/cpdd', cdf=cdf, mins=mins,
                    maxs=maxs, avgs=avgs, stds=stds)
