"""Maximum particle density evolution.

Collects and outputs the maximum particle density as a function of time.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

# Collect .athdf outputs, init sim consts.
athinput = athena_read.athinput('../athinput.si')
Omega = athinput['problem']['omega']         # local Keplerian ang. freq.
T = 2*np.pi/Omega                            # orbital period
outputs = sorted(list(Path('../athdf').glob(athinput['job']['problem_id'] +
                                            '.out2.*.athdf')))
times, rhopmax = [], []                      # sim out times, max dust dens.

for output in outputs:                       # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time'] / T)
    rhopmax.append(np.amax(data['rhop']))
    print('{:3.1f}% done.'.format(100*len(rhopmax)/len(outputs)))

np.savez_compressed('../output/growth', times=np.asarray(times),
                    rhopmax=np.asarray(rhopmax))
