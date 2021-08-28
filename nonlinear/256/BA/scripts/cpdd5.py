"""Cumulative Particle-Density Distribution.

Computes and outputs the cumulative particle-density distribution of the SI
run over all dust density snapshots during the saturated, turbulent state.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

def makesubdir(name):
    if not os.path.exists(name):
        os.makedirs(name)

# Collect .athdf inputs, outputs; init sim consts. & grid
athinput = athena_read.athinput('../athinput.si')
Omega = athinput['problem']['omega']         # local Keplerian ang. freq.
tau_s = athinput['particles']['taus0']*Omega # dimensionless stopping time
epsilon = athinput['problem']['epsilon']     # avg. BG dust/gas ρ-ratio
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                            '.out2.*.athdf')))
rhopslist = []                               # dust densities

# Load & process saturated-state data into memory
i_sat = 375                                  # 1st index of sat. state
sat_outputs = outputs[i_sat:]
for output in sat_outputs:
    data = athena_read.athdf(output)
    rhopslist.append(data['rhop'] / epsilon) # convert to ρ_p / ρ_p0
rhops = np.asarray(rhopslist)                # no averaging
flat = np.ravel(rhops)
rhops = np.sort(flat)
# define prob. w/ total # cells & total # snapshots
cdf = np.linspace(1, 0 , rhops.size, endpoint=False)

# Save CPDD plotting data
makesubdir('../plots') # create file output directory
np.savez_compressed('../plots/cpdd', rhops=rhops, cdf=cdf)
