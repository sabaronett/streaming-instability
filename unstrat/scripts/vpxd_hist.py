#!/usr/bin/env python3
#==============================================================================
# vpxd_hist.py
#
# Computes and saves the probability density function histogram of zonal radial
# velocities, weighted by particle densities, of the dust during the saturated,
# turbulent state of the SI.
#
# Author: Stanley A. Baronett
# Created: 2022-03-09
# Last Modified: 2022-03-09
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

t_sat, n_bins = float(sys.argv[1]), int(sys.argv[2]) # t_sat in code unit [T]
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt'] # time between vp1 outputs
i_sat = int(t_sat/dt)          # output index of saturated state
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out1.*.athdf')))
sat_outputs = outputs[i_sat:]  # slice saturated state
vp1s, rhops = [], []

print(f'Compiling data...', flush=True)

for i,output in enumerate(sat_outputs):
    data = athena_read.athdf(output)
    vp1s = np.append(vp1s, data['vp1'].flatten())
    rhops = np.append(rhops, data['rhop'].flatten())

    print('  {:.0%}'.format(i/len(sat_outputs)), flush=True)

print('  100%\nComputing histogram...', flush=True)

hist, bin_edges = np.histogram(vp1s, bins=n_bins, weights=rhops, density=True)

print('... Done.\nSaving...', flush=True)
np.savez_compressed('output/vpxd_hist', hist=hist, bin_edges=bin_edges)
