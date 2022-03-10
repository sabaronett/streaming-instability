#!/usr/bin/env python3
#==============================================================================
# vpxd_hist.py
#
# Computes and saves the probability density function histogram and central
# moments of the zonal radial velocity distribution, weighted by particle
# densities, of the dust during the saturated, turbulent state of the SI.
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
from scipy.stats import skew, kurtosis

t_sat, n_bins = float(sys.argv[1]), int(sys.argv[2]) # t_sat in code unit [T]
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt'] # time between vp1 outputs
i_sat = int(t_sat/dt)          # output index of saturated state
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out1.*.athdf')))
sat_outputs = outputs[i_sat:]  # slice saturated state
vpxs, rhops = [], []

print(f'Compiling data...', flush=True)

for i,output in enumerate(sat_outputs):
    data = athena_read.athdf(output)
    vpxs.append(data['vp1'])
    rhops.append(data['rhop'])
    print('  {:.0%}'.format(i/len(sat_outputs)), flush=True)

print('  100%\nComputing histogram and central moments...', flush=True)
hist, bin_edges = np.histogram(vpxs, bins=n_bins, weights=rhops, density=True)
wavg = np.average(vpxs, weights=rhops)
wdist = vpxs*rhops
wvar = np.var(wdist)
wskew = skew(wdist, axis=None)
wkurt = kurtosis(wdist, axis=None)
print('... Done.\nSaving...', flush=True)
np.savez_compressed('output/vpxd_hist', hist=hist, bin_edges=bin_edges,
                    wavg=wavg, wvar=wvar, wskew=wskew, wkurt=wkurt)
