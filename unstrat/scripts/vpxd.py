#!/usr/bin/env python3
#==============================================================================
# vpxd.py
#
# Computes and saves the probability density function histogram and central
# moments of the radial velocity distribution of dust particles during the
# saturated, turbulent state of the SI.
#
# Author: Stanley A. Baronett
# Created: 2022-03-09
# Updated: 2022-07-08
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path
from scipy.stats import skew, kurtosis

# Collect Athena++ inputs, outputs, and sim constants
t_sat, bins = float(sys.argv[1]), int(sys.argv[2]) # t_sat in code unit [T]
athinput = athena_read.athinput('athinput.si')
Pi = athinput['problem']['duy0']                   # radial pressure gradient
etav_K = Pi*athinput['hydro']['iso_sound_speed']   # scaling factor
dt = athinput['outputp']['dt']                     # time between vpx outputs
i_sat = int(t_sat/dt)                              # sat state output index
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out1.*.athdf')))
sat_outputs = outputs[i_sat:]                      # slice saturated state
vpxs = []

print(f'Compiling data outputs...', flush=True)

for i, output in enumerate(sat_outputs):
    data = athena_read.particles(str(output))
    vpxs.append(data['vpx'])
    print('\t{:.0%}'.format(i/len(sat_outputs)), flush=True)

print('\t100%\nComputing histogram and central moments...', flush=True)
vpxs = np.asarray(vpxs)/etav_K                     # normalize for plot
hist, bin_edges = np.histogram(vpxs, bins=bins, density=True)
mean = np.mean(vpxs)
var = np.var(vpxs)
skewness = skew(vpxs, axis=None)
kurt = kurtosis(vpxs, axis=None)
print('... Done.\nSaving...', flush=True)
np.savez_compressed('output/vpxd2', hist=hist, bin_edges=bin_edges,
                    mean=mean, var=var, skew=skewness, kurt=kurt)
