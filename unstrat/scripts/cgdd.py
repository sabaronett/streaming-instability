#!/usr/bin/env python3
#==============================================================================
# cgdd.py
#
# Computes and saves the time-averaged cumulative gas density distribution
# of the SI run during the saturated, turbulent state, including time-varying
# minima, maxima, and standard deviations.
#
# Author: Stanley A. Baronett
# Created: 2022-08-02
# Updated: 2022-08-03
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
t_sat = float(sys.argv[1])
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
i_sat = int(t_sat / dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out1.*.athdf')))
sat_outputs = outputs[i_sat:]
rhogs = []

print(f'Compiling data...', flush=True)

for i, output in enumerate(sat_outputs):
    data = athena_read.athdf(output)
    sort = np.sort(data['rho'], axis=None)
    rhogs.append(sort)
    print('  {:.2%}'.format(i/len(sat_outputs)), flush=True)

print('  100%\nComputing statistical quantities...', flush=True)
rhogs = np.asarray(rhogs)
mins = np.amin(rhogs, axis=0)
maxs = np.amax(rhogs, axis=0)
avgs = np.average(rhogs, axis=0)
stds = np.std(rhogs, axis=0)
cdf = np.linspace(1, 0, avgs.size, endpoint=False)

np.savez_compressed('output/cgdd', cdf=cdf, mins=mins, maxs=maxs, avgs=avgs,
                    stds=stds)
