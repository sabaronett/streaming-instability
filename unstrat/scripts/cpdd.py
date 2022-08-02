#!/usr/bin/env python3
#==============================================================================
# cpdd.py
#
# Computes and saves the time-averaged cumulative particle density distribution
# of the SI run during the saturated, turbulent state, including time-varying
# minima, maxima, and standard deviations in logarithmic space.
#
# Author: Stanley A. Baronett
# Created: 2022-02-08
# Updated: 2022-08-02
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
t_sat = float(sys.argv[1])
athinput = athena_read.athinput('athinput.si')
dt = athinput['output2']['dt']
i_sat = int(t_sat / dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out2.*.athdf')))
sat_outputs = outputs[i_sat:]
rhops = []

print(f'Compiling data...', flush=True)

for i, output in enumerate(sat_outputs):
    data = athena_read.athdf(output)
    sort = np.sort(data['rhop'], axis=None)
    rhops.append(sort)
    print('  {:.2%}'.format(i/len(sat_outputs)), flush=True)

print('  100%\nComputing statistical quantities...', flush=True)
rhops = np.asarray(rhops)
mins = np.amin(rhops, axis=0)
maxs = np.amax(rhops, axis=0)
avgs = np.average(rhops, axis=0)
log = np.log(rhops)
std_log = np.std(log, axis=0)
stds = np.exp(std_log)
cdf = np.linspace(1, 0, avgs.size, endpoint=False)

np.savez_compressed('output/cpdd', cdf=cdf, mins=mins, maxs=maxs, avgs=avgs,
                    stds=stds)
