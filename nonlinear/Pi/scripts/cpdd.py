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
# Last Modified: 2022-02-10
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Passed arguments
# path = sys.argv[1]                           # Relative path to run
t_sat = float(sys.argv[1])                   # / T

# Collect .athdf inputs, outputs, sim consts.
athinput = athena_read.athinput('../athinput.si')
epsilon = athinput['problem']['epsilon']     # avg. BG dust/gas ρ-ratio
dt = athinput['output2']['dt']               # time between rhop outputs
i_sat = int(t_sat / dt)                      # output index of sat. state
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"]+
                                               '.out2.*.athdf')))
rhops = []                                   # sorted dust density snapshots

# Load & process saturated-state data into memory
sat_outputs = outputs[i_sat:]
for output in sat_outputs:
    data = athena_read.athdf(output)
    temp = data['rhop'].flatten() / epsilon  # flatten & convert
    rhops.append(np.sort(temp))              # sort
    prog = 100*len(rhops)/len(sat_outputs)   # percentage done
    print('{:3.1f}%% done.'.format(prog), flush=True)

# Find min., max., avg. of each ordered rhop over saturated state
mins = np.amin(rhops, axis=0)
maxs = np.amax(rhops, axis=0)
avgs = np.average(rhops, axis=0)

# Find std. dev. computed in log space
lnrhops = np.log(rhops)
stds = np.exp(np.std(lnrhops, axis=0))
cdf = np.linspace(1, 0 , mins.size, endpoint=False)

np.savez_compressed('../output/cpdd', cdf=cdf, mins=mins, maxs=maxs,
                    avgs=avgs, stds=stds)
