#!/usr/bin/env python3
#==============================================================================
# vxs.py
#
# Computes and saves the time-averaged drift velocity distribution of the SI
# run during the saturated, turbulent state.
#
# Author: Stanley A. Baronett
# Created: 2022-03-02
# Last Modified: 2022-03-02
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

t_sat = float(sys.argv[1])                # / T
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']            # time between vp1 outputs
i_sat = int(t_sat / dt)                   # output index of sat. state
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out1.*.athdf')))
vp1s = []                                 # sorted dust density snapshots
sat_outputs = outputs[i_sat:]             # slice saturated state

for output in sat_outputs:
    data = athena_read.athdf(output)
    temp = data['vp1'].flatten()
    vp1s.append(np.sort(temp))            # sort
    prog = 100*len(vp1s)/len(sat_outputs) # percentage progress

    print('{:3.1f}%% done.'.format(prog), flush=True)

avgs = np.average(vp1s, axis=0)           # avg. ordered vx

np.savez_compressed('output/vxs', avgs=avgs)
