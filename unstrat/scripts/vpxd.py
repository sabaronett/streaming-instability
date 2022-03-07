#!/usr/bin/env python3
#==============================================================================
# vpxd.py
#
# Compiles and saves the flattened, zonal radial velocities and particle
# densities of the dust during the saturated, turbulent state of the SI.
#
# Author: Stanley A. Baronett
# Created: 2022-03-02
# Last Modified: 2022-03-06
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

t_sat = float(sys.argv[1])     # [T]
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt'] # time between vp1 outputs
i_sat = int(t_sat/dt)          # output index of saturated state
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out1.*.athdf')))
sat_outputs = outputs[i_sat:]  # slice saturated state
vp1s, rhops = [], []

for i,output in enumerate(sat_outputs):
    data = athena_read.athdf(output)
    vp1s = np.append(vp1s, data['vp1'].flatten())
    rhops = np.append(rhops, data['rhop'].flatten())

    print(f'{i/len(sat_outputs):.1%} done.')

np.savez_compressed('output/vpxd', vp1s=vp1s, rhops=rhops)
