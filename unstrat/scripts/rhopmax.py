#!/usr/bin/env python3
#==============================================================================
# rhopmax.py
#
# Finds and saves the time-varying maximum particle density of the SI run.
#
# Author: Stanley A. Baronett
# Created: 2023-11-05
# Updated: 2023-11-05
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
athinput = athena_read.athinput('athinput.si')
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out2.*.athdf')))
rhopmax, t = [], []

print(f'Compiling data...', flush=True)

for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    rhopmax.append(athdf['rhop'].max)
    t.append(athdf['Time'])
    print('  {:.2%}'.format(i/len(outputs)), flush=True)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/time_series_rhopmax', rhopmax=rhopmax, t=t)
