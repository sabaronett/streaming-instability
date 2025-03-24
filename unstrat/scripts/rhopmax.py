#!/usr/bin/env python3
#==============================================================================
# rhopmax.py
#
# Finds and saves the time-varying maximum particle density of the SI run.
#
# Author: Stanley A. Baronett
# Created: 2023-11-05
# Updated: 2025-03-23
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
import os
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
athinput = athena_read.athinput('athinput.si')
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                         '.out2.*.athdf')))
rhopmax, t = [], []

print(f'Compiling data...', flush=True)

for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    rhopmax.append(np.max(athdf['rhop']))
    t.append(athdf['Time'])
    print('  {:.2%}'.format(i/len(outputs)), flush=True)

print(f'  Done.\nSaving results...', flush=True)
if not os.path.exists('npz'):
    os.makedirs('npz')  # create directory as needed
np.savez_compressed('npz/time_series_rhopmax', rhopmax=rhopmax, t=t)
