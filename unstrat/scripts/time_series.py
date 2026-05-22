#!/usr/bin/env python3
#===============================================================================
# time_series.py
#
# Compiles a time series of SI run properties from snapshots; currently
# includes:
#   - max(rhop)
#
# Author: Stanley A. Baronett
# Created: 2026-05-22
# Updated: 2026-05-22
#===============================================================================
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
maxrhop = []
t = []

print(f'Compiling data...', flush=True)

for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    maxrhop.append(np.max(athdf['rhop']))
    t.append(athdf['Time'])
    print('  {:.2%}'.format(i/len(outputs)), flush=True)

print(f'  Done.\nSaving results...', flush=True)
if not os.path.exists('npz'):
    os.makedirs('npz')
np.savez_compressed('npz/time_series', maxrhop=maxrhop, t=t)
