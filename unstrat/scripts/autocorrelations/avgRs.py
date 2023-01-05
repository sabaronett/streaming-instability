#!/usr/bin/env python3
#==============================================================================
# avgRs.py
#
# Compute and save the time average and variability of normalized auto-
# correlation of snapshot of the dust and gas density fields.
#
# Author: Stanley A. Baronett
# Created: 2022-12-26
# Updated: 2023-01-05
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path
from scipy import fftpack

# Collect Athena++ inputs, outputs, and sim constants
t_sat = float(sys.argv[1])
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
i_sat = int(t_sat/dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"] +
                                         '.out1.*.athdf')))
outputs = outputs[i_sat:]
epsilon = athinput['problem']['epsilon']
Rps, Rgs = [], []

print(f'Computing sptial correlations and stacking outputs...', flush=True)
for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)

    # Process gas
    drho = athdf['rho'][0] - 1
    ft = fftpack.fft2(drho)
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    Rgs.append(shift)

    # Process dust
    drho = athdf['rhop'][0] - epsilon
    ft = fftpack.fft2(drho)
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    Rps.append(shift)

    print(f'  {(i + 1)/len(outputs):3.0%}', flush=True)

print('  Done.\nComputing statistics...', flush=True)
avgRgs, avgRps = np.average(Rgs, axis=0), np.average(Rps, axis=0)
stdRgs, stdRps = np.std(Rgs, axis=0), np.std(Rps, axis=0)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/avgRs', avgRgs=avgRgs, avgRps=avgRps,
                    stdRgs=stdRgs, stdRps=stdRps)
