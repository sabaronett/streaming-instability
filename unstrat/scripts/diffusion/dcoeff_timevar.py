#!/usr/bin/env python3
#==============================================================================
# dcoeff_timevar.py
#
# Compute the time average and standard deviation over the saturation state of
# the diffusion coefficient normalized to the product of the sound speed and
# gas scale height, c_s*H_g.
#
# Author: Stanley A. Baronett
# Created: 2023-01-16
# Updated: 2023-01-17
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np

# Collect Athena++ I/O and sim constants
print(f'Loading data...', end='', flush=True)
athinput = athena_read.athinput('athinput.si')
dt = athinput['outputp']['dt']
c_s = athinput['hydro']['iso_sound_speed']
disp = np.load('dat/disp.npz')
dirs = ['dxp', 'dyp']
print(' Done.')

print(f'Copmuting statistics...', end='', flush=True)
for dir in dirs:
    variances = np.var(disp[dir], axis=1)
    dvariances = variances - np.roll(variances, 1)
    dvariances = np.delete(dvariances, 0)
    Ds = dvariances/dt/c_s
    if dir == 'dxp':
        avgDx = np.average(Ds)
        stdDx = np.std(Ds)
    else:
        avgDz = np.average(Ds)
        stdDz = np.std(Ds)
print(' Done.')

print(f'Saving results...', end='', flush=True)
kwds = dict(avgDx = avgDx)
kwds['stdDx'], kwds['avgDz'], kwds['stdDz'] = stdDx, avgDz, stdDz
np.savez_compressed('npz/dcoeff_timevar', **kwds)
print(' Done.')
