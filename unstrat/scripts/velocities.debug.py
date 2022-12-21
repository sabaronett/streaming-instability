#!/usr/bin/env python3
#==============================================================================
# velocities.debug.py
#
# NOTE: This version is to debug bins with a non-zero average but zero standard
# deviation and only computes dP/du_z.
#
# Computes and saves the time-averaged probability density function of gas
# velocities at saturation, including the time variability.
# NOTE: Velocities are normalized to the product of the radial pressure
# gradient and sound speed, Pi*c_s.
#
# Author: Stanley A. Baronett
# Created: 2022-12-21
# Updated: 2022-12-21
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path
from scipy import stats

# Collect Athena++ inputs, outputs, and sim constants
last = int(sys.argv[1])                                 # last no. of snapshots
bins = int(sys.argv[2])
lim = float(sys.argv[3])
bin_edges = np.linspace(-lim, lim, num=(bins + 1))
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
c_s = athinput['hydro']['iso_sound_speed']
Pi = athinput['problem']['duy0']
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                        '.out1.*.athdf')))
outputs = outputs[-last:]
uz_stack, rho_stack = [], []
uz_hists = []

print(f'Stacking outputs...', flush=True)
for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    uz_stack.append(athdf['vel2']/Pi/c_s)
    rho_stack.append(athdf['rho'])
    print(f'  {(i + 1)/len(outputs):3.0%}', flush=True)

print(f'  Done.\nComputing velocity histograms...', flush=True)
uz_stack = np.asarray(uz_stack)                              # Change to np.stack
for i in range(uz_stack.shape[0]):
    hist, bin_edges = np.histogram(uz_stack[i], bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    uz_hists.append(hist)
    print(f'  {(i + 1)/uz_stack.shape[0]:3.0%}', flush=True)

print('  Done.\nComputing velocity statistics...', flush=True)
# uz_hists = np.stack(uz_hists)
bin_avg_uzs = np.average(uz_hists, axis=0)
bin_std_uzs = np.std(uz_hists, axis=0)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/velocities.debug', bin_edges=bin_edges,
                    uz_stack=uz_stack,
                    rho_stack=rho_stack,
                    uz_hists=uz_hists,
                    bin_avg_uzs=bin_avg_uzs,
                    bin_std_uzs=bin_std_uzs,)
