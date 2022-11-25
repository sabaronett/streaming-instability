#!/usr/bin/env python3
#==============================================================================
# gas_vel.py
#
# Computes and saves the time-averaged probability density function of gas
# velocities at saturation, including the time variability of each bin.
# NOTE: Velocities are normalized to the product of the radial pressure
# gradient and sound speed, Pi*c_s.
#
# Author: Stanley A. Baronett
# Created: 2022-11-25
# Updated: 2022-11-25
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
t_sat = float(sys.argv[1])
bins = int(sys.argv[2])
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
c_s = athinput['hydro']['iso_sound_speed']
Pi = athinput['problem']['duy0']
i_sat = int(t_sat/dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                        '.out1.*.athdf')))
outputs = outputs[i_sat:]
ux_stack, uz_stack, rho_stack = [], [], []
ux_hists, uz_hists = [], []

print(f'Stacking outputs...', flush=True)
for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    uxs, uzs, rhos = athdf['vel1']/Pi/c_s, athdf['vel2']/Pi/c_s, athdf['rho']
    ux_stack.append(uxs)
    uz_stack.append(uzs)
    rho_stack.append(rhos)
    print(f'  {i/len(outputs):3.0%}', flush=True)

ux_stack, uz_stack = np.asarray(ux_stack), np.asarray(uz_stack)
print(f'Finding bin edges...', flush=True)
edge = np.abs(ux_stack.min())
if np.abs(ux_stack.max()) > edge: edge = np.abs(ux_stack.max())
if np.abs(uz_stack.min()) > edge: edge = np.abs(uz_stack.min())
if np.abs(uz_stack.max()) > edge: edge = np.abs(uz_stack.max())
bin_edges = np.linspace(-edge, edge, num=bins)

print(f'  100%\nComputing histograms...', flush=True)
for i in range(ux_stack.shape[0]):
    hist, bin_edges = np.histogram(ux_stack[i], bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    ux_hists.append(hist)
    hist, bin_edges = np.histogram(uz_stack[i], bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    uz_hists.append(hist)
    print(f'  {i/ux_stack.shape[0]:3.0%}', flush=True)

print('  100%\nComputing statistical quantities...', flush=True)
avg_uxs, avg_uzs = np.average(ux_hists, axis=0), np.average(uz_hists, axis=0)
std_uxs, std_uzs = np.std(ux_hists, axis=0), np.std(uz_hists, axis=0)

np.savez_compressed('npz/gas_vel', avg_uxs=avg_uxs, avg_uzs=avg_uzs,
                    std_uxs=std_uxs, std_uzs=std_uzs, bin_edges=bin_edges)
