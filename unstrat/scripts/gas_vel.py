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
# Updated: 2022-11-30
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
t_sat = float(sys.argv[1])
bins = int(sys.argv[2])
ulim = float(sys.argv[3])
bin_edges = np.linspace(-ulim, ulim, num=bins)
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
    print(f'  {(i + 1)/len(outputs):3.0%}', flush=True)

print(f'  Done.\nComputing histograms...', flush=True)
ux_stack, uz_stack = np.asarray(ux_stack), np.asarray(uz_stack)
for i in range(ux_stack.shape[0]):
    hist, bin_edges = np.histogram(ux_stack[i], bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    ux_hists.append(hist)
    hist, bin_edges = np.histogram(uz_stack[i], bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    uz_hists.append(hist)
    print(f'  {(i + 1)/ux_stack.shape[0]:3.0%}', flush=True)

print('Computing statistics...', flush=True)
avg_uxs, avg_uzs = np.average(ux_stack), np.average(uz_stack)
std_uxs, std_uzs = np.std(ux_stack), np.std(uz_stack)
bin_avg_uxs = np.average(ux_hists, axis=0)
bin_avg_uzs = np.average(uz_hists, axis=0)
bin_std_uxs = np.std(ux_hists, axis=0)
bin_std_uzs = np.std(uz_hists, axis=0)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/gas_vel', bin_edges=bin_edges,
                    avg_uxs=avg_uxs, avg_uzs=avg_uzs,
                    std_uxs=std_uxs, std_uzs=std_uzs,
                    bin_avg_uxs=bin_avg_uxs, bin_avg_uzs=bin_avg_uzs,
                    bin_std_uxs=bin_std_uxs, bin_std_uzs=bin_std_uzs)
