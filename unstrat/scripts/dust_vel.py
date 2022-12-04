#!/usr/bin/env python3
#==============================================================================
# dust_vel.py
#
# Computes and saves the total average, standard deviation, and time-averaged
# probability density function of dust velocities at saturation, including the
# time variability and average and standard deviation of the particle density
# of each bin.
# NOTE: Velocities are normalized to the product of the radial pressure
# gradient and sound speed, Pi*c_s.
#
# Author: Stanley A. Baronett
# Created: 2022-12-01
# Updated: 2022-12-04
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path
from scipy import stats

# Collect Athena++ inputs, outputs, and sim constants
t_sat = float(sys.argv[1])
bins = int(sys.argv[2])
vlim = float(sys.argv[3])
bin_edges = np.linspace(-vlim, vlim, num=bins)
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
c_s = athinput['hydro']['iso_sound_speed']
Pi = athinput['problem']['duy0']
i_sat = int(t_sat/dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                        '.out1.*.athdf')))
outputs = outputs[i_sat:]
vx_stack, vz_stack, rhop_stack = [], [], []
vx_hists, vz_hists = [], []

print(f'Stacking outputs...', flush=True)
for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    vx_stack.append(athdf['vp1']/Pi/c_s)
    vz_stack.append(athdf['vp2']/Pi/c_s)
    rhop_stack.append(athdf['rhop'])
    print(f'  {(i + 1)/len(outputs):3.0%}', flush=True)

print(f'  Done.\nComputing velocity histograms...', flush=True)
vx_stack, vz_stack = np.asarray(vx_stack), np.asarray(vz_stack)
for i in range(vx_stack.shape[0]):
    hist, bin_edges = np.histogram(vx_stack[i], bins=bin_edges, density=True,
                                   weights=rhop_stack[i])
    vx_hists.append(hist)
    hist, bin_edges = np.histogram(vz_stack[i], bins=bin_edges, density=True,
                                   weights=rhop_stack[i])
    vz_hists.append(hist)
    print(f'  {(i + 1)/vx_stack.shape[0]:3.0%}', flush=True)

print(f'  Done.\nComputing density histograms...', flush=True)
vx_flat, vz_flat = vx_stack.ravel(), vz_stack.ravel()
rhop_flat = np.asarray(rhop_stack).ravel()
bin_avg_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='mean', bins=bin_edges)
bin_avg_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='mean', bins=bin_edges)
bin_std_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='std', bins=bin_edges)
bin_std_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='std', bins=bin_edges)

print('  Done.\nComputing statistics...', flush=True)
avg_vxs = np.average(vx_stack, weights=rhop_stack)
avg_vzs = np.average(vz_stack, weights=rhop_stack)
std_vxs = np.sqrt(np.average((vx_stack - avg_vxs)**2, weights=rhop_stack))
std_vzs = np.sqrt(np.average((vz_stack - avg_vzs)**2, weights=rhop_stack))
bin_avg_vxs = np.average(vx_hists, axis=0)
bin_avg_vzs = np.average(vz_hists, axis=0)
bin_std_vxs = np.std(vx_hists, axis=0)
bin_std_vzs = np.std(vz_hists, axis=0)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/dust_vel', bin_edges=bin_edges,
                    avg_vxs=avg_vxs, avg_vzs=avg_vzs,
                    std_vxs=std_vxs, std_vzs=std_vzs,
                    bin_avg_vxs=bin_avg_vxs, bin_avg_vzs=bin_avg_vzs,
                    bin_std_vxs=bin_std_vxs, bin_std_vzs=bin_std_vzs,
                    bin_avg_rhopxs=bin_avg_rhopxs,
                    bin_avg_rhopzs=bin_avg_rhopzs,
                    bin_std_rhopxs=bin_std_rhopxs,
                    bin_std_rhopzs=bin_std_rhopzs)
