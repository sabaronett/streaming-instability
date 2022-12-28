#!/usr/bin/env python3
#==============================================================================
# velocities.py
#
# Computes and saves the total average, standard deviation, and time-averaged
# probability density function of gas and dust velocities at saturation,
# including the time variability in logarithmic space and binned statistics of
# the gas or particle density distribution within each bin.
# NOTE: Velocities are normalized to the product of the radial pressure
# gradient and sound speed, Pi*c_s.
#
# Author: Stanley A. Baronett
# Created: 2022-12-01
# Updated: 2022-12-28
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
lim = float(sys.argv[3])
bin_edges = np.linspace(-lim, lim, num=(bins + 1))
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
c_s = athinput['hydro']['iso_sound_speed']
Pi = athinput['problem']['duy0']
i_sat = int(t_sat/dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                        '.out1.*.athdf')))
outputs = outputs[i_sat:]
# Gas
ux_stack, uz_stack, rho_stack = [], [], []
ux_hists, uz_hists = [], []
# Dust
vx_stack, vz_stack, rhop_stack = [], [], []
vx_hists, vz_hists = [], []

print(f'Stacking outputs...', flush=True)
for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    # Gas
    ux_stack.append(athdf['vel1']/Pi/c_s)
    uz_stack.append(athdf['vel2']/Pi/c_s)
    rho_stack.append(athdf['rho'])
    # Dust
    vx_stack.append(athdf['vp1']/Pi/c_s)
    vz_stack.append(athdf['vp2']/Pi/c_s)
    rhop_stack.append(athdf['rhop'])
    print(f'  {(i + 1)/len(outputs):3.0%}', flush=True)

print(f'  Done.\nComputing velocity histograms...', flush=True)
for i, ux in enumerate(ux_stack):
    # Gas
    hist, bin_edges = np.histogram(ux, bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    ux_hists.append(hist)
    hist, bin_edges = np.histogram(uz_stack[i], bins=bin_edges, density=True,
                                   weights=rho_stack[i])
    uz_hists.append(hist)
    # Dust
    hist, bin_edges = np.histogram(vx_stack[i], bins=bin_edges, density=True,
                                   weights=rhop_stack[i])
    vx_hists.append(hist)
    hist, bin_edges = np.histogram(vz_stack[i], bins=bin_edges, density=True,
                                   weights=rhop_stack[i])
    vz_hists.append(hist)
    print(f'  {(i + 1)/len(ux_stack):3.0%}', flush=True)

print('  Done.\nComputing velocity statistics...', flush=True)
# Gas
ux_stack, uz_stack = np.stack(ux_stack), np.stack(uz_stack)
avg_uxs = np.average(ux_stack, weights=rho_stack)
avg_uzs = np.average(uz_stack, weights=rho_stack)
std_uxs = np.sqrt(np.average((ux_stack - avg_uxs)**2, weights=rho_stack))
std_uzs = np.sqrt(np.average((uz_stack - avg_uzs)**2, weights=rho_stack))
bin_avg_uxs = np.average(ux_hists, axis=0)
bin_avg_uzs = np.average(uz_hists, axis=0)
bin_std_uxs = np.std(ux_hists, axis=0)
bin_std_uzs = np.std(uz_hists, axis=0)
bin_log_std_uxs = np.exp(np.std(np.log(ux_hists), axis=0))
bin_log_std_uzs = np.exp(np.std(np.log(uz_hists), axis=0))
# Dust
vx_stack, vz_stack = np.stack(vx_stack), np.stack(vz_stack)
avg_vxs = np.average(vx_stack, weights=rhop_stack)
avg_vzs = np.average(vz_stack, weights=rhop_stack)
std_vxs = np.sqrt(np.average((vx_stack - avg_vxs)**2, weights=rhop_stack))
std_vzs = np.sqrt(np.average((vz_stack - avg_vzs)**2, weights=rhop_stack))
bin_avg_vxs = np.average(vx_hists, axis=0)
bin_avg_vzs = np.average(vz_hists, axis=0)
bin_std_vxs = np.std(vx_hists, axis=0)
bin_std_vzs = np.std(vz_hists, axis=0)
bin_log_std_vxs = np.exp(np.std(np.log(vx_hists), axis=0))
bin_log_std_vzs = np.exp(np.std(np.log(vz_hists), axis=0))

print(f'  Done.\nComputing density statistics within each velocity bin...',
      flush=True)
# Gas
ux_flat, uz_flat = ux_stack.ravel(), uz_stack.ravel()
rho_flat = np.asarray(rho_stack).ravel()
indices = np.where(rho_flat == 0)[0]
ux_flat, uz_flat = np.delete(ux_flat, indices), np.delete(uz_flat, indices)
rho_flat = np.delete(rho_flat, indices)
bin_avg_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='mean', bins=bin_edges)
bin_avg_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='mean', bins=bin_edges)
bin_std_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='std', bins=bin_edges)
bin_std_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='std', bins=bin_edges)
bin_med_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='median', bins=bin_edges)
bin_med_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='median', bins=bin_edges)
bin_cnt_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='count', bins=bin_edges)
bin_cnt_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='count', bins=bin_edges)
bin_sum_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='sum', bins=bin_edges)
bin_sum_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='sum', bins=bin_edges)
bin_min_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='min', bins=bin_edges)
bin_min_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='min', bins=bin_edges)
bin_max_rhoxs, bin_edges, binnumnber = stats.binned_statistic(ux_flat,
    rho_flat, statistic='max', bins=bin_edges)
bin_max_rhozs, bin_edges, binnumnber = stats.binned_statistic(uz_flat,
    rho_flat, statistic='max', bins=bin_edges)
# Dust
vx_flat, vz_flat = vx_stack.ravel(), vz_stack.ravel()
rhop_flat = np.asarray(rhop_stack).ravel()
indices = np.where(rhop_flat == 0)[0]
vx_flat, vz_flat = np.delete(vx_flat, indices), np.delete(vz_flat, indices)
rhop_flat = np.delete(rhop_flat, indices)
log_rhop_flat = np.log(rhop_flat)
bin_avg_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='mean', bins=bin_edges)
bin_avg_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='mean', bins=bin_edges)
bin_std_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='std', bins=bin_edges)
bin_std_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='std', bins=bin_edges)
bin_log_std_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    log_rhop_flat, statistic='std', bins=bin_edges)
bin_log_std_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    log_rhop_flat, statistic='std', bins=bin_edges)
bin_med_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='median', bins=bin_edges)
bin_med_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='median', bins=bin_edges)
bin_cnt_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='count', bins=bin_edges)
bin_cnt_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='count', bins=bin_edges)
bin_sum_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='sum', bins=bin_edges)
bin_sum_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='sum', bins=bin_edges)
bin_min_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='min', bins=bin_edges)
bin_min_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='min', bins=bin_edges)
bin_max_rhopxs, bin_edges, binnumnber = stats.binned_statistic(vx_flat,
    rhop_flat, statistic='max', bins=bin_edges)
bin_max_rhopzs, bin_edges, binnumnber = stats.binned_statistic(vz_flat,
    rhop_flat, statistic='max', bins=bin_edges)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/velocities', bin_edges=bin_edges,
                    # Gas
                    avg_uxs=avg_uxs, avg_uzs=avg_uzs,
                    std_uxs=std_uxs, std_uzs=std_uzs,
                    bin_avg_uxs=bin_avg_uxs, bin_avg_uzs=bin_avg_uzs,
                    bin_std_uxs=bin_std_uxs, bin_std_uzs=bin_std_uzs,
                    bin_log_std_uxs=bin_log_std_uxs,
                    bin_log_std_uzs=bin_log_std_uzs,
                    bin_avg_rhoxs=bin_avg_rhoxs,
                    bin_avg_rhozs=bin_avg_rhozs,
                    bin_std_rhoxs=bin_std_rhoxs,
                    bin_std_rhozs=bin_std_rhozs,
                    bin_med_rhoxs=bin_med_rhoxs,
                    bin_med_rhozs=bin_med_rhozs,
                    bin_cnt_rhoxs=bin_cnt_rhoxs,
                    bin_cnt_rhozs=bin_cnt_rhozs,
                    bin_sum_rhoxs=bin_sum_rhoxs,
                    bin_sum_rhozs=bin_sum_rhozs,
                    bin_min_rhoxs=bin_min_rhoxs,
                    bin_min_rhozs=bin_min_rhozs,
                    bin_max_rhoxs=bin_max_rhoxs,
                    bin_max_rhozs=bin_max_rhozs,
                    # Dust
                    avg_vxs=avg_vxs, avg_vzs=avg_vzs,
                    std_vxs=std_vxs, std_vzs=std_vzs,
                    bin_avg_vxs=bin_avg_vxs, bin_avg_vzs=bin_avg_vzs,
                    bin_std_vxs=bin_std_vxs, bin_std_vzs=bin_std_vzs,
                    bin_log_std_vxs=bin_log_std_vxs,
                    bin_log_std_vzs=bin_log_std_vzs,
                    bin_avg_rhopxs=bin_avg_rhopxs,
                    bin_avg_rhopzs=bin_avg_rhopzs,
                    bin_std_rhopxs=bin_std_rhopxs,
                    bin_std_rhopzs=bin_std_rhopzs,
                    bin_log_std_rhopxs=bin_log_std_rhopxs,
                    bin_log_std_rhopzs=bin_log_std_rhopzs,
                    bin_med_rhopxs=bin_med_rhopxs,
                    bin_med_rhopzs=bin_med_rhopzs,
                    bin_cnt_rhopxs=bin_cnt_rhopxs,
                    bin_cnt_rhopzs=bin_cnt_rhopzs,
                    bin_sum_rhopxs=bin_sum_rhopxs,
                    bin_sum_rhopzs=bin_sum_rhopzs,
                    bin_min_rhopxs=bin_min_rhopxs,
                    bin_min_rhopzs=bin_min_rhopzs,
                    bin_max_rhopxs=bin_max_rhopxs,
                    bin_max_rhopzs=bin_max_rhopzs)
