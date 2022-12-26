#!/usr/bin/env python3
#==============================================================================
# gas_densities.py
#
# Computes and saves the time-averaged probability density function of gas
# densities at saturation, including the time variability, in linear and
# logarithmic spaces, and various binned statistics.
#
# Author: Stanley A. Baronett
# Created: 2022-12-26
# Updated: 2022-12-26
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
bin_edges = np.linspace((1 - lim), (1 + lim), num=(bins + 1))
athinput = athena_read.athinput('athinput.si')
dt = athinput['output1']['dt']
i_sat = int(t_sat/dt)
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
                                        '.out1.*.athdf')))
outputs = outputs[i_sat:]
rho_stack, rho_hists = [], []

print(f'Stacking outputs...', flush=True)
for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    rho_stack.append(athdf['rho'])
    print(f'  {(i + 1)/len(outputs):3.0%}', flush=True)

print(f'  Done.\nComputing densitiy histograms...', flush=True)
for i, rho in enumerate(rho_stack):
    hist, bin_edges = np.histogram(rho, bins=bin_edges, density=True)
    rho_hists.append(hist)
    print(f'  {(i + 1)/len(rho_stack):3.0%}', flush=True)

print('  Done.\nComputing density statistics...', flush=True)
bin_avg_rhos = np.average(rho_hists, axis=0)
bin_std_rhos = np.std(rho_hists, axis=0)
bin_log_std_rhos = np.exp(np.std(np.log(rho_hists), axis=0))
rho_flat = np.stack(rho_stack).ravel()
bin_mean_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='mean', bins=bin_edges)
bin_std2_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='std', bins=bin_edges)
bin_med_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='median', bins=bin_edges)
bin_cnt_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='count', bins=bin_edges)
bin_sum_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='sum', bins=bin_edges)
bin_min_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='min', bins=bin_edges)
bin_max_rhos, bin_edges, binnumnber = stats.binned_statistic(rho_flat,
    rho_flat, statistic='max', bins=bin_edges)

print(f'  Done.\nSaving results...', flush=True)
np.savez_compressed('npz/gas_densities', bin_edges=bin_edges,
                    bin_avg_rhos=bin_avg_rhos,
                    bin_std_rhos=bin_std_rhos,
                    bin_log_std_rhos=bin_log_std_rhos,
                    bin_mean_rhos=bin_mean_rhos,
                    bin_std2_rhos=bin_std2_rhos,
                    bin_med_rhos=bin_med_rhos,
                    bin_cnt_rhos=bin_cnt_rhos,
                    bin_sum_rhos=bin_sum_rhos,
                    bin_min_rhos=bin_min_rhos,
                    bin_max_rhos=bin_max_rhos)
