#!/usr/bin/env python3
#==============================================================================
# AB_Rs_rad-prof-log-bins-Pi.py
#
# Plot with logarithmically spaced bins azimuthally-averaged radial profiles,
# scaled by Pi, of the normalized autocorrelations of final snapshots of the
# dust and gas density fields across a range of radial pressure gradients for
# case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-28
# Updated: 2022-10-28
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack, stats

def norms(xv, zv, pole):
    rv = np.zeros((len(zv), len(xv)))
    for i, z in enumerate(zv):
        for j, x in enumerate(xv):
            rv[i][j] = np.sqrt((x - pole[0])**2 + (z - pole[1])**2)
    return rv

fig, axs = plt.subplots(2, sharex=True, figsize=(3.15, 4))
workdir = '../../..'
case = 'AB'
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = 2048

for i, Pi in enumerate(Pis):
    # Collect parameters
    print(f'{case}/{Pi[0]}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi[0]}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    c_s = athinput['hydro']['iso_sound_speed']
    etar = float(Pi[0])*c_s
    data = athena_read.athdf(f'{path}/athdf/SI.out1.00100.athdf')
    xv, zv = data['x1v'], data['x2v']
    x0, z0 = int(len(xv)/2), int(len(zv)/2)
    pole = (xv[x0], zv[z0])
    rv = norms(xv, zv, pole).ravel()
    indices = np.where(rv > xv[-1])[0]
    indices = np.append(indices, np.where(rv == 0)[0])
    rv = np.delete(rv, indices)
    r0 = xv[x0]
    base = np.sqrt(2)
    leftmost_edge = r0/np.sqrt(base)
    num = int(np.sqrt(res))
    bin_edges = leftmost_edge*np.logspace(0, num, num=(num + 1), base=base)

    # Process dust
    ft = fftpack.fft2(data['rhop'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    log = np.log10(shift).ravel()
    log = np.delete(log, indices)
    dust_means, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(),
        log.ravel(), statistic='mean', bins=bin_edges)
    dust_stds, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(),
        log.ravel(), statistic='std', bins=bin_edges)
    dust_highs = dust_means + dust_stds
    dust_lows = dust_means - dust_stds

    # Process gas
    diff = data['rho'][0] - 1
    ft = fftpack.fft2(diff)
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    shift = np.delete(shift, indices)
    gas_means, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(),
        shift.ravel(), statistic='mean', bins=bin_edges)
    gas_stds, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(),
        shift.ravel(), statistic='std', bins=bin_edges)
    gas_highs = gas_means + gas_stds
    gas_lows = gas_means - gas_stds
    bin_counts, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(),
        rv.ravel(), statistic='count', bins=bin_edges)

    # Plot histograms
    axs[0].stairs(dust_means, bin_edges/float(Pi[0]), baseline=float('-inf'),
                  color=Pi[1], lw=1.5, label=Pi[0])
    axs[0].stairs(dust_highs, bin_edges/float(Pi[0]), baseline=dust_lows,
                  fill=True, color=Pi[1], alpha=0.2)
    axs[1].stairs(gas_means, bin_edges/float(Pi[0]), baseline=float('-inf'),
                  color=Pi[1], lw=1.5)
    axs[1].stairs(gas_highs, bin_edges/float(Pi[0]), baseline=gas_lows,
                  fill=True, color=Pi[1], alpha=0.2)
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(title=r'$\Pi$')
axs[0].set(ylabel=r'$\log\mathcal{R}_\mathrm{p}$')
axs[1].set(xscale='log', xlabel=r'$r/(\Pi H_\mathrm{g})$',
           ylabel=r'$\mathcal{R}_\mathrm{g}^\prime$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_Rs_rad-prof-log-bins-Pi.pdf', bbox_inches='tight',
            pad_inches=0.01)
