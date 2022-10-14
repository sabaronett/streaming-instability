#!/usr/bin/env python3
#==============================================================================
# AB_Rs_rad-prof.py
#
# Plot azimuthally-averaged radial profiles of the normalized autocorrelations
# of final snapshots of the dust and gas density fields across a range of
# radial pressure gradients for case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-14
# Updated: 2022-10-14
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack, stats

def norms(xv, zv):
    rv = np.zeros((len(zv), len(xv)))
    for i, z in enumerate(zv):
        for j, x in enumerate(xv):
            rv[i][j] = np.sqrt(x**2 + z**2)
    return rv

fig, axs = plt.subplots(2, sharex=True, figsize=(3.15, 6), dpi=300)
workdir = '../../..'
case = 'AB'
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = 2048
bins = 100

for i, Pi in enumerate(Pis):
    # Collect parameters and plot densities
    print(f'{case}/{Pi[0]}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi[0]}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    c_s = athinput['hydro']['iso_sound_speed']
    etar = float(Pi[0])*c_s
    data = athena_read.athdf(f'{path}/athdf/SI.out1.00100.athdf')
    xv, zv = data['x1v'], data['x2v']
    rv = norms(xv, zv)
    ft = fftpack.fft2(data['rhop'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    log = np.log10(shift)
    dust_means, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(), log.ravel(),
        statistic='mean', bins=bins)
    dust_stds, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(), log.ravel(),
        statistic='std', bins=100)
    dust_highs = dust_means + dust_stds
    dust_lows = dust_means - dust_stds
    ft = fftpack.fft2(data['rho'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    offset = (shift - 1)*1e11
    gas_means, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(), offset.ravel(),
        statistic='mean', bins=bins)
    gas_stds, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(), offset.ravel(),
        statistic='std', bins=bins)
    gas_highs = gas_means + gas_stds
    gas_lows = gas_means - gas_stds
    bin_counts, bin_edges, binnumnber = stats.binned_statistic(rv.ravel(), rv.ravel(),
        statistic='count', bins=bins)

    axs[0].stairs(dust_means, bin_edges/etar, baseline=float('-inf'), color=Pi[1],
                  lw=1.5, label=Pi[0])
    axs[0].stairs(dust_highs, bin_edges/etar, baseline=dust_lows, fill=True,
                                 color=Pi[1], alpha=0.2)
    axs[1].stairs(gas_means, bin_edges/etar, baseline=float('-inf'), color=Pi[1], lw=1.5)
    axs[1].stairs(gas_highs, bin_edges/etar, baseline=gas_lows, fill=True,
                                 color=Pi[1], alpha=0.2)
    axs[2].stairs(bin_counts, bin_edges/etar, color=Pi[1], lw=1.5)
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(title=r'$\Pi$')
axs[0].set(ylabel=r'$\log\left(\int_0^{2\pi}\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}\mathrm{d}\phi\right)$')
axs[1].set(yscale='symlog',
           ylabel=r'$\left(\int_0^{2\pi}\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}\mathrm{d}\phi-1\right)\times10^{11}$')
axs[2].set(xscale='log', yscale='log', xlabel=r'$r^\prime/(\eta r)$', 
           ylabel=r'Bin Count')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_Rs_rad-prof.pdf', bbox_inches='tight', pad_inches=0.01)
