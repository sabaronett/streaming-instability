#!/usr/bin/env python3
#==============================================================================
# plot_vpxd_pdat.py
#
# Computes, plots, and saves the probability density function histogram of the
# radial velocity distribution of the dust particle data (.dat) during the
# saturated, turbulent state of the SI.
#
# Author: Stanley A. Baronett
# Created: 2022-06-10
# Updated: 2022-06-10
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
fig, axs = plt.subplots(2, 1, figsize=(6,8), dpi=150)
runs = ['AB', 'BA']
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = '2048'
n_bins = 100

for i,ax in enumerate(axs.flat):
    for Pi in Pis:
        path = '../%s/%s/%s/'%(runs[i], Pi[0], res)
        outputs = sorted(list(Path(path+'dat').glob('SI.pout.*.dat')))
        time, pdata = athena_read.particles(outputs[0])  # equilibrium
        vpx0 = np.average(pdata['vpx'])
        time, pdata = athena_read.particles(outputs[-1]) # last snapshot

        ax.axvline(vpx0, color=Pi[1], ls='--',
                   label=r'$v_{\mathrm{p},x,0},\,0\,T$')
        ax.hist(pdata['vpx'], bins=n_bins, density=True, histtype='step',
                color=Pi[1], label=r'{:s}, {:.0f} $T$'.format(Pi[0], time))

    ax.grid()
    ax.legend(title=r'$\Pi,\,t$')
    ax.minorticks_on()
    ax.set(ylabel=r'd$f$ / d$v_{\mathrm{p},x}$')
    ax.text(0.03, 0.97, runs[i], ha='left', va='top', size='x-large',
            transform=ax.transAxes)
    ax.tick_params(which='both', top=True, right=True)

ax.set(xlabel=r'$v_{\mathrm{p},x}$ $[H/T]$')
plt.savefig('output/vpxd_pdat.pdf', bbox_inches='tight', pad_inches=0.01)
