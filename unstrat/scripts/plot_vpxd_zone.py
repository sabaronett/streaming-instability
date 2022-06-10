#!/usr/bin/env python3
#==============================================================================
# plot_vpxd_zone.py
#
# Computes, plots, and saves the probability density function histogram of the
# zonal radial velocity distribution, weighted by particle densities, of the
# dust during the saturated, turbulent state of the SI.
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
arrays, avgvpxs = [[], []], []
t_sats = [2.0, 20.0]

for i,ax in enumerate(axs.flat):
    for Pi in Pis:
        path = '../%s/%s/%s/'%(runs[i], Pi[0], res)
        athinput = athena_read.athinput(path+'athinput.si')
        i_sat = int(t_sats[i]/athinput['output1']['dt'])
        if runs[i]=='AB' and Pi[0]=='0.01':
            i_sat = int(4.0/athinput['output1']['dt'])
        hst = athena_read.hst(path+'output/SI.hst')
        Np = athinput['problem']['npx1']*athinput['problem']['npx2']\
             *athinput['problem']['npx3']
        vpx0 = hst['vp1'][0]/Np
        outputs = sorted(list(Path(path+'athdf').glob(
                  athinput["job"]["problem_id"]+'.out1.*.athdf')))
        output = outputs[-1] # last snapshot
        data = athena_read.athdf(output)
        vpxs = data['vp1'].flatten()
        rhops = data['rhop'].flatten()
        t = data['Time']

        arrays[0].append(runs[i])
        arrays[1].append(Pi[0])
        avgvpxs.append(np.average(vpxs, weights=rhops))
        ax.hist(vpxs, bins=n_bins, density=True, weights=rhops, color=Pi[1],
                histtype='step', label=r'{:s}, {:.0f} $T$'.format(Pi[0], t))
        ax.axvline(vpx0, color=Pi[1], ls='--',
                   label=r'$v_{\mathrm{p},x,0},\,0\,T$')

    ax.grid()
    ax.legend(title=r'$\Pi,\,t$')
    ax.minorticks_on()
    ax.set(ylabel=r'd$f_\mathrm{zone}$ / d$v_{\mathrm{p},x}$')
    ax.text(0.03, 0.97, runs[i], ha='left', va='top', size='x-large',
            transform=ax.transAxes)
    ax.tick_params(which='both', top=True, right=True)

ax.set(xlabel=r'$v_{\mathrm{p},x}$ $[H/T]$')
plt.savefig('output/vpxd_zone.pdf', bbox_inches='tight', pad_inches=0.01)
