#!/usr/bin/env python3
#==============================================================================
# vpxd_mini.py
#
# Plots the density-weighted, zonal radial velocity distributions of the dust
# during the saturated, turbulent state of the SI across multiple runs.
#
# Author: Stanley A. Baronett
# Created: 2022-03-07
# Last Modified: 2022-03-07
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# Passed arguments
last = int(sys.argv[1])                        # only process last many outputs

fig, axs = plt.subplots(2, 1, figsize=(6,8), dpi=300)
runs = ['AB', 'BA']
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = '2048'
xlims = [(-1.0, 1.0), (-2.0, 2.0)]
n_bins = 50
arrays, avgvpxs = [[], []], []
t_sats = [2.0, 20.0]

for i,ax in enumerate(axs.flat):
    for Pi in Pis:
        path = '%s/%s/%s/'%(runs[i], Pi[0], res)
        print(f'Processing '+path+'...')
        athinput = athena_read.athinput(path+'athinput.si')
        i_sat = int(t_sats[i]/athinput['output1']['dt'])
        if runs[i]=='AB' and Pi[0]=='0.01':
            i_sat = int(4.0/athinput['output1']['dt'])
        etav_K = float(Pi[0])*athinput['hydro']['iso_sound_speed']
        hst = athena_read.hst(path+'output/SI.hst')
        Np = athinput['problem']['npx1']*athinput['problem']['npx2']\
             *athinput['problem']['npx3']
        vpx0 = hst['vp1'][0]/Np/etav_K
        outputs = sorted(list(Path(path+'athdf').glob(
                  athinput["job"]["problem_id"]+'.out1.*.athdf')))
        sat_outputs = outputs[-last:]
        data = athena_read.athdf(sat_outputs[0])
        t = data['Time']
        label = r'{:s}, {:.1f} $T$'.format(Pi[0], t)
        if runs[i]=='BA':
            label = r'{:s}, {:.0f} $T$'.format(Pi[0], t)
        vpxs, rhops = [], []

        for j,output in enumerate(sat_outputs):
            data = athena_read.athdf(output)
            vpxs = np.append(vpxs, data['vp1'].flatten())
            rhops = np.append(rhops, data['rhop'].flatten())
            print('  {:.0%} done'.format(j/len(sat_outputs)))

        print('  Done.')
        arrays[0].append(runs[i])
        arrays[1].append(Pi[0])
        avgvpxs.append(np.average(vpxs, weights=rhops)/etav_K)
        ax.hist(vpxs/etav_K, bins=n_bins, density=True, weights=rhops,
                histtype='step', label=label)

    ax.axvline(vpx0, c='black', ls='--', label=r'$v_{\mathrm{p},x,0},\,0\,T$')
    ax.grid()
    ax.legend(loc='upper right', title=r'$\Pi,\,t\geq$')
    ax.minorticks_on()
    ax.set(xlim=xlims[i], ylabel=r'd$f_\mathrm{zone}$ / d$v_{\mathrm{p},x}$')
    ax.text(0.03, 0.97, runs[i], ha='left', va='top', size='xx-large',
            transform=ax.transAxes)
    ax.tick_params(which='both', top=True, right=True)

ax.set(xlabel=r'$v_{\mathrm{p},x}$ / $(\eta v_\mathrm{K})$')
plt.savefig(f'scripts/figs/vpxd_mini_{last:02}.pdf', bbox_inches='tight', pad_inches=0.01)

tuples = list(zip(*arrays))
names = ['Case', '$\Pi$']
index = pd.MultiIndex.from_tuples(tuples, names=names)
df = pd.DataFrame({r'$\langle v_{\textrm{p},x} \rangle$':avgvpxs}, index=index)
print('=====TABLE=====\n', pd)
