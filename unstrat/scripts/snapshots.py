#!/usr/bin/env python3
#==============================================================================
# snapshots.py
#
# Plot suite of final snapshots of the dust and gas density field.
#
# Author: Stanley A. Baronett
# Created: 2022-07-25
# Updated: 2022-07-25
#==============================================================================
import sys
sys.path.insert(0, '/home/stanley/bitbucket/ccyang/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import time

start = time.time()
fig, axs = plt.subplots(2, 4, figsize=(8, 4), dpi=300)
workdir = '../unstrat'
case = sys.argv[1]
Pis = ['0.01', '0.02', '0.05', '0.10']
res = '2048'
vmin, vmax = 1e-1, 1e1                           # CPDD > 85%; CPDD < 5%
if sys.argv[2]: vmin = float(sys.argv[2])
if sys.argv[3]: vmax = float(sys.argv[3])

for i, Pi in enumerate(Pis):
    # Collect data, parameters
    path = f'{workdir}/{case}/{Pi}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    c_s = athinput['hydro']['iso_sound_speed']
    Omega = athinput['problem']['omega']
    H = c_s/Omega                                # gas scale height
    data = athena_read.athdf(f'{path}/athdf/SI.out1.00100.athdf')
    xf, zf = data['x1f']/H, data['x2f']/H
    t = data['Time']

    # Collect and plot densities
    dust = np.clip(data['rhop'][0], vmin, vmax)
    gas  = np.clip(data['rho'][0], vmin, vmax)
    rhops = axs[0][i].pcolormesh(xf, zf, dust, norm=colors.LogNorm())
    rhogs = axs[1][i].pcolormesh(xf, zf, gas, norm=colors.LogNorm())

    axs[0][i].set(title=f'$\Pi={Pi:s}$')

for ax in axs.flat:
    ax.label_outer()
    ax.minorticks_on()
    ax.set(aspect='equal', xlabel=r'$x$ / $H$', ylabel=r'$z$ / $H$')
    ax.tick_params(axis='both', which='both', top=True, right=True)

cb_rhop = fig.colorbar(rhops, ax=axs.flat)
cb_rhop.set_label(r'$\rho$ / $\rho_\mathrm{g,0}$')
fig.suptitle(f'{case}, $t=${t:.0f} $T$')
plt.savefig('figs/ab_snaps.png', dpi=1000, bbox_inches='tight', pad_inches=0.01)

# Print walltime used
wt = time.time() - start
wth = wt//3600
wtm = (wt - wth*3600)//60
wts = wt%60
print(f'Walltime Used : {wt:.0f} s ({wth:02.0f}:{wtm:02.0f}:{wts:02.0f})')
