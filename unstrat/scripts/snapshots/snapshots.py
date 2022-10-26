#!/usr/bin/env python3
#==============================================================================
# snapshots.py
#
# Plot suite of final snapshots of the dust and gas density fields across a
# range of radial pressure gradients for a single case of a dimensionless
# stopping time and total dust-to-gas mass ratio.
#
# Author: Stanley A. Baronett
# Created: 2022-07-25
# Updated: 2022-10-26
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import numpy as np

mpl.rcParams["axes.formatter.offset_threshold"] = 2
fig, axs = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(7, 4.5))
workdir = '../..'
case = 'AB'
Pis = ['0.01', '0.02', '0.05', '0.10']
res = '2048'
vmin_p, vmax_p = 1e-1, 1e1 # default for AB (CPDD > 85%; CPDD < 5%)
out = 100                  # default for AB (last out)

# Check for and override with user-passed arguments
if len(sys.argv) > 1:
    case = sys.argv[1]
    res = sys.argv[2]
    vmin_p = float(sys.argv[3])
    vmax_p = float(sys.argv[4])
    out = int(sys.argv[5])


for i, Pi in enumerate(Pis):
    # Collect parameters and plot densities
    path = f'{workdir}/{case}/{Pi}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    c_s = athinput['hydro']['iso_sound_speed']
    Omega = athinput['problem']['omega']
    H_g = c_s/Omega
    data = athena_read.athdf(f'{path}/athdf/SI.out1.{out:05}.athdf')
    xf, zf = data['x1f']/H_g, data['x2f']/H_g
    t = data['Time']
    clip = np.clip(data['rhop'][0], vmin_p, vmax_p)
    rhops = axs[0][i].pcolormesh(xf, zf, clip, norm=colors.LogNorm(),
                                 cmap='plasma')
    rhogs = axs[1][i].pcolormesh(xf, zf, data['rho'][0])

    # Add and format dust color bars, titles, and x-axis labels
    cb_rhop = fig.colorbar(rhops, ax=axs[0][i], location='top')
    axs[0][i].set_title(f'$\Pi={float(Pi)}$', pad=40)
    axs[0][i].set(aspect='equal')
    axs[1][i].set(xlabel=r'$x/H_\mathrm{g}$', aspect='equal')

    # Add and format gas color bars
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-1, 1))
    cb_rhog = fig.colorbar(rhogs, ax=axs[1][i], location='top',
                           format=formatter)

for ax in axs.flat:
    ax.label_outer()
    ax.minorticks_on()
    ax.set(aspect='equal')
    ax.tick_params(axis='both', which='both', top=True, right=True)
    ax.tick_params(axis='x', labelrotation=45)

# Format and save figure
axs[0][0].text(-0.65, 1.31, r'$\rho_\mathrm{p}/\rho_\mathrm{g,0}$', ha='left',
               va='top', transform=axs[0][0].transAxes)
axs[1][0].text(-0.65, 1.31, r'$\rho_\mathrm{g}/\rho_\mathrm{g,0}$', ha='left',
               va='top', transform=axs[1][0].transAxes)
axs[0][0].set(ylabel=r'$z/H_\mathrm{g}$')
axs[1][0].set(ylabel=r'$z/H_\mathrm{g}$')
plt.savefig(f'figs/{case}-{res}_snaps.png', dpi=1000, bbox_inches='tight',
            pad_inches=0.01)
