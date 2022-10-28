#!/usr/bin/env python3
#==============================================================================
# AB_Rs_pow-spect_gas.py
#
# Plot the power spectrum of the normalized autocorrelations of snapshots of
# the dust and gas density fields across a range of radial pressure gradients
# for case AB.
# NOTE: In this version, no offset is applied to the gas autocorrelation.
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
from scipy import fftpack

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
    # Collect parameters and plot densities
    print(f'{case}/{Pi[0]}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi[0]}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    c_s = athinput['hydro']['iso_sound_speed']
    # etar = float(Pi[0])*c_s
    data = athena_read.athdf(f'{path}/athdf/SI.out1.00100.athdf')
    xv, zv = data['x1v'], data['x2v']
    left, right = xv[-1]/res, np.sqrt(8)*xv[-1]
    x0, z0 = int(len(xv)/2), int(len(zv)/2)
    pole = (xv[x0], zv[z0])
    rv = norms(xv, zv, pole)

    # Process dust
    ft = fftpack.fft2(data['rhop'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    log = np.log10(shift)
    # Process gas
    ft = fftpack.fft2(data['rho'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)

    # Plot spectrum
    axs[0].scatter(rv, log, s=0.1, color=Pi[1], label=Pi[0])
    axs[1].scatter(rv, shift, s=0.1, color=Pi[1])
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(loc='upper right', title=r'$\Pi$')
axs[0].set(ylabel=r'$\log\mathcal{R}_\mathrm{p}$')
axs[1].set(xlim=(left, right), xscale='log', xlabel=r'$r/H_\mathrm{g}$',
           ylabel=r'$\mathcal{R}_\mathrm{g}$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_Rs_pow-spect_gas.png', dpi=1000, bbox_inches='tight',
            pad_inches=0.01)
