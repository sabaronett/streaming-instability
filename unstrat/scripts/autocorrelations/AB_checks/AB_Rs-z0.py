#!/usr/bin/env python3
#==============================================================================
# AB_Rs-z0.py
#
# Plot 1D cuts (at z=0) of normalized autocorrelations of snapshots of the dust
# and gas density fields across a range of radial pressure gradients for case
# AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-10
# Updated: 2022-10-13
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

fig, axs = plt.subplots(2, sharex=True, figsize=(3.15, 4), dpi=300)
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
    etar = float(Pi[0])*c_s
    data = athena_read.athdf(f'{path}/athdf/SI.out1.00100.athdf')
    xv, zv = data['x1v'], data['x2v']
    z0 = int(len(zv)/2)
    ft = fftpack.fft2(data['rhop'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    log = np.log10(shift)
    ft = fftpack.fft2(data['rho'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    offset = (shift-1)*1e9

    axs[0].semilogx(xv/etar, log[z0], color=Pi[1], label=Pi[0])
    axs[1].semilogx(xv/etar, offset[z0], color=Pi[1])
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(loc='upper right', title=r'$\Pi$')
axs[0].set(ylabel=r'$\log\left[\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}(z=0)\right]$')
axs[1].set(yscale='symlog', xlabel=r'$x/(\eta r)$',
           ylabel=r'$\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}(z=0)\times10^{-9}+1$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_Rs-z0.pdf', bbox_inches='tight', pad_inches=0.01)
