#!/usr/bin/env python3
#==============================================================================
# AB_power_spectrum.py
#
# Plot the power spectrum of the normalized autocorrelations of snapshots of
# the dust and gas density fields across a range of radial pressure gradients
# for case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-10
# Updated: 2022-10-10
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy import fftpack

def norms (xf, zf):
    size = int(len(zf)*len(xf))
    rf = np.empty(size)
    for i,z in enumerate(zf):
        for j,x in enumerate(zf):
            rf[i+j] = np.sqrt(x**2 + z**2)
    return rf

fig, axs = plt.subplots(2, sharex=True, figsize=(3.15, 4))
workdir = '../..'
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
    xf, zf = data['x1f']/etar, data['x2f']/etar
    center = int(len(xf)/2)
    xf, zf = np.delete(xf, center), np.delete(zf, center)
    rf = norms(xf, zf)
    ft = fftpack.fft2(data['rhop'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    ravel = shift.ravel()
    log = np.log10(ravel)
    axs[0].scatter(rf, log, s=0.1, color=Pi[1], label=Pi[0])

    ft = fftpack.fft2(data['rho'][0])
    ac = fftpack.ifft2(ft*np.conjugate(ft)).real
    norm = ac/ac[0][0]
    shift = fftpack.fftshift(norm)
    ravel = shift.ravel()
    offset = (ravel - 1)*1e8
    axs[1].scatter(rf, offset, s=0.1, color=Pi[1])
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(title=r'$\Pi$')
axs[0].set(ylabel=r'$\log(\int\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}\mathrm{d}r)$',
           title='Power Spectrum')
axs[1].set(xlabel=r'$x/(\eta r)$', xscale='log',
           ylabel=r'$\int\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}\mathrm{d}r\times10^{-8}+1$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_Rs_power-spectrum.png', dpi=1000,
            bbox_inches='tight', pad_inches=0.01)
