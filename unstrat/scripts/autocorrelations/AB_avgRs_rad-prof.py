#!/usr/bin/env python3
#==============================================================================
# AB_avgRs_rad-prof.py
#
# Plot azimuthally-averaged radial profiles of the time-averaged normalized
# autocorrelations of snapshots of the dust and gas density fields across a
# range of radial pressure gradients for case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-09
# Updated: 2022-10-09
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from pathlib import Path
from scipy import fftpack

def norms (xf, zf):
    rf = np.empty((len(zf), len(xf)))
    for i,z in enumerate(zf):
        for j,x in enumerate(zf):
            rf[i][j] = np.sqrt(x**2 + z**2)
    return rf

def rad_prof(data, center):
    """
    https://stackoverflow.com/a/21242776/13149558
    """
    y, x = np.indices(data.shape)
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(int)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    profile = tbin / nr
    return profile

fig, axs = plt.subplots(2, sharex=True, figsize=(3.15, 4), dpi=150)
workdir = '../..'
case = 'AB'
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = 2048
t_sat = 4 # [T]

# Check for and override with user-passed arguments
if len(sys.argv) > 1:
    res = int(sys.argv[1])
    t_sat = float(sys.argv[2])

for i, Pi in enumerate(Pis):
    # Collect parameters and plot densities
    print(f'{case}/{Pi[0]}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi[0]}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    outputs = sorted(list(Path(f'{path}/athdf').glob(\
        athinput['job']['problem_id']+'.out1.*.athdf')))
    dt = athinput['output1']['dt']
    i_sat  = int(t_sat/dt)
    outputs = outputs[i_sat:]
    c_s = athinput['hydro']['iso_sound_speed']
    etar = float(Pi[0])*c_s
    data = athena_read.athdf(outputs[0])
    xf, zf = data['x1f']/etar, data['x2f']/etar
    center = int(len(xf)/2)
    xf, zf = np.delete(xf, center), np.delete(zf, center)
    rf = norms(xf, zf)
    rs = rad_prof(rf, [center, center])
    Rps = np.empty((len(outputs), res, res))
    Rgs = np.empty((len(outputs), res, res))

    for j, output in enumerate(outputs):
        data = athena_read.athdf(output)
        ft = fftpack.fft2(data['rhop'][0])
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        log = np.log10(shift)
        Rps[j] = log
        ft = fftpack.fft2(data['rho'][0])
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        offset = (shift - 1)*1e8
        Rgs[j] = offset
        print(f'\t{j/len(outputs):.0%}', flush=True)

    avgRp = np.average(Rps, axis=0)
    avgRg = np.average(Rgs, axis=0)
    Rp_prof = rad_prof(avgRp, [center, center])
    Rg_prof = rad_prof(avgRg, [center, center])
    axs[0].plot(rs, Rp_prof, color=Pi[1], label=Pi[0])
    axs[1].plot(rs, Rg_prof, color=Pi[1])
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(title=r'$\Pi$')
axs[0].set(ylabel=r'$\log(\int\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}\mathrm{d}r)$',
           title='Radial Profile')
axs[1].set(xlabel=r'$x/(\eta r)$', xscale='log',
           ylabel=r'$\int\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}\mathrm{d}r\times10^{-8}+1$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_autocorrelations.png', dpi=1000,
            bbox_inches='tight', pad_inches=0.01)
