#!/usr/bin/env python3
#==============================================================================
# BA_avgRs.py
#
# Plot time-averaged normalized autocorrelations of snapshots of the dust and
# gas density fields across a range of radial pressure gradients for case BA.
#
# Author: Stanley A. Baronett
# Created: 2022-09-28
# Updated: 2022-11-08
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path
from scipy import fftpack

mpl.rcParams["axes.formatter.offset_threshold"] = 2
fig, axs = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(9.32, 5.6))
workdir = '../..'
case = 'BA'
Pis = ['0.01', '0.02', '0.05', '0.10']
res = 2048
t_sat = 150 # [T]

# Check for and override with user-passed arguments
if len(sys.argv) > 1:
    res = int(sys.argv[1])
    t_sat = float(sys.argv[2])

for i, Pi in enumerate(Pis):
    # Collect parameters
    print(f'{case}/{Pi}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    outputs = sorted(list(Path(f'{path}/athdf').glob(\
        athinput['job']['problem_id']+'.out1.*.athdf')))
    dt = athinput['output1']['dt']
    i_sat  = int(t_sat/dt)
    outputs = outputs[i_sat:]
    c_s = athinput['hydro']['iso_sound_speed']
    Omega = athinput['problem']['omega']
    epsilon = athinput['problem']['epsilon']
    H_g = c_s/Omega
    data = athena_read.athdf(outputs[0])
    xv, zv = data['x1v']/H_g, data['x2v']/H_g
    Rps = np.empty((len(outputs), res, res))
    Rgs = np.empty((len(outputs), res, res))

    for j, output in enumerate(outputs):
        data = athena_read.athdf(output)

        # Process dust
        diff = data['rhop'][0] - epsilon
        ft = fftpack.fft2(diff)
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        clip = np.clip(shift, float('-inf'), 0.1)
        Rps[j] = clip

        # Process gas
        diff = data['rho'][0] - 1
        ft = fftpack.fft2(diff)
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        Rgs[j] = shift
        print(f'\t{j/len(outputs):.0%}', flush=True)

    avgRp = np.average(Rps, axis=0)
    avgRg = np.average(Rgs, axis=0)
    mesh_p = axs[0][i].pcolormesh(xv, zv, avgRp, cmap='plasma')
    mesh_g = axs[1][i].pcolormesh(xv, zv, avgRg)

    # Add and format color bars, titles, and x-axis labels
    cb_rhop = fig.colorbar(mesh_p, ax=axs[0][i], location='top')
    cb_rhog = fig.colorbar(mesh_g, ax=axs[1][i], location='top')
    axs[0][i].set_title(f'$\Pi={float(Pi)}$', pad=42)
    axs[0][i].set(aspect='equal')
    axs[1][i].set(xlabel=r'$x/H_\mathrm{g}$', aspect='equal')
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.label_outer()
    ax.minorticks_on()
    ax.set(aspect='equal')
    ax.tick_params(axis='both', which='both', top=True, right=True)
    ax.tick_params(axis='x', labelrotation=45)

# Format and save figure
axs[0][0].text(-0.42, 1.27,
               r'$\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}$',
               ha='left', va='top', transform=axs[0][0].transAxes)
axs[1][0].text(-0.42, 1.27,
               r'$\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}$',
               ha='left', va='top', transform=axs[1][0].transAxes)
axs[0][0].set(ylabel=r'$z/H_\mathrm{g}$')
axs[1][0].set(ylabel=r'$z/H_\mathrm{g}$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_avgRs.png', dpi=800, bbox_inches='tight',
            pad_inches=0.01)
