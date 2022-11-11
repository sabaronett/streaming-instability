#!/usr/bin/env python3
#==============================================================================
# AB_avgRs-z0-x0.py
#
# Plot 1D cuts, at z = 0 and x = 0, of time-averaged, normalized auto-
# correlations of snapshots of the dust and gas density fields across a range
# of radial pressure gradients for case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-30
# Updated: 2022-11-11
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import fftpack

fig, axs = plt.subplots(2, sharex=True, figsize=(4.45, 4.32))
workdir = '../../..'
case = 'AB'
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = 2048
t_sat = 5 # [T]

# Check for and override with user-passed arguments
if len(sys.argv) > 1:
    res = int(sys.argv[1])
    t_sat = float(sys.argv[2])

for i, Pi in enumerate(Pis):
    # Collect parameters
    print(f'{case}/{Pi[0]}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi[0]}/{res}'
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
    x0, z0 = len(xv)//2, len(zv)//2
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
        Rps[j] = shift

        # Process gas
        diff = data['rho'][0] - 1
        ft = fftpack.fft2(diff)
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        Rgs[j] = shift
        print(f'\t{j/len(outputs):.0%}', flush=True)

    # Compute avg., std dev., and plot 1D slices
    avgRp = np.average(Rps, axis=0)
    stdRp = np.std(Rps, axis=0)
    avgRg = np.average(Rgs, axis=0)
    stdRg = np.std(Rgs, axis=0)
    
    axs[0].semilogx(xv, avgRp[z0], color=Pi[1])
    axs[0].fill_betweenx(xv, avgRp[z0]/stdRp[z0], avgRp[z0]*stdRp[z0],
                         color=Pi[1], ec=None, alpha=0.2)
    axs[0].semilogx(xv, avgRp[:, x0], color=Pi[1], ls='--')
    axs[0].fill_betweenx(xv, avgRp[:, x0]/stdRp[:, x0],
                         avgRp[:, x0]*stdRp[:, x0], color=Pi[1], ec=None,
                         alpha=0.2)
    axs[1].semilogx(xv, avgRg[x0], color=Pi[1], label=Pi[0])
    axs[1].fill_betweenx(xv, avgRg[z0]/stdRg[z0], avgRg[z0]*stdRg[z0],
                         color=Pi[1], ec=None, alpha=0.2)
    axs[1].semilogx(xv, avgRg[:, z0], color=Pi[1], ls='--')
    axs[1].fill_betweenx(xv, avgRg[:, x0]/stdRg[:, x0],
                         avgRg[:, x0]*stdRg[:, x0], color=Pi[1], ec=None,
                         alpha=0.2)
    print(f'\tdone.', flush=True)

    # Plot ghost points for colorless line style and add legends
    ls_dust, = axs[0].semilogx([], [], color='tab:gray',
                               label=r'$\mathrm{R}_{\rho\rho}(z=0$)')
    ls_gas,  = axs[0].semilogx([], [], color='tab:gray', ls='--',
                               label=r'$\mathrm{R}_{\rho\rho}(x=0)$')
    axs[0].legend(handles=[ls_dust, ls_gas], loc='upper right')
    axs[1].legend(loc='lower left', title=r'$\Pi$')
    
for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].set(ylim=(-0.05, 1.05),
           ylabel=r'$\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}$')
axs[1].set(ylim=(-0.05, 1.05), xscale='log', xlabel=r'$(x,z)/H_\mathrm{g}$',
           ylabel=r'$\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_avgRs-z0-x0.pdf', bbox_inches='tight',
            pad_inches=0.01)
