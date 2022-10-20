#!/usr/bin/env python3
#==============================================================================
# AB_avgRs_rad-prof.py
#
# Plot with logarithmically-spaced bins azimuthally-averaged radial profiles of
# the time-averaged normalized autocorrelations of snapshots of the dust and
# gas density fields across a range of radial pressure gradients for case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-09
# Updated: 2022-10-20
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import fftpack, stats

def norms(xv, zv, pole):
    rv = np.zeros((len(zv), len(xv)))
    for i, z in enumerate(zv):
        for j, x in enumerate(xv):
            rv[i][j] = np.sqrt((x - pole[0])**2 + (z - pole[1])**2)
    return rv

fig, axs = plt.subplots(2, sharex=True, figsize=(3.15, 4))
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
    # Collect parameters
    print(f'{case}/{Pi[0]}: Processing...', flush=True)
    path = f'{workdir}/{case}/{Pi[0]}/{res}'
    athinput = athena_read.athinput(f'{path}/athinput.si')
    outputs = sorted(list(Path(f'{path}/athdf').glob(
        athinput['job']['problem_id']+'.out1.*.athdf')))
    dt = athinput['output1']['dt']
    i_sat  = t_sat//dt
    outputs = outputs[i_sat:]
    c_s = athinput['hydro']['iso_sound_speed']
    etar = float(Pi[0])*c_s
    data = athena_read.athdf(outputs[0])
    xv, zv = data['x1v'], data['x2v']
    left, right = xv[-1]/res, 1.32*xv[-1]
    x0, z0 = len(xv)//2, len(zv)//2
    pole = (xv[x0], zv[z0])
    rv = norms(xv, zv, pole).ravel()
    indices = np.where(rv > xv[-1])[0]
    indices = np.append(indices, np.where(rv == 0)[0])
    rv = np.delete(rv, indices)
    r0 = xv[x0]
    base = np.sqrt(2)
    leftmost_edge = r0/np.sqrt(base)
    num = int(np.sqrt(res))
    bin_edges = leftmost_edge*np.logspace(0, num, num=(num + 1), base=base)
    rvs = np.empty((len(outputs), res, res))
    Rps = np.empty((len(outputs), res, res))
    Rgs = np.empty((len(outputs), res, res))

    for j, output in enumerate(outputs):
        data = athena_read.athdf(output)
        rvs[j] = rv
        # Process dust
        ft = fftpack.fft2(data['rhop'][0])
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        log = np.log10(shift).ravel()
        Rps[j] = np.delete(log, indices)
        # Process gas
        ft = fftpack.fft2(data['rho'][0])
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        offset = (shift.ravel() - 1)*1e11
        Rgs[j] = np.delete(offset, indices)
        print(f'\t{j/len(outputs):.0%}', flush=True)

    # Bin dust
    dust_means, bin_edges, binnumnber = stats.binned_statistic(rvs.ravel(),
        Rps.ravel(), statistic='mean', bins=bin_edges)
    dust_stds, bin_edges, binnumnber = stats.binned_statistic(rvs.ravel(),
        Rps.ravel(), statistic='std', bins=bin_edges)
    dust_highs = dust_means + dust_stds
    dust_lows = dust_means - dust_stds
    # Bin gas
    gas_means, bin_edges, binnumnber = stats.binned_statistic(rvs.ravel(),
        Rgs.ravel(), statistic='mean', bins=bin_edges)
    gas_stds, bin_edges, binnumnber = stats.binned_statistic(rvs.ravel(),
        Rgs.ravel(), statistic='std', bins=bin_edges)
    gas_highs = gas_means + gas_stds
    gas_lows = gas_means - gas_stds
    # Plot results
    axs[0].stairs(dust_means, bin_edges, baseline=float('-inf'),
                  color=Pi[1], lw=1.5, label=Pi[0])
    axs[0].stairs(dust_highs, bin_edges, baseline=dust_lows, fill=True,
                color=Pi[1], alpha=0.2)
    axs[1].stairs(gas_means, bin_edges, baseline=float('-inf'),
                  color=Pi[1], lw=1.5)
    axs[1].stairs(gas_highs, bin_edges, baseline=gas_lows, fill=True,
                  color=Pi[1], alpha=0.2)
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].legend(title=r'$\Pi$')
axs[0].set(ylabel=r'$\log\left(\int_0^{2\pi}\mathrm{R}_{\rho_\mathrm{p}\rho_\mathrm{p}}\mathrm{d}\phi\right)$')
axs[1].set(xlim=(left, right), xscale='log', yscale='symlog',
           ylabel=r'$\left(\int_0^{2\pi}\mathrm{R}_{\rho_\mathrm{g}\rho_\mathrm{g}}\mathrm{d}\phi-1\right)\times10^{11}$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_avgRs_rad-prof.pdf', bbox_inches='tight',
            pad_inches=0.01)
