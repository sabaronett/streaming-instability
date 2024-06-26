#!/usr/bin/env python3
#==============================================================================
# AB_avgRs_rad-prof.py
#
# Plot azimuthally-averaged radial profiles, logarithmically binned and scaled
# by Pi, of the time-averaged normalized autocorrelations of snapshots of the
# dust and gas density fields across a range of radial pressure gradients for
# case AB.
#
# Author: Stanley A. Baronett
# Created: 2022-10-09
# Updated: 2022-12-15
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

fig, axs = plt.subplots(2, sharex=True, figsize=(4.45, 4.32))
workdir = '../..'
case = 'AB'
Pis = [['0.01', 'tab:red'], ['0.02', 'tab:orange'],
       ['0.05', 'tab:green'], ['0.10', 'tab:blue']]
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
    outputs = sorted(list(Path(f'{path}/athdf').glob(
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
    Rps = np.empty((len(outputs), len(rv)))
    Rgs = np.empty((len(outputs), len(rv)))

    for j, output in enumerate(outputs):
        data = athena_read.athdf(output)

        # Process dust
        diff = data['rhop'][0] - epsilon
        ft = fftpack.fft2(diff)
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        Rps[j] = np.delete(shift, indices).ravel()
        
        # Process gas
        diff = data['rho'][0] - 1
        ft = fftpack.fft2(diff)
        ac = fftpack.ifft2(ft*np.conjugate(ft)).real
        norm = ac/ac[0][0]
        shift = fftpack.fftshift(norm)
        Rgs[j] = np.delete(shift, indices).ravel()
        print(f'\t{j/len(outputs):.0%}', flush=True)

    # Bin dust
    avgRp = np.average(Rps, axis=0)
    dust_means, bin_edges, binnumnber = stats.binned_statistic(rv, avgRp,
        statistic='mean', bins=bin_edges)
    dust_stds, bin_edges, binnumnber = stats.binned_statistic(rv, avgRp,
        statistic='std', bins=bin_edges)
    dust_highs = dust_means + dust_stds
    dust_lows = dust_means - dust_stds

    # Bin gas
    avgRg = np.average(Rgs, axis=0)
    gas_means, bin_edges, binnumnber = stats.binned_statistic(rv, avgRg,
        statistic='mean', bins=bin_edges)
    gas_stds, bin_edges, binnumnber = stats.binned_statistic(rv, avgRg,
        statistic='std', bins=bin_edges)
    gas_highs = gas_means + gas_stds
    gas_lows = gas_means - gas_stds

    # Plot histograms
    axs[0].stairs(dust_means, bin_edges/float(Pi[0]), baseline=float('-inf'),
                  color=Pi[1], lw=1.5)
    axs[0].stairs(dust_highs, bin_edges/float(Pi[0]), baseline=dust_lows,
                  fill=True, color=Pi[1], alpha=0.2)
    axs[1].stairs(gas_means, bin_edges/float(Pi[0]), baseline=float('-inf'),
                  color=Pi[1], lw=1.5, label=Pi[0])
    axs[1].stairs(gas_highs, bin_edges/float(Pi[0]), baseline=gas_lows,
                  fill=True, color=Pi[1], alpha=0.2)
    print(f'\tdone.', flush=True)

for ax in axs.flat:
    ax.grid()
    ax.label_outer()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True, right=True)

# Format and save figure
axs[0].set(ylim=(-0.05, 1.05), ylabel=r'$\mathcal{R}_\mathrm{p}$')
axs[1].legend(loc='lower left', title=r'$\Pi$')
axs[1].set(ylim=(-0.05, 1.05), xscale='log', xlabel=r'$r/(\Pi H_\mathrm{g})$',
           ylabel=r'$\mathcal{R}_\mathrm{g}$')
plt.subplots_adjust(hspace=0)
plt.savefig(f'figs/{case}_avgRs_rad-prof.pdf', bbox_inches='tight',
            pad_inches=0.01)
