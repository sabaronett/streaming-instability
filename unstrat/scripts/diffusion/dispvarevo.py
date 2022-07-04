#!/usr/bin/env python3
#==============================================================================
# dispvarevo.py
#
# Plot the variance of the particle displacement distribution over time and
# overplot the diffusion coefficient fit for comparison.
#
# Author: Stanley A. Baronett
# Created: 2022-07-04
# Updated: 2022-07-04
#==============================================================================
import sys
sys.path.insert(0, '/home/stanley/github/sabaronett/local/athena/athena-dust/vis/python')
import athena_read
import numpy as np
import matplotlib.pyplot as plt

t0 = float(sys.argv[1])
dtmax = int(sys.argv[2])
tfit = np.asarray([0, dtmax])
dirs = ['x', 'z']
fig, axs = plt.subplots(2, 1, figsize=(6, 8), dpi=150)

athinput = athena_read.athinput('athinput.si')
Pi = athinput['problem']['duy0']
taus = athinput['particles']['taus0']*athinput['problem']['omega']
epsilon = athinput['problem']['epsilon']
res = athinput['mesh']['nx1']
disp = np.load('dat/disp.npz')
dcoeffs = np.load(f'output/dcoeff-{str(dtmax)}.npz')
sat = 'sat'
if round(taus, 1) == 0.1 and epsilon == 1.0: case = 'AB'
elif round(taus, 1) == 1.0 and epsilon == 0.2: case = 'BA'
else: case = '?'
title = f'{case}, $\Pi={Pi}$ ({res}$^2$), $t_0=t_\mathrm{{{sat}}}={t0:.1f}\,T$'

axs[0].set(title=title)
axs[1].set(xlabel='$\Delta t\,/\,T$')

for i, dir in enumerate(dirs):
    if dir == 'x':
        vars = np.var(disp['dxp'], axis=1)
        dcoeff = dcoeffs['dpx'][0]
    else:
        vars = np.var(disp['dyp'], axis=1)
        dcoeff = dcoeffs['dpy'][0]
    
    axs[i].grid()
    axs[i].plot(disp['t'], vars, label=f'$\sigma_{{{dir}}}^2(t)$')
    axs[i].plot(tfit, 2*dcoeff*tfit, '--', label=f'$2D_{{{dir}}}$')
    axs[i].legend()
    axs[i].set(ylabel=f'$\sigma_{{{dir}}}^2$')

fig.subplots_adjust(hspace=0)
plt.savefig('output/dispvarevo.pdf', bbox_inches='tight', pad_inches=0.01)
