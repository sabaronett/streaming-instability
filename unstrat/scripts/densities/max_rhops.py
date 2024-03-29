#!/usr/bin/env python3
#==============================================================================
# max_rhop.py
#
# Plot the maximum dust density evolution.
#
# Author: Stanley A. Baronett
# Created: 2022-09-27
# Updated: 2022-11-03
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
from pathlib import Path

fig, axs = plt.subplots(2, 1, figsize=(7, 4.5))
workdir = '../..'
cases = ['AB', 'BA']
Pis = [['0.01', 'tab:blue'], ['0.02', 'tab:green'],
       ['0.05', 'tab:orange'], ['0.10', 'tab:red']]
res = '2048'

for i, case in enumerate(cases):
    for Pi in Pis:
        path = f'{workdir}/{case}/{Pi[0]}/{res}'
        athinput = athena_read.athinput(f'{path}/athinput.si')
        outputs = sorted(list(Path(f'{path}/athdf').glob(\
            athinput['job']['problem_id']+'.out1.*.athdf')))
        times, max_rhops = [], []

        for output in outputs:
            athdf = athena_read.athdf(output)
            times.append(athdf['Time'])
            max_rhop = athdf['rhop'].max()
            max_rhops.append(max_rhop)

        axs[i].semilogy(times, max_rhops, color=Pi[1], label=Pi[0])
        print(f'{case}/{Pi[0]} done.', flush=True)

# Format subplots
for ax in axs.flat:
    ax.grid()
    ax.minorticks_on()
    ax.set(ylabel=r'$\max{\rho_\mathrm{p}}/\rho_\mathrm{g,0}$')
    ax.tick_params(which='both', top=True, right=True)

axs[0].legend(loc='lower right', title=r'$\Pi$')
axs[0].set(title='AB')
axs[1].set(title='BA')
axs[1].set(xlabel=r'$t/T$')

plt.savefig('figs/max_rhops.pdf', bbox_inches='tight', pad_inches=0.01)
