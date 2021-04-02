"""Plot max. & min., dust & gas, density growth over time.

Leave one blank line.  The rest of this docstring should contain an
overall description of the module or program.  Optionally, it may also
contain a brief description of exported classes and functions and/or usage
examples.

  Typical usage example:

  foo = ClassFoo()
  bar = foo.FunctionBar()
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Collect .athdf outputs, init sim consts. and grid
athinput = athena_read.athinput('../athinput.si.nas')
outputs = sorted(list(Path('../athdf').glob(athinput['job']['problem_id'] +
                                        '.out1.*.athdf')))
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
T = 2*np.pi/Omega                          # orbital period
epsilon = athinput['problem']['epsilon']   # avg. dust/gas ρ-ratio in BG state
times = []                                 # sim output times
rhopmax, rhopmin = [], []                  # max/min dust densities
rhogmax, rhogmin = [], []                  # max/min gas densities

for output in outputs:                     # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time'] / T)
    temp = np.amax(data['rhop']) - epsilon # difference w/ epsilon
    rhopmax.append(temp)
    temp = epsilon - np.amin(data['rhop'])
    rhopmin.append(temp)
    temp = np.amax(data['rho']) - 1
    rhogmax.append(temp)
    temp = 1 - np.amin(data['rho'])
    rhogmin.append(temp)

fig, axs = plt.subplots(2, 2, sharex='col', figsize=(17,7),
                        gridspec_kw={'hspace': 0})

# upper-left subplot
axs[0,0].set_title(r'Max. Densities', size='x-large')
axs[0,0].set_ylabel(r'$\rho_{p,max} - \epsilon$', size='large')
axs[0,0].semilogy(times, rhopmax)
# lower-left subplot
axs[1,0].set_ylabel(r'$\rho_{g,max} - 1$', size='large')
axs[1,0].set_xlabel('Orbital Period', size='large')
axs[1,0].semilogy(times, rhogmax)
# upper-right subplot
axs[0,1].set_title(r'Min. Densities', size='x-large')
axs[0,1].set_ylabel(r'$\epsilon - \rho_{p,min}$', size='large')
axs[0,1].semilogy(times, rhopmin)
# lower-left subplot
axs[1,1].set_ylabel(r'$1 - \rho_{g,min}$', size='large')
axs[1,1].set_xlabel('Orbital Period', size='large')
axs[1,1].semilogy(times, rhogmin)

for i,ax in enumerate(axs.flat):
    ax.grid()

plt.savefig('../plots/density_growth.pdf', bbox_inches='tight')
