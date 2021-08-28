"""Cumulative Particle-Density Distribution.

Computes, plots, and outputs the time-averaged cumulative particle-density
distribution of the SI run during the saturated, turbulent state, including
time-varying minima, maxima, and standard deviations in logarithmic space.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

def writetxt(x, y, path='data.txt'):
    """
    Function to output series data to a two-column text file.
        
    Parameters
    ----------
    x : numpy.ndarray
        Array to be written out to the first column of the data file.
    y : numpy.ndarray
        Array to be written out to the second column of the data file.
    path : str
        Path and filename of the data file to be outputted.
        Default path set to working directory and filename "data.txt"
    """
    with open(path, 'w') as f: # will overwrite existing file
        for i in range(x.size):
            f.write('{:.16e}\t{:.16e}\n'.format(x[i], y[i]))

def makesubdir(name):
    if not os.path.exists(name):
        os.makedirs(name)

# Collect .athdf inputs, outputs; init sim consts. & grid
athinput = athena_read.athinput('../athinput.si')
Omega = athinput['problem']['omega']         # local Keplerian ang. freq.
tau_s = athinput['particles']['taus0']*Omega # dimensionless stopping time
epsilon = athinput['problem']['epsilon']     # avg. BG dust/gas ρ-ratio
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                            '.out2.*.athdf')))
rhopslist = []                               # dust densities

# Load & process saturated-state data into memory
i_sat = 375                                  # 1st index of sat. state
sat_outputs = outputs[i_sat:]
for output in sat_outputs:
    data = athena_read.athdf(output)
    rhopslist.append(data['rhop'] / epsilon) # convert to ρ_p / ρ_p0
rhops = np.asarray(rhopslist)                # no averaging
flat = np.ravel(rhops)
rhops = np.sort(flat)
# define prob. w/ total # cells & total # snapshots
cdf = np.linspace(1, 0 , rhops.size, endpoint=False)

# CPDD
fig, ax = plt.subplots(figsize=(5,5))
ax.set_title('Cumulative Particle-Density Distributions \n'\
             +r'($\tau_s={:.1f},\,\epsilon={:.1f}$)'
             .format(tau_s, epsilon), size='x-large')
ax.set_xlabel(r'$\rho_p$ / $\langle \rho_p \rangle$', size='large')
ax.set_ylabel(r'P$(>\rho_p)$', size='large')
ax.loglog(rhops, cdf)
ax.set_xlim(0.01, 1000)
ax.set_ylim(1e-10, 1)
ax.grid()

# Save plotting data
makesubdir('../plots') # create file output directory
plt.savefig('../plots/cpdd4.pdf', bbox_inches='tight', pad_inches=0.01)
writetxt(rhops, cdf, '../plots/cpdd4.txt')
