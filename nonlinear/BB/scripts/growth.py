"""Plot maximum particle density evolution.

Plots and outputs the maximum particle density as a function of time.
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
    x : array_like
        Array to be written out to the first column of the data file.
    y : array_like
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

# Collect .athdf outputs, init sim consts.
athinput = athena_read.athinput('../athinput.si')
Omega = athinput['problem']['omega']         # local Keplerian ang. freq.
tau_s = athinput['particles']['taus0']*Omega # dimensionless stopping time
epsilon = athinput['problem']['epsilon']     # avg. BG dust/gas ρ-ratio
T = 2*np.pi/Omega                            # orbital period
outputs = sorted(list(Path('../athdf').glob(athinput['job']['problem_id'] +
                                            '.out2.*.athdf')))
times, rhopmax = [], []                      # sim out times, max dust dens.

for output in outputs:                       # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time'] / T)
    temp = np.amax(data['rhop'])
    rhopmax.append(temp)

# Plot
fig, ax = plt.subplots(figsize=(6,5))
ax.set_title('Maximum Particle Density Evolution \n'\
             +r'($\tau_s={:.1f},\,\epsilon={:.1f}$)'
             .format(tau_s, epsilon), size='x-large')
ax.set_ylabel(r'$\rho_{p,max} / \rho_{g0}$', size='large')
ax.set_xlabel(r'$t$ / $T$', size='large')
ax.semilogy(times, rhopmax)
ax.grid()

# Save figure and plotting data
makesubdir('../plots') # create file output directory
plt.savefig('../plots/growth.pdf', bbox_inches='tight',
            pad_inches=0.01)
writetxt(times, rhopmax, '../plots/growth.txt')
