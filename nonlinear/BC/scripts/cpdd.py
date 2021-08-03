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
rhops = []                                   # dust densities

# Load & process saturated-state data into memory
i_sat = 1000                                 # 1st index of sat. state
sat_outputs = outputs[i_sat:]
for output in sat_outputs:
    data = athena_read.athdf(output)
    temp = data['rhop'].flatten() / epsilon  # flatten & convert
    rhops.append(np.sort(temp))              # sort

# Find min., max., avg. of each ordered rhop over saturated state
min_rhops = np.amin(rhops, axis=0)
max_rhops = np.amax(rhops, axis=0)
avg_rhops = np.average(rhops, axis=0)

# Find std. dev. computed in log space
lnrhops = np.log(rhops)
sigmas = np.exp(np.std(lnrhops, axis=0))
cdf = np.linspace(1, 0 , min_rhops.size, endpoint=False)

# CPDD
fig, ax = plt.subplots(figsize=(5,5))
ax.set_title('Cumulative Particle-Density Distributions \n'\
             +r'($\tau_s={:.1f},\,\epsilon={:.1f}$)'
             .format(tau_s, epsilon), size='x-large')
ax.set_xlabel(r'$\rho_p$ / $\langle \rho_p \rangle$', size='large')
ax.set_ylabel(r'P$(>\rho_p)$', size='large')
ax.loglog(avg_rhops, cdf, label=r'$\mu$')
ax.fill_betweenx(cdf, min_rhops, max_rhops, alpha=0.2, label='[min., max.]')
ax.fill_betweenx(cdf, avg_rhops/sigmas, avg_rhops*sigmas, alpha=0.4,
                 label=r'$[\sigma^{-1}\mu,\,\sigma\mu]$')
ax.set_xlim(0.1, 1000)
ax.set_ylim(1e-5, 1)
ax.legend(loc='lower left')
ax.grid()

# Save figure and plotting data
makesubdir('../plots') # create file output directory
plt.savefig('../plots/cpdd.pdf', bbox_inches='tight',
            pad_inches=0.01)
writetxt(min_rhops, cdf, '../plots/cpdd_min.txt')
writetxt(max_rhops, cdf, '../plots/cpdd_max.txt')
writetxt(avg_rhops, cdf, '../plots/cpdd_avg.txt')
writetxt(sigmas, cdf, '../plots/cpdd_std.txt')
