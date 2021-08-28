"""Plot evolution of max. dust density.

"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Collect .athdf outputs, init sim consts. and grid
athinput = athena_read.athinput('../athinput.si')
outputs = sorted(list(Path('../athdf').glob(athinput['job']['problem_id'] +
                                        '.out2.*.athdf')))
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
T = 2*np.pi/Omega                          # orbital period
epsilon = athinput['problem']['epsilon']   # avg. dust/gas ρ-ratio in BG state
times = []                                 # sim output times
rhopmax = []                               # max dust densities

for output in outputs:                     # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time'] / T)
    temp = np.amax(data['rhop']) - epsilon # difference w/ epsilon
    rhopmax.append(temp)

fig, ax = plt.subplots(figsize=(17,7))

# upper-left subplot
ax.set_title(r'Max. Densities', size='x-large')
ax.set_ylabel(r'$\rho_{p,max} - \epsilon$', size='large')
ax.set_xlabel(r'$t$ / $T$', size='large')
ax.semilogy(times, rhopmax)
ax.grid()

plt.savefig('../plots/rhop_evo.pdf', bbox_inches='tight')
