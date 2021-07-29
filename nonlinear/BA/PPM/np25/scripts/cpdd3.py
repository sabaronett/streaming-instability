"""Cumulative Particle-Density Distribution.

Computes and plots the time-averaged cumulative particle-density distribution
of the SI run during the saturated, turbulent state.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import numpy as np
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

# Collect .athdf inputs, outputs; init sim consts. & grid
athinput = athena_read.athinput('../athinput.si')
nx1 = athinput['mesh']['nx1']                # num. radial zones
nx2 = athinput['mesh']['nx2']                # num. vertical zones
nx3 = athinput['mesh']['nx3']                # num. azimuthal zones
zones = nx1*nx2*nx3                          # total num. of zones
c_s = athinput['hydro']['iso_sound_speed']   # sound speed
Omega = athinput['problem']['omega']         # local Keplerian angular frequency
tau_s = athinput['particles']['taus0']*Omega # dimensionless stopping time
epsilon = athinput['problem']['epsilon']     # avg. dust/gas ρ-ratio in BG state
Np_tot = athinput['problem']['npx1']\
    *athinput['problem']['npx2']\
    *athinput['problem']['npx3']             # total number of particles
Np = Np_tot/nx1/nx2/nx3                      # theo avg num particles per cell
H = c_s / Omega                              # gas scale height
T = 2*np.pi/Omega                            # orbital period
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                        '.out2.*.athdf')))
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f'] / H, data['x2f'] / H
rhops = []                                   # dust densities

# Load & process saturated-state data into memory
i_sat = 375                                  # 1st index of sat. state
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
mus = np.average(lnrhops, axis=0)
sigmas = np.std(lnrhops, axis=0)
logstd_rhops = np.exp(mus + 0.5*np.square(sigmas))\
    *np.sqrt(np.exp(np.square(sigmas)) - 1)
cdf = np.linspace(1, 0 , min_rhops.size, endpoint=False)

# CPDD
fig, ax = plt.subplots(figsize=(6,5))
# ax.set_aspect('equal')
ax.set_title('Cumulative Particle-Density Distributions', size='x-large')
ax.set_xlabel(r'$\rho_p$ / $\langle \rho_p \rangle$', size='large')
ax.set_ylabel(r'P$(>\rho_p)$', size='large')
ax.loglog(avg_rhops, cdf, label='Average')
ax.fill_betweenx(cdf, min_rhops, max_rhops, alpha=0.2, label='Min./Max.')
ax.fill_betweenx(cdf, avg_rhops-logstd_rhops, avg_rhops+logstd_rhops, alpha=0.4,
                 label=r'$\pm\sigma$')
ax.plot([1e-2, 1e3], [1e-1, 1e-1], '--', color='black')
ax.set_xlim(0.1, 1000)
ax.set_ylim(1e-5, 1)
ax.legend(loc='lower left')
ax.grid()

# Save figure and plotting data
plt.savefig('../plots/CPDD2_np{:.0f}.pdf'.format(Np), bbox_inches='tight',
            pad_inches=0.01)
writetxt(min_rhops, cdf, '../plots/CPDD_np{:.0f}_min.txt'.format(Np))
writetxt(max_rhops, cdf, '../plots/CPDD_np{:.0f}_max.txt'.format(Np))
writetxt(avg_rhops, cdf, '../plots/CPDD_np{:.0f}_avg.txt'.format(Np))
writetxt(logstd_rhops, cdf, '../plots/CPDD_np{:.0f}_logstd.txt'.format(Np))
