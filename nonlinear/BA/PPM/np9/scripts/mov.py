"""Animate Dust Density.

Compiles and saves an MPEG-4 animation of the dust particle density evolution.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from math import floor
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
from pathlib import Path
from scipy import integrate

# Collect Athena++ inputs & outputs
athinput = athena_read.athinput('../athinput.si')
tlim = athinput['time']['tlim']            # max. simulation time
nx1 = athinput['mesh']['nx1']              # num. radial zones
nx2 = athinput['mesh']['nx2']              # num. vertical zones
nx3 = athinput['mesh']['nx3']              # num. azimuthal zones
c_s = athinput['hydro']['iso_sound_speed'] # sound speed
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
Np_tot = athinput['problem']['npx1']\
    *athinput['problem']['npx2']\
    *athinput['problem']['npx3']           # total number of particles
Np = Np_tot/nx1/nx2/nx3                    # theo avg num particles per cell
epsilon = athinput['problem']['epsilon']   # avg. dust/gas ρ-ratio in BG state
H = c_s/Omega                              # gas scale height
T = 2*np.pi/Omega                          # orbital period
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                        '.out2.*.athdf')))
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f']/H, data['x2f']/H
times, rhops, rhopmaxs = [], [], []        # times, dust densities & maximums

for output in outputs:                     # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time']/T)
    rhops.append(data['rhop'][0])          # [0] flattens 3D array
    rhopmaxs.append(np.amax(data['rhop']))

# Calc time-avg rhopmax during saturated state & set cmap log min/max
t_sat = 15                                 # determined graphically
i_sat = floor(len(times)*t_sat/tlim)
vmin = epsilon/Np                          # quantized minimum
vmax = integrate.simps(rhopmaxs[i_sat:], times[i_sat:])/(tlim-t_sat)

# Initialize first frame
clipped = np.clip(rhops[0], vmin, vmax)
fig, ax = plt.subplots(dpi=225)
ax.set_aspect('equal')
ax.set_title('$t={:.2f}$ / $T$'.format(times[0]))
ax.set_xlabel('$x$ / $H_g$')
ax.set_ylabel('$z$ / $H_g$')
img = ax.pcolormesh(xf, zf, clipped, norm=colors.LogNorm(vmin, vmax))
cb = plt.colorbar(img)
cb.set_label(r'$\rho_p$ / $\rho_0$')

def animate(i):
    """Update frame.

    Args:
        i: Frame number.
    """
    ax.set_title('$t={:.2f}$ / $T$'.format(times[i]))
    clipped = np.clip(rhops[i].ravel(), vmin, vmax) # flattens, clips array
    img.set_array(clipped)
    img.set_clim(vmin, vmax)

# Compile and save animation
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title='Dust Density', artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path'] = '/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=14500)
anim.save('../movies/rhop_full.mp4', writer=writer)
