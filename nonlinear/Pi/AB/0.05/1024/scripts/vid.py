#!/usr/bin/env python3
#==============================================================================
# vid.py
#
# Compiles and saves an MPEG-4 animation of the dust particle density.
#
# Author: Stanley A. Baronett
# Created: 2022-02-08
# Last Modified: 2022-02-08
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
from pathlib import Path

# Set cmap log min/max, dpi
vmin = 1e-1 # CPDD > 70%
vmax = 1e1  # CPDD ~ 1%
dpi = 225   # 1080p

# Collect Athena++ inputs & outputs
athinput = athena_read.athinput('../athinput.si')
c_s = athinput['hydro']['iso_sound_speed'] # sound speed
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
H = c_s/Omega                              # gas scale height
T = 2*np.pi/Omega                          # orbital period
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                            '.out2.*.athdf')))
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f']/H, data['x2f']/H
times, rhops = [], []                      # times, dust densities

for output in outputs:                     # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time']/T)
    rhops.append(data['rhop'][0])          # [0] flattens 3D array

# Initialize first frame
clipped = np.clip(rhops[0], vmin, vmax)
fig, ax = plt.subplots(dpi=dpi)
ax.set_aspect('equal')
ax.set_title('$t={:.2f}$ / $T$'.format(times[0]))
ax.set_xlabel('$x$ / $H$')
ax.set_ylabel('$z$ / $H$')
img = ax.pcolormesh(xf, zf, clipped, norm=colors.LogNorm(vmin, vmax))
cb = plt.colorbar(img)
cb.set_label(r'$\rho_\mathrm{p}$ / $\rho_\mathrm{g,0}$')

def animate(i):
    """Update frame.

    Args:
        i: Frame number.
    """
    ax.set_title('$t={:.2f}$ / $T$'.format(times[i]))
    clipped = np.clip(rhops[i].ravel(), vmin, vmax) # flattens, clips array
    img.set_array(clipped)
    img.set_clim(vmin, vmax)
    print('Frame {:5d}'.format(i))

# Compile and save animation
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title='Dust Density', artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path'] = '/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save('../video/rhop.mp4', writer=writer)
