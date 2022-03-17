#!/usr/bin/env python3
#==============================================================================
# vid.py
#
# Creates an MPEG-4 animation of the dust particle density.
#
# Author: Stanley A. Baronett
# Created: 2022-02-08
# Last Modified: 2022-03-17
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
from pathlib import Path

# Collect Athena++ I/O, set plotting configurations
athinput = athena_read.athinput('athinput.si')
run, vmin, vmax = 'AB', 1e-2, 1e1          # CPDD > 70%; CPDD < 1%
if athinput['problem']['epsilon'] == 0.2:
    run, vmin, vmax = 'BA', 1e-2, 1e0      # CPDD*ε > 85%; CPDD < 5%
res, dpi = athinput['mesh']['nx1'], 450    # 2160p default
if res < 2048: dpi = 225                   # 1080p for lower resolution runs
c_s = athinput['hydro']['iso_sound_speed'] # sound speed
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
Pi = athinput['problem']['duy0']           # radial pressure gradient
H = c_s/Omega                              # gas scale height
T = 2*np.pi/Omega                          # orbital period
outputs = sorted(list(Path('athdf').glob(athinput["job"]["problem_id"]+
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
ax.set(aspect='equal', title=f'$t={times[0]:.2f}$ / $T$',
       xlabel='$x$ / $H$', ylabel='$z$ / $H$')
img = ax.pcolormesh(xf, zf, clipped, norm=colors.LogNorm(vmin, vmax))
cb = plt.colorbar(img)
cb.set_label(r'$\rho_\mathrm{p}$ / $\rho_\mathrm{g,0}$')

def animate(i):
    """Update frame.

    Args:
        i: Frame number.
    """
    ax.set_title(f'{run:s}, $\Pi=${Pi:.2f}, $t={times[i]:.2}$ / $T$')
    clipped = np.clip(rhops[i].ravel(), vmin, vmax) # flatten, clip array
    img.set_array(clipped)
    img.set_clim(vmin, vmax)
    print(f'  frame {i:4n}', flush=True)            # print frame progess

# Compile and save animation
print('Processing frames...', flush=True)
Pi = athinput['problem']['duy0']
title = '%s-Pi%s-%s Dust Density'%(run, Pi, res)
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title=title, artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path']='/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save('%s-Pi%s-%s_rhop.mp4'%(run, Pi, res), writer=writer)
print('Done.\nVideo saved.', flush=True)
