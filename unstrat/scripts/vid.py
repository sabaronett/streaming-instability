#!/usr/bin/env python3
#==============================================================================
# vid.py
#
# Create an MPEG-4 animation of the dust density field.
#
# Author: Stanley A. Baronett
# Created: 2022-02-08
# Updated: 2022-09-04
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
import numpy as np
from pathlib import Path

# Collect Athena++ I/O, set plotting configurations
athinput = athena_read.athinput('athinput.si')
case, vmin, vmax = 'AB', 1e-1, 1e1
if athinput['problem']['epsilon'] == 0.2:
    case, vmin, vmax = 'BA', 2e-2, 2e0
res, dpi = athinput['mesh']['nx1'], 450    # 2160p default
if res < 2048: dpi = 225                   # 1080p for lower resolution runs
c_s = athinput['hydro']['iso_sound_speed']
Omega = athinput['problem']['omega']
Pi = athinput['problem']['duy0']
Hg = c_s/Omega
T = 2*np.pi/Omega
outputs = sorted(list(Path('athdf').glob(athinput['job']['problem_id']+
                                         '.out2.*.athdf')))
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f']/Hg, data['x2f']/Hg
times, rhops = [], []

for output in outputs:
    data = athena_read.athdf(output)
    times.append(data['Time']/T)
    rhops.append(data['rhop'][0])

# Initialize first frame
clip = np.clip(rhops[0], vmin, vmax)
fig, ax = plt.subplots(dpi=dpi)
ax.set(aspect='equal', title=f'{case}, $\Pi=${Pi}, $t={times[0]:.2f}T$',
       xlabel='$x/H_\mathrm{g}$', ylabel='$z/H_\mathrm{g}$')
img = ax.pcolormesh(xf, zf, clip, norm=colors.LogNorm(vmin, vmax), cmap='plasma')
cb = plt.colorbar(img)
cb.set_label(r'$\rho_\mathrm{p}/\rho_\mathrm{g,0}$')

def animate(i):
    """
    Update frame.

    Parameters
    ----------
        i : int
            Frame number.
    """
    ax.set_title(f'{case}, $\Pi=${Pi}, $t={times[i]:.2f}T$')
    clip = np.clip(rhops[i].ravel(), vmin, vmax)
    img.set_array(clip)
    img.set_clim(vmin, vmax)
    print(f'  frame {i:4n}', flush=True)

# Compile and save animation
print('Processing frames...', flush=True)
title = f'{case}-Pi{Pi:.2f}-{res}'
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title=(title+' Dust Density'), artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path']='/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save(title+'_rhop.mp4', writer=writer)
print('Done.\nVideo saved.', flush=True)
