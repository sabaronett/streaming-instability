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
res, dpi = athinput['mesh']['nx1'], 450    # 2160p default
if res < 2048: dpi = 225                   # 1080p for lower resolution cases
fig, ax = plt.subplots(dpi=dpi)
case = 'AB'
if athinput['problem']['epsilon'] == 0.2: case = 'BA'
Pi = athinput['problem']['duy0']           # radial pressure gradient
vmin, vmax = 1e-1, 1e1                     # default for AB
if case == 'BA': vmin, vmax = 2e-2, 2e0
c_s = athinput['hydro']['iso_sound_speed']
Omega = athinput['problem']['omega']
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
clips = np.clip(rhops, vmin, vmax)
frame = clips[0]
ax.set(aspect='equal', title=f'{case:s}, $\Pi=${Pi}, $t={times[0]:.2f}T$',
       xlabel=r'$x$ / $H_\mathrm{g}$', ylabel=r'$z$ / $H_\mathrm{g}$')
img = ax.pcolormesh(xf, zf, frame, norm=colors.LogNorm(), cmap='plasma')
cb = plt.colorbar(img)
cb.set_label(r'$\rho_\mathrm{p}$ / $\rho_\mathrm{g,0}$')

def animate(i):
    """
    Update frame.

    Parameters
    ----------
        i : int
            Frame number.
    """
    ax.set_title(f'{case:s}, $\Pi=${Pi}, $t={times[i]:.2f}T$')
    frame = clips[i]
    img.set_array(frame)
    print(f'  frame {i:4n}', flush=True)   # print frame progess

# Compile and save animation
print('Processing frames...', flush=True)
title = f'{case}-Pi{Pi:.2f}-{res}'
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title=(title+' Dust Density'), artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path']='/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save(title+'_rhop.mp4', writer=writer)
print('Done.\nVideo saved.', flush=True)
