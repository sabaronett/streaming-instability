"""Animate Dust Density.

Compiles and saves an MPEG-4 animation of the dust particle density.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from math import floor
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
import os
from pathlib import Path

# Set cmap log min/max
vmin = 0.04
vmax = 40

# Collect Athena++ inputs & outputs
athinput = athena_read.athinput('../athinput.si')
c_s = athinput['hydro']['iso_sound_speed'] # sound speed
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
H = c_s/Omega                              # gas scale height
T = 2*np.pi/Omega                          # orbital period
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"]+
                                            '.out1.*.athdf')))
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f']/H, data['x2f']/H
times, rhops = [], []                      # times, dust densities

for output in outputs:                     # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time']/T)
    rhops.append(data['rhop'][0])          # [0] flattens 3D array

# Initialize first frame
clipped = np.clip(rhops[0], vmin, vmax)
fig, ax = plt.subplots(dpi=120, figsize=(16,9))
img = ax.pcolormesh(xf, zf[8:2168], clipped[8:2168,:], cmap='afmhot',
                    norm=colors.LogNorm(vmin, vmax), shading='auto')
ax.set(xticks=[], yticks=[], frame_on=False)
time_text = ax.text(0.93, -0.5, '', color='white', fontsize=24, ha='right')
ax.arrow(-0.78, 0.5, -0.15, 0, width=0.01, length_includes_head=True,
         shape='right', fc='white', ec='None')
ax.text(-0.89, 0.45, 'star', color='white', size=24)
ax.arrow(-0.93, -0.5, 0.15, 0, width=0.01, length_includes_head=True,
         shape='right', fc='white', ec='None')
ax.text(-0.83, -0.47, r'$r$', color='white', size=24)
ax.arrow(-0.93, -0.5, 0, 0.15, width=0.01, length_includes_head=True,
         shape='left', fc='white', ec='None')
ax.text(-0.905, -0.39, r'$z$', color='white', size=24)
fig.tight_layout(pad=0)
fig.patch.set_facecolor('black')

def init():
    time_text.set_text('{:2.0f} years'.format(times[0]))
    return time_text

def animate(i):
    """Update frame.

    Args:
        i: Frame number.
    """
    time_text.set_text('{:2.0f} years'.format(times[i]))
    ax.arrow(-0.78, 0.5, -0.15, 0, width=0.01, length_includes_head=True,
         shape='right', fc='white', ec='None')
    ax.text(-0.89, 0.45, 'star', color='white', size=24)
    ax.arrow(-0.93, -0.5, 0.15, 0, width=0.01, length_includes_head=True,
            shape='right', fc='white', ec='None')
    ax.text(-0.83, -0.47, r'$r$', color='white', size=24)
    ax.arrow(-0.93, -0.5, 0, 0.15, width=0.01, length_includes_head=True,
            shape='left', fc='white', ec='None')
    ax.text(-0.905, -0.39, r'$z$', color='white', size=24)
    clipped = np.clip(rhops[i][8:2168,:].ravel(), vmin, vmax) # flattens, clips array
    img.set_array(clipped)
    img.set_clim(vmin, vmax)
    print('Frame {:3d}'.format(i))

# Compile and save animation
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(times),
                               repeat=False)
metadata = dict(title='Dust Density', artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path']='/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=30, metadata=metadata, bitrate=-1)
anim.save('../video/HECC_Annual_Report_SI.mp4', writer=writer, savefig_kwargs={'facecolor':'black'})