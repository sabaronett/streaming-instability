"""Animate Dust Density.

Compiles and saves an MPEG-4 animation of the dust particle density.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# Collect .athdf outputs, init sim consts. and grid
athinput = athena_read.athinput('../athinput.si')
outputs = sorted(list(Path('../athdf').glob(athinput["job"]["problem_id"] +
                                        '.out2.*.athdf')))
c_s = athinput['hydro']['iso_sound_speed'] # sound speed
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
H = c_s/Omega                              # gas scale height
T = 2*np.pi/Omega                          # orbital period
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f']/H, data['x2f']/H
times = []                                 # sim output times
rhops = []                                 # particle density

for output in outputs:                     # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time']/T)
    rhops.append(data['rhop'][0])          # [0] effectively flattens 3D array

# Initialize first frame
fig, ax = plt.subplots(dpi=225)
ax.set_aspect('equal')
ax.set_title('$t = {:.3f}$ / $T$'.format(times[0]))
ax.set_xlabel('$x$ / $H_g$')
ax.set_ylabel('$z$ / $H_g$')
img = ax.pcolormesh(xf, zf, rhops[0])
cb = plt.colorbar(img)
cb.set_label(r'$\rho_p$ / $\rho_0$')

def animate(i):
    """Update frame.

    Args:
        i: Frame number.
    """
    ax.set_title('$t={:.3f}$'.format(times[i]))
    img.set_array(rhops[i].ravel()) # flatten 2D array to 1D array
    img.set_clim(rhops[i].min(), rhops[i].max())

# Compile and save animation
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title='Dust Density', artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path'] = '/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=14500)
anim.save('../movies/rhop_full.mp4', writer=writer)