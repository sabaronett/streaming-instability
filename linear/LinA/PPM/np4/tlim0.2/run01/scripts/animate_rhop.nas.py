"""Compile and export animations of Athena++ simulation variables.

Leave one blank line.  The rest of this docstring should contain an
overall description of the module or program.  Optionally, it may also
contain a brief description of exported classes and functions and/or usage
examples.

  Typical usage example:

  foo = ClassFoo()
  bar = foo.FunctionBar()
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
                                        '.out1.*.athdf')))
c_s = athinput['hydro']['iso_sound_speed'] # sound speed
Omega = athinput['problem']['omega']       # local Keplerian angular frequency
H = c_s / Omega                            # gas scale height
T = 2*np.pi/Omega                          # orbital period
data = athena_read.athdf(outputs[0])
xf, zf = data['x1f'] / H, data['x2f'] / H
times = []                                 # sim output times
rhops = []                                 # particle density

for output in outputs:               # load all data into memory
    data = athena_read.athdf(output)
    times.append(data['Time'] / T)
    rhops.append(data['rhop'][0])    # [0] effectively flattens 3D array

# Initialize first frame
fig, ax = plt.subplots(dpi=200)
ax.set_aspect('equal')
ax.set_title('$t={:.3f}$'.format(times[0]))
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
writer = animation.FFMpegWriter(fps=30, metadata=metadata, bitrate=14500)
anim.save('../movies/rhop.nas.mp4', writer=writer)