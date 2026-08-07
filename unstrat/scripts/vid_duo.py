#!/usr/bin/env python3
#===============================================================================
# vid.py
#
# Create an MPEG-4 animation of the gas and dust density fields
#
# Author: Stanley A. Baronett
# Created: 2026-08-07
# Updated: 2026-08-07
#===============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
from matplotlib import animation
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

# Load data and compute normalization factors
athinput = athena_read.athinput('athinput.si')
case, vmin, vmax = 'AB', 0.1, 10
if athinput['problem']['epsilon'] == 0.2:
    case, vmin, vmax = 'BA', 0.02, 2
if athinput['job']['problem_id'] == 'As':
    case, vmin, vmax = 'As', 0.01, 1
    n_px = int(athinput['problem']['npx1']/athinput['mesh']['nx1'])
    n_pz = int(athinput['problem']['npx2']/athinput['mesh']['nx2'])
    n_p = n_px*n_pz
    drhog = 1e-3
res, dpi = athinput['mesh']['nx1'], 450    # 2160p (4K) default
if res < 2048: dpi = 225                   # 1080p for lower resolution runs
# elif res == 4096: dpi = 900                # 4320p (8K)
c_s = athinput['hydro']['iso_sound_speed']
Omega = athinput['problem']['omega']
Pi = athinput['problem']['duy0']
Hg = c_s/Omega
T = 2*np.pi/Omega
outputs = sorted(list(Path('athdf').glob(athinput['job']['problem_id']+
                                         '.out1.*.athdf')))
data = athena_read.athdf(outputs[0])
xv, zv = data['x1v']/Hg, data['x2v']/Hg
times, rhogs, rhops = [], [], []

for output in outputs:
    data = athena_read.athdf(output)
    times.append(data['Time']/T)
    rhogs.append(data['rho'][0])
    rhops.append(data['rhop'][0])

# Initialize first frame
fig, axs = plt.subplots(1, 2, sharey=True, dpi=dpi)

# Gas density field
img_rhog = axs[0].pcolormesh(xv, zv, rhogs[0], cmap='RdBu_r', vmin=1-drhog,
                              vmax=1+drhog)
mpl.rcParams["axes.formatter.offset_threshold"] = 2
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-1, 1))
cb_rhog = fig.colorbar(img_rhog, ax=axs[0], location='top', format=formatter)
cb_rhog.set_label(r'$\rho_\mathrm{g}/\rho_\mathrm{g,0}$')

# Dust density field
img_rhop = axs[1].pcolormesh(xv, zv, np.clip(rhops[0], vmin, vmax),
                              norm=colors.LogNorm(), cmap='cividis')
cb_rhop = fig.colorbar(img_rhop, ax=axs[1], location='top')
cb_rhop.set_label(r'$\rho_\mathrm{p}/\rho_\mathrm{g,0}$')

# Format plot
fig.suptitle(rf'Model {case} ({res}$^2$, $n_\mathrm{{p}}=${n_p}),'\
             +f'$t/T={times[0]:.1f}$')
axs[0].set_xticks([-0.1, -0.05, 0, 0.05, 0.1])
axs[0].set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
axs[0].set(ylabel=r'$z/H_\mathrm{g}$')
for ax in axs:
    ax.set(aspect='equal', xlabel=r'$x/H_\mathrm{g}$')

def animate(i):
    """
    Update frame.

    Parameters
    ----------
        i : int
            Frame number.
    """
    fig.suptitle(rf'Model {case} ({res}$^2$, $n_\mathrm{{p}}=${n_p}),'\
                 +rf'$t/T={times[i]:.1f}$')
    img_rhog.set_array(rhogs[i].ravel())
    img_rhop.set_array(rhops[i].ravel())
    print(f'  frame {i:4n}', flush=True)

# Compile and save animation
print('Processing frames...', flush=True)
if case == 'As':
    title = f'{case}-{res}-np{n_p}'
else:
    title = f'{case}-Pi{Pi:.2f}-{res}'
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title=(title+' Dust Density'), artist='Stanley A. Baronett')
plt.rcParams['animation.ffmpeg_path']='/nasa/pkgsrc/sles12/2018Q3/bin/ffmpeg3'
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save(title+'_rhog_rhop.mp4', writer=writer)
print('Done.\nVideo saved.', flush=True)
