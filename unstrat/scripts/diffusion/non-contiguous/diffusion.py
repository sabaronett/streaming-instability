#!/usr/bin/env python3
#==============================================================================
# diffusion.py
#
# Compute and save the dust diffusion coefficient.
#
# Author: Stanley A. Baronett, Chao-Chin Yang
# Created: 2021-05-10
# Last Modified: 2023-01-10
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
from dust import Diffusion, MultiSpecies
import numpy as np

# Get the stopping times, solid-to-gas ratios, and no. of outputs to process
athinput = athena_read.athinput('athinput.si')
# tausmin, tausmax = 10**par.logtausmin, 10**par.logtausmax
# ms = MultiSpecies(par=par)
taus = athinput['problem']['omega']*athinput['particles']['taus0']
epsilon = athinput['problem']['epsilon']
nout = 3         # Minimum number needed to calculate the diffusion coefficient
if sys.argv[1]: nout = int(sys.argv[1])

# Find the displacement and the diffusion coefficient
diff = Diffusion(athinput=athinput, nout=nout)
dpx, dpy, dpz = diff.coefficient()
np.savez_compressed(f'output/dcoeff', taus=taus, dpx=dpx, dpy=dpy, dpz=dpz)
