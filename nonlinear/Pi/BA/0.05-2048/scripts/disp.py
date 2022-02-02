#!/usr/bin/env python3
#=======================================================================
# disp.py
#
# Computes the maximum Euclidean particle displacement between outputs.
#
# Author: Stanley A. Baronett
# Created: 2021-02-01
#=======================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path

athinput = athena_read.athinput('../athinput.si')
outputs = sorted(list(Path('../dat').glob\
                           (athinput['job']['problem_id']+'.pout.*.dat')))
# prepare first displacement calc
first = '../dat/SI.pout.00000.dat'
time, plast = athena_read.particles(first)
disp_max = 0.0

for output in outputs:
    time, pdata = athena_read.particles('%s'%output)
    disp = np.sqrt((pdata['xp']-plast['xp'])**2+(pdata['yp']-plast['yp'])**2\
                   +(pdata['zp']-plast['zp'])**2) # Euclidean distance
    plast = pdata                                 # prepare for next calc
    temp = np.amax(disp)

    if (temp > disp_max):
        disp_max = temp                           # update max displacement

print('Maximum particle displacement: {:0.2f} H'.format(disp_max))
