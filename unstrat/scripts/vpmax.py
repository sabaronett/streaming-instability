#!/usr/bin/env python3
#=======================================================================
# vpmax.py
#
# Prints the maximum radial and vertical particle velocities among
# particle output (.pout) data.
#
# Author: Stanley A. Baronett
# Created: 2022-06-09
# Updated: 2022-06-10
#=======================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np

cases = ['AB', 'BA']
Pis = ['0.01', '0.02', '0.05', '0.10']
res = 2048

print('Run\t\tv_{p,x,max} [H/T]\t\tv_{p,z,max} [H/T]')
print('=======================================================================')

for case in cases:
    for Pi in Pis:
        run = f'{case:s}/{Pi:s}'
        output = '../'+run+f'/{res:n}/dat/SI.pout.00001.dat'
        time, pdata = athena_read.particles(output)
        vpxmax, vpzmax = np.amax(pdata['vpx']),  np.amax(pdata['vpz'])

        print(f'{run:s}\t\t{vpxmax:.16e}\t\t{vpzmax:.16e}')
