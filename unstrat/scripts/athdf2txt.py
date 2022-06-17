#!/usr/bin/env python3
#==============================================================================
# athdf2txt.py
#
# Reads and saves gas (density and velocity) data to a text file (.txt).
#
# Author: Stanley A. Baronett
# Created: 2022-06-17
# Updated: 2022-06-17
#==============================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np

output = str(sys.argv[1])
print('Loading gas data...')
data = athena_read.athdf(output)
print(f'Done.\nSaving data to .txt...')
np.savetxt('rhog.txt', data['rho'][0])
np.savetxt('vgx.txt', data['vel1'][0])
np.savetxt('vgz.txt', data['vel2'][0])
np.savetxt('vgy.txt', data['vel3'][0])
print(f'Done.')
