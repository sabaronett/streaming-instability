#!/usr/bin/env python3
#==============================================================================
# athdf2txt.py
#
# Reads and saves gas (density and velocity) data to a text file (.txt).
#
# Author: Stanley A. Baronett
# Created: 2022-06-17
# Updated: 2024-02-13
#==============================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np

output = str(sys.argv[1])
print('Loading gas data...')
data = athena_read.athdf(output)
time = data['Time']
print(f'Done.\nSaving data to .txt...')
np.savetxt(f'rhog_{time:.1f}T.txt', data['rho'][0])
np.savetxt(f'vgx_{time:.1f}T.txt', data['vel1'][0])
np.savetxt(f'vgz_{time:.1f}T.txt', data['vel2'][0])
np.savetxt(f'vgy_{time:.1f}T.txt', data['vel3'][0])
print(f'Done.')
