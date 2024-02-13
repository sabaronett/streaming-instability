#!/usr/bin/env python3
#==============================================================================
# athdf_dust2txt.py
#
# Reads and saves dust (density and velocity) data to a text file (.txt).
#
# Author: Stanley A. Baronett
# Updated: 2024-02-13
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
np.savetxt(f'rhop_{time:.1f}T.txt', data['rhop'][0])
np.savetxt(f'vpx_{time:.1f}T.txt', data['vp1'][0])
np.savetxt(f'vpz_{time:.1f}T.txt', data['vp2'][0])
np.savetxt(f'vpy_{time:.1f}T.txt', data['vp3'][0])
print(f'Done.')
