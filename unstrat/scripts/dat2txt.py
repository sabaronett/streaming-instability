#!/usr/bin/env python3
#==============================================================================
# dat2txt.py
#
# Reads and saves particle (position and velocity) data to a text file (.txt).
#
# Author: Stanley A. Baronett
# Created: 2022-06-02
# Updated: 2022-06-17
#==============================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np

output = str(sys.argv[1])
print('Loading particle data...')
time, pdata = athena_read.particles(output)
fname = output[0:-4]+'.txt'
print(f'Done.\nSaving to {fname:s}...')
np.savetxt(fname, pdata)
print(f'Done.')
