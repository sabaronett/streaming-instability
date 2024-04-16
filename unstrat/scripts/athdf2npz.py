#!/usr/bin/env python3
#==============================================================================
# athdf2npz.py
#
# Reads Athena++ HDF5 outputs and saves them to compressed NumPy pickles (.npz)
# NOTE: This currently assumes 2D data (see `data[f'var'][0]` on Line )
#
# Author: Stanley A. Baronett
# Created: 2024-04-15
# Updated: 2024-04-15
#==============================================================================

import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
import os
from pathlib import Path

print('Loading data...', flush=True)
out_no = int(sys.argv[1])                    # <output#> in athinput.*
t_max = float(sys.argv[2])                   # [orbital period T]
athinput = athena_read.athinput('athinput.si')
problem_id = athinput['job']['problem_id']
dt = athinput[f'output{out_no}']['dt']       # [T]
i_max = int(t_max/dt) + 1                    # account for initial output
outputs = sorted(list(Path('athdf').glob(f'{problem_id}.out{out_no}.*.athdf')))
os.makedirs('npz', exist_ok=True)

for i, output in enumerate(outputs[:i_max]):
    print('  {:.2%}'.format(i/len(outputs[:i_max])), flush=True)
    data = athena_read.athdf(output)
    kwds = {'Time': data['Time']}

    for j, arg in enumerate(sys.argv[3:]):
        kwds.update({f'{str(arg)}': data[f'{str(arg)}']})
    
    np.savez_compressed(f'npz/{problem_id}.out.{i:05d}', **kwds)

print(f'Done.')
