"""Dust Density Autocorrelation Stack (ACS)

Compute autocorrelation of dust density and sum over multiple snapshots.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path
from scipy import fftpack

# Collect .athdf outputs, init sim consts.
athinput = athena_read.athinput('../athinput.si')
outputs = sorted(list(Path('../athdf').glob(athinput['job']['problem_id'] +
                                            '.out2.019[5-9]*.athdf')))
data = athena_read.athdf(outputs[0])
stack = np.zeros_like(data['rhop'][0])

for output in outputs:
    dataFT = fftpack.fft2(data['rhop'][0])
    dataAC = fftpack.ifft2(dataFT*np.conjugate(dataFT)).real
    shifted = fftpack.fftshift(dataAC)
    stack = np.add(stack, shifted)

np.savez_compressed('autocorrelation/acs-{:d}'.format(len(outputs)),acs=stack)
