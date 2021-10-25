"""Dust Density Autocorrelation Stack (ACS)

Compute autocorrelation of dust density and sum over multiple snapshots.
"""
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
import numpy as np
from pathlib import Path
from scipy import fft

# Collect .athdf outputs, init sim consts.
athinput = athena_read.athinput('../athinput.si')
outputs = sorted(list(Path('../athdf').glob(athinput['job']['problem_id'] +
                                            '.out2.019[5-9]*.athdf')))
data = athena_read.athdf(outputs[0])
stack = np.zeros_like(data['rhop'][0])

for output in outputs:
    dataFT = fft.fft2(data['rhop'][0], workers=16)
    dataAC = fft.ifft2(dataFT*np.conjugate(dataFT), workers=16).real
    shifted = fft.fftshift(dataAC)
    stack = np.add(stack, shifted)

np.savez_compressed('autocorrelation/acs-{:d}'.format(len(outputs)),acs=stack)
