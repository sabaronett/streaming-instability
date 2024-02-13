# Problem 1C

Same as Problem 1A but with a resolution of $1024^2$ (also 1 particle per cell, i.e., 1,048,576 particles).

Snapshots taken at the same times.

*Objective: compare behavior of drag force backreaction with resolution.*


## Directory

### `athdf/`

Contains snapshots in [Athena++'s own HDF5 format](https://github.com/PrincetonUniversity/athena/wiki/HDF5-Format).
The particle density field (as mapped to the gas grid via the standard particle-mesh method assignment) can be retrieved using [Athena++'s Python reader](https://github.com/PrincetonUniversity/athena/wiki/Reading-Data-into-Python) (also included in parent folder) and this sample script:
```python
import sys
sys.path.insert(0, '../athena/vis/python')
import athena_read

data = athena_read.athdf('../athdf/SI.out1.00005.athdf')
data['rhop'][0]
```


### `dat/`

Contains snapshots in binary format.
The positions of particles can be retrieved using the modified Python reader for Athena++ in the parent folder and this sample script:
```python
import athena_read

time, pdata = athena_read.particles('../dat/SI.pout.00005.dat')
pdata['xp']
pdata['yp']
pdata['zp']
```
where `'yp'` and `'zp'` keys correspond to vertical and azimuthal coordinates in the reference frame of the shearing box, respectively.


### `performance.txt`

The entire job output to console, containing performance metrics at the end, i.e.,
```
zone-cycles/cpu_second = 7.3153490690240011e+07
```
