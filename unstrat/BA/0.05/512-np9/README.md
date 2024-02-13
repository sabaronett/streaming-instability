# Problem 1B

Same as problem 1A but with 9 particles per cell (2,359,296 particles).

Snapshots taken at the same times.

*Objective: compare convergence with number of particles per grid cell.*


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
zone-cycles/cpu_second = 7.0013032591338446e+06
```
