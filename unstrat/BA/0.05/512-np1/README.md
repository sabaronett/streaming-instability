# Problem 1A

Unstratified monodisperse streaming instability.
This is run BA of Johansen \& Youdin (2007, ApJ, 662, 627).

The box should have dimensions $L_x=L_z = 2H \times 2H$, resolution is $512\times 512$, Stokes number ${\rm St} = 1$, dust-to-gas ratio $\varepsilon=0.2$, number of particles: 1 particle per cell (262,144 particles).

Snapshots taken at 5, 10, 20, 50, 100 orbits (units of $2\pi/\Omega$).
Snapshots should contain particle density and particle positions.

*Objective: compare cumulative distribution function, maximum density, density field, dispersion.*


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
zone-cycles/cpu_second = 3.1177433222655527e+07
```
