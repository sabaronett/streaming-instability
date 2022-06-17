|Author | Stanley A. Baronett|
|-------|--------------------|
|Created| 2022-06-17         |
|Updated| 2022-06-17         |

# Run 'AB/0.05'
  - Originally run name from Johansen & Youdin (2007)
    - St = 0.1
    - Dust-to-gas density ratio, epsilon = 1.0
    - Radial pressure gradient, Pi = 0.02
    - 2D (axisymmetric), unstratified shearing box
    - L_x x L_z = 0.1 H x 0.1 H
  - Resolution: 2048 x 2048 cells
  - Code units:
    - Length:  H     (gas scale height)
    - Time:    T     (orbital period)
    - Density: rho_g (gas density)
  - Dust density evolution video: https://youtu.be/vumcOoYBP_k

## File Contents
The following files are complete grid snapshot of the gas at t = 10 T (orbits)
- 'rhog.txt': gas density
- 'vgx.txt': radial gas velocity
- 'vgy.txt': azimuthal gas velocity
- 'vgz.txt': vertical gas velocity
