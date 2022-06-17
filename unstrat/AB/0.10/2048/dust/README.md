|Author | Stanley A. Baronett|
|-------|--------------------|
|Created| 2022-06-17         |
|Updated| 2022-06-17         |

# Run 'AB/0.05'
  - Originally run name from Johansen & Youdin (2007)
    - St = 0.1
    - Dust-to-gas density ratio, epsilon = 1.0
    - Radial pressure gradient, Pi = 0.1
    - 2D (axisymmetric), unstratified shearing box
    - L_x x L_z = 0.1 H x 0.1 H
  - Resolution: 2048 x 2048 cells
  - Code units:
    - Length:  H     (gas scale height)
    - Time:    T     (orbital period)
    - Density: rho_g (gas density)
  - Dust density evolution video: https://youtu.be/dSfsSFgFE0o

## File Contents
- 'SI.pout.00001.txt' (ASCII)
  - Lagrangian particle data snapshot of the dust at t = 10 T (orbits)
  - Columns:
    1. xp
    2. yp
    3. zp
    4. vpx
    5. vpy
    6. vpz
