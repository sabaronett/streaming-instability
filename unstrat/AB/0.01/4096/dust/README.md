|Author | Stanley A. Baronett|
|-------|--------------------|
|Created| 2024-02-12         |
|Updated| 2024-02-13         |

# 'AB/0.01/4096'
  - Original parameters from Johansen & Youdin (2007)
    - $\tau_s = 0.1$
    - Dust-to-gas density ratio $\epsilon = 1.0$
    - Radial pressure gradient $\Pi = 0.01$
    - 2D axisymmetric unstratified shearing box
    - $L_x \times L_z = 0.1^2 H_g^2$
  - Resolution: 4096 x 4096 cells
  - Code units:
    - Length:  $H_g$    (gas scale height)
    - Time:    $T$      (orbital period)
    - Density: $\rho_g$ (gas density)

## File Contents
The following dust related snapshot is appended with the simulation time at which it was taken (e.g., $t = 20.0T$).
- 'dust_20.0T.txt' (ASCII)
  - Lagrangian particle data snapshot
  - Columns:
    1. xp
    2. zp
    3. yp
    4. vpx
    5. vpz
    6. vpy
