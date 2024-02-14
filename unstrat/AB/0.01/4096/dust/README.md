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
The following dust related snapshots are appended with the simulation time at which they were taken (e.g., $t = 20.0T$)
- 'rhop_20.0T.txt': dust density
- 'vpx_20.0T.txt': radial dust velocity
- 'vpy_20.0T.txt': azimuthal dust velocity
- 'vpz_20.0T.txt': vertical dust velocity
