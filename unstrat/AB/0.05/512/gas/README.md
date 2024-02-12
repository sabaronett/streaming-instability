|Author | Stanley A. Baronett|
|-------|--------------------|
|Created| 2024-02-12         |
|Updated| 2024-02-12         |

# 'AB/0.05/512'
  - Original parameters from Johansen & Youdin (2007)
    - $\tau_s = 0.1$
    - Dust-to-gas density ratio $\epsilon = 1.0$
    - Radial pressure gradient $\Pi = 0.01$
    - 2D axisymmetric unstratified shearing box
    - $L_x \times L_z = 0.1^2 H_g^2$
  - Resolution: 512 x 512 cells
  - Code units:
    - Length:  $H_g$    (gas scale height)
    - Time:    $T$      (orbital period)
    - Density: $\rho_g$ (gas density)

## File Contents
The following gas related snapshots are taken at $t = 10T$
- 'rhog.txt': gas density
- 'vgx.txt': radial gas velocity
- 'vgy.txt': azimuthal gas velocity
- 'vgz.txt': vertical gas velocity
