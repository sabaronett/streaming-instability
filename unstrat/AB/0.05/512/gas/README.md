| Author  | Stanley A. Baronett |
|---------|---------------------|
| Created | 2022-07-07          |
| Updated | 2022-07-07          |

# Run 'AB/0.05'
  - Originally run name from Johansen & Youdin (2007)
  - 2D (axisymmetric), unstratified shearing box
  - Parameters:
    | $\tau_s$ | $\epsilon$ | $\Pi$ | $L_x \times L_z$ | Res.   |
    |----------|------------|-------|------------------|--------|
    | 0.1      | 1.0        | 0.05  | $0.1^2\,[H^2]$   | $512^2$|
  - Code units:
    | Dimension | Unit                   |
    |-----------|------------------------|
    | Length    | $H$ (gas scale height) |
    | Time      | $T$ (orbital period)   |
    | Density   | $\rho_g$ (gas density) |

## File Contents
The following files come from a complete grid snapshot of the gas at t = 10 T (orbits):
| Filename | Contents               |
|----------|------------------------|
| rhog.txt | gas density            |
| vgx.txt  | radial gas velocity    |
| vgy.txt  | azimuthal gas velocity |
| vgz.txt  | vertical gas velocity  |
