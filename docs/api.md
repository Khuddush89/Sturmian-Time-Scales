# 🛠️ Model API

Import numerical functions from `sturmian.model`.
Inputs are NumPy-compatible one-dimensional node arrays and edge coefficients.
Coefficient arrays broadcast to `len(t) - 1` edges.

| Function | Returns / purpose |
| :--- | :--- |
| `harmonic(n)` | First `n` harmonic-sum nodes |
| `alternating(n)` | `n` nodes with alternating gaps `0.25`, `0.55` |
| `gauge(t, k0, k1)` | Edge `nu` and normalized nodal gauge `E` |
| `pencil(t, alpha, ...)` | Stiffness diagonal, off-diagonal, mass diagonal, gauge |
| `spectrum(t, alpha, vectors=False, ...)` | Ascending Dirichlet eigenvalues |
| `spectrum(..., vectors=True)` | Eigenvalues, reduced vectors `Y`, original vectors `X` |
| `propagate(t, alpha, ...)` | Reduced nodal states `Y`, `U`, crossing indices, nodal-zero indices |
| `count_zeros(t, alpha, q)` | Crossing and nodal-zero indices using positive rescaling |
| `dense_matrix(s, alpha, L=1)` | Exact `2 × 2` dense-block transfer |
| `gap_matrix(s, alpha, g=0.5)` | `2 × 2` jump transfer for one gap |
| `pulse_transfer(s, alpha)` | Transfer through the three-block pulse domain |
| `first_pulse_root(alpha)` | First positive pulse Dirichlet parameter |
| `pulse_mesh(m)` | Pulse domain with `m` subdivisions per dense block |
| `pulse_profile(alpha, ...)` | Dense-block nodes and reduced/original profiles |
| `hybrid_exact(alpha, ...)` | Jump nodes and DOP853 dense-block state solutions |
| `hybrid_mesh_error(alpha, m)` | Maximum state error on the reference hybrid domain |

`pencil` accepts scalar or edge values for `p`, `q`, `r`, `k0`, and `k1`.
It requires at least three nodes, finite coefficients, positive `p` and `r`,
positive `k0`, and a positive gauge condition.

```python
import numpy as np
from sturmian.model import pencil, spectrum

t = np.array([0.0, 0.25, 0.8, 1.0])
lam, y, x = spectrum(t, alpha=0.75, q=0.4, vectors=True)
_, _, mass, _ = pencil(t, alpha=0.75, q=0.4)
assert np.allclose(y.T @ (mass[:, None] * y), np.eye(len(lam)))
```

Use `rescale=True` in `propagate` only for signs. Rescaled states no longer
represent physical amplitudes. Gauge arrays can underflow on long domains;
`count_zeros` avoids using the gauge to decide signs.

The verification and plotting engines run sequentially and use per-run output
state. For simultaneous runs, use separate processes and output directories.
