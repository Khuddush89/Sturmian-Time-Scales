# 📐 Mathematical conventions

## Positive gauge

The proportional conformable derivative is

$$x^{\Delta_\alpha}=\kappa_1x+\kappa_0x^\Delta.$$

For graininess $\mu$, assume $\kappa_0>0$ and
$\nu=\kappa_0-\mu\kappa_1>0$. The positive gauge
$E_0=e_{-\kappa_1/\kappa_0}(\cdot,t_0)$ gives $x=E_0y$ and

$$\mathcal L_\alpha(E_0y)=\nu E_0\left[(p\nu y^\Delta)^\Delta
+\frac{q}{\kappa_0}y^\sigma\right].$$

Thus the reduced coefficients are $P=p\nu$ and $Q=q/\kappa_0$.
Positivity of the gauge preserves nodal zeros and signs across gaps.

## Finite grids

For nodes $t_0<\cdots<t_N$, the supplied grid is the time scale itself.
It is not an interpolation of a continuous domain. The second-order equation
is imposed with the final two nodes removed, on indices $0,\ldots,N-2$.

The Dirichlet unknowns are $y(t_1),\ldots,y(t_{N-1})$.
An auxiliary terminal flux extension supports integration by parts; it does
not impose another equation for the endpoint value of the solution.

The solver forms $Ay=\lambda By$ with a symmetric tridiagonal stiffness
matrix and positive diagonal mass matrix. It returns ascending eigenvalues;
reduced eigenvectors satisfy $Y^TBY=I$.

## Generalized zeros

- Count a nodal zero at its node.
- Count a strict sign crossing across a gap at the gap's left endpoint.
- Use only gaps wholly contained in the supplied domain or reporting horizon.
- Exclude the initial nodal zero from shooting counts.
- Use positive rescaling of the reduced state for long-horizon sign counts.

Exact equality detects nodal zeros in `count_zeros`; no hidden zero tolerance
is applied. Near-zero floating-point values need a separately justified
convention when extending the examples.

## Canonical gains

The examples commonly use $\kappa_0=\alpha$ and $\kappa_1=1-\alpha$,
with $0<\alpha\le1$. Every gap must satisfy
$\alpha>\mu/(1+\mu)$; equality is inadmissible.

For the pulse time scale $[0,1]\cup[1.5,2.5]\cup[3,4]$, gaps have
length $1/2$, so admissibility requires $\alpha>1/3$.
Dense blocks use exact transfer matrices; gaps use separate jump matrices.

## Figures and proof scope

Solid curves occur only on dense blocks. Dotted connectors across gaps are
visual guides. Curves against $\alpha$ are parameter curves, not trajectories
on a time scale.

The 535 assertions verify finite examples, identities, spectral computations,
and selected high-precision calculations. They do not prove general theorems
or determine infinite-horizon oscillation by finite numerical sampling.
