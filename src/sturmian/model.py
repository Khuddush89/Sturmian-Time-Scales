"""Reproducible finite-grid and hybrid computations for proportional conformable time-scale equations.

All finite differences are taken on the supplied time scale. A flux value at
the terminal node is an auxiliary extension, not an additional equation for y.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq


def harmonic(n):
    return np.cumsum(1.0 / np.arange(1, n + 1))


def alternating(n):
    return np.r_[0.0, np.cumsum(np.resize([0.25, 0.55], n - 1))]


def _steps(t, minimum=2):
    """Require a finite, strictly increasing one-dimensional time scale."""
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or len(t) < minimum or not np.all(np.isfinite(t)):
        raise ValueError(f"The time scale must contain at least {minimum} finite nodes")
    h = np.diff(t)
    if np.any(h <= 0):
        raise ValueError("Time-scale nodes must be strictly increasing")
    return h


def gauge(t, k0, k1):
    h = _steps(t)
    k0 = np.broadcast_to(k0, h.shape)
    k1 = np.broadcast_to(k1, h.shape)
    if not np.all(np.isfinite(k0)) or np.any(k0 <= 0) or not np.all(np.isfinite(k1)):
        raise ValueError("Gains must be finite and k0 must be positive")
    nu = k0 - h * k1
    if np.any(nu <= 0):
        raise ValueError("The positive-gauge condition is not satisfied")
    return nu, np.r_[1.0, np.cumprod(nu / k0)]


def pencil(t, alpha, p=1.0, q=0.0, r=1.0, k0=None, k1=None):
    h = _steps(t, minimum=3)
    k0 = np.broadcast_to(alpha if k0 is None else k0, h.shape)
    k1 = np.broadcast_to(1 - alpha if k1 is None else k1, h.shape)
    p, q, r = [np.broadcast_to(z, h.shape) for z in (p, q, r)]
    if any(not np.all(np.isfinite(z)) for z in (p, q, r)) or np.any(p <= 0) or np.any(r <= 0):
        raise ValueError("Coefficients must be finite, with p and r positive")
    nu, E = gauge(t, k0, k1)
    c = p * nu / h
    diagonal = c[:-1] + c[1:] - h[:-1] * q[:-1] / k0[:-1]
    offdiagonal = -c[1:-1]
    mass = h[:-1] * r[:-1] / k0[:-1]
    return diagonal, offdiagonal, mass, E


def spectrum(t, alpha, p=1.0, q=0.0, r=1.0, vectors=False, k0=None, k1=None):
    d, e, b, E = pencil(t, alpha, p, q, r, k0, k1)
    result = eigh_tridiagonal(d / b, e / np.sqrt(b[:-1] * b[1:]), eigvals_only=not vectors)
    if not vectors:
        return result
    lam, V = result
    Y = V / np.sqrt(b[:, None])
    X = E[1:-1, None] * Y
    return lam, Y, X


def propagate(t, alpha, p=1.0, q=0.0, initial=(0.0, 1.0), rescale=False):
    """Reduced state (y,u); rescale=True is for sign counts, not amplitudes."""
    h = _steps(t)
    p, q = [np.broadcast_to(z, h.shape) for z in (p, q)]
    nu = alpha - h * (1 - alpha)
    if not 0 < alpha <= 1 or np.any(nu <= 0):
        raise ValueError("Nonadmissible order")
    y, u = initial
    Y, U, crossings, actual = [y], [u], [], []
    for j, hh in enumerate(h):
        yy = y + hh * u / (p[j] * nu[j])
        uu = u - hh * q[j] * yy / alpha
        if yy == 0:
            actual.append(j + 1)
        elif (y > 0 and yy < 0) or (y < 0 and yy > 0):
            crossings.append(j)
        if rescale:
            s = max(abs(yy), abs(uu))
            if s:
                yy, uu = yy / s, uu / s
        y, u = yy, uu
        Y.append(y)
        U.append(u)
    return np.asarray(Y), np.asarray(U), crossings, actual


def count_zeros(t, alpha, q):
    """Streaming positive rescaling: no amplitude or gauge underflow.

    A nodal zero is counted at its node. A crossing is counted at its left
    endpoint, and only edges wholly contained in the supplied grid are used.
    The initial nodal zero is excluded.
    """
    h = _steps(t)
    if not 0 < alpha <= 1:
        raise ValueError("Nonadmissible order")
    q = np.broadcast_to(q, h.shape)
    y, u = 0.0, 1.0
    crossings, actual = [], []
    for j, hh in enumerate(h):
        nu = alpha - hh * (1 - alpha)
        if nu <= 0:
            raise ValueError("Nonadmissible order")
        yy = y + hh * u / nu
        uu = u - hh * q[j] * yy / alpha
        if yy == 0:
            actual.append(j + 1)
        elif (y > 0 and yy < 0) or (y < 0 and yy > 0):
            crossings.append(j)
        s = max(abs(yy), abs(uu))
        y, u = yy / s, uu / s
    return crossings, actual


def dense_matrix(s, alpha, L=1.0):
    """Real block transfer, including negative spectral parameter and s=0."""
    if s == 0:
        return np.array([[1.0, L / alpha], [0.0, 1.0]])
    root = np.sqrt(abs(s))
    z = root * L / alpha
    if s > 0:
        return np.array([[np.cos(z), np.sin(z) / root], [-root * np.sin(z), np.cos(z)]])
    return np.array([[np.cosh(z), np.sinh(z) / root], [root * np.sinh(z), np.cosh(z)]])


def gap_matrix(s, alpha, g=0.5):
    nu = alpha - g * (1 - alpha)
    if not 0 < alpha <= 1 or g < 0 or nu <= 0:
        raise ValueError("The positive-gauge condition is not satisfied across the gap")
    return np.array([[1.0, g / nu], [-g * s / alpha, 1.0 - g * g * s / (alpha * nu)]])


def pulse_transfer(s, alpha):
    D, J = dense_matrix(s, alpha), gap_matrix(s, alpha)
    return D @ J @ D @ J @ D


def first_pulse_root(alpha):
    """Smallest positive root in s=lambda+q, verified by nodal tracking."""
    return brentq(
        lambda s: pulse_transfer(s, alpha)[0, 1], 0.0, np.pi**2 * alpha**2 / 9.0, xtol=5e-15
    )


def pulse_mesh(m):
    return np.r_[np.linspace(0, 1, m + 1), np.linspace(1.5, 2.5, m + 1), np.linspace(3, 4, m + 1)]


def pulse_profile(alpha, s=None, samples=401):
    if s is None:
        s = first_pulse_root(alpha)
    state = np.array([0.0, 1.0])
    E = 1.0
    blocks = []
    for left in [0.0, 1.5, 3.0]:
        tt = np.linspace(left, left + 1, samples)
        yy = np.array([(dense_matrix(s, alpha, z - left) @ state)[0] for z in tt])
        xx = E * np.exp(-(1 - alpha) * (tt - left) / alpha) * yy
        blocks.append((tt, yy, xx))
        state = dense_matrix(s, alpha) @ state
        E *= np.exp(-(1 - alpha) / alpha)
        if left < 3:
            state = gap_matrix(s, alpha) @ state
            E *= (alpha - 0.5 * (1 - alpha)) / alpha
    return blocks


def hybrid_exact(alpha, T=14.0, samples=301):
    """Original conformable state, x(0)=0, x^Delta(0)=0.1."""
    x, v = 0.0, 0.1 * alpha
    nodes = {0.0: np.array([x, v])}
    blocks = []

    def rhs(t, S):
        x, v = S
        p, q = 1 / (1 + t), 1 / (1 + t * t)
        return [(v / p - (1 - alpha) * x) / alpha, (-q * x - (1 - alpha) * v) / alpha]

    def jump(t, x, v):
        p, q = 1 / (1 + t), 1 / (1 + t * t)
        xx = x + (v / p - (1 - alpha) * x) / alpha
        vv = v + (-q * xx - (1 - alpha) * v) / alpha
        return xx, vv

    for k in range(int(T // 3) + 1):
        left = float(3 * k)
        if left >= T:
            break
        x, v = jump(left, x, v)
        nodes[left + 1] = np.array([x, v])
        stop = min(left + 2, T)
        sol = solve_ivp(
            rhs,
            (left + 1, stop),
            [x, v],
            method="DOP853",
            rtol=2e-13,
            atol=2e-14,
            dense_output=True,
        )
        tt = np.linspace(left + 1, stop, samples)
        blocks.append((tt, sol.sol(tt), sol.sol))
        x, v = sol.y[:, -1]
        nodes[stop] = np.array([x, v])
        if stop >= T:
            break
        x, v = jump(left + 2, x, v)
        nodes[left + 3] = np.array([x, v])
    return nodes, blocks


def hybrid_mesh_error(alpha, m, T=14.0):
    nodes, blocks = hybrid_exact(alpha, T)
    grid = np.sort(
        np.unique(
            np.r_[
                np.arange(0, T + 1, 3),
                *[np.linspace(3 * k + 1, 3 * k + 2, m + 1) for k in range(5)],
            ]
        )
    )
    x, v = 0.0, 0.1 * alpha
    error = 0.0
    for j, h in enumerate(np.diff(grid)):
        t = grid[j]
        p, q = 1 / (1 + t), 1 / (1 + t * t)
        x += h * (v / p - (1 - alpha) * x) / alpha
        v += h * (-q * x - (1 - alpha) * v) / alpha
        tt = float(grid[j + 1])
        reference = nodes.get(tt)
        if reference is None:
            for ts, _, fun in blocks:
                if ts[0] <= tt <= ts[-1]:
                    reference = fun(tt)
                    break
        if reference is None:
            raise ValueError("Mesh node is not on the reference time scale")
        error = max(error, float(np.max(np.abs([x, v] - reference))))
    return error
