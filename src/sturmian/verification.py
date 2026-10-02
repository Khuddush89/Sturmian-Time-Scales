"""Independent audit of identities and numerical statements in the manuscript.

Run: sturmian verify --output outputs/reference
Results: data/*.csv and reports/verification_results.json.
Numerical checks illustrate the identities; they are not proofs of theorems.
"""

from __future__ import annotations

import csv
import importlib.metadata
import json
import platform
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .model import (
    alternating,
    count_zeros,
    first_pulse_root,
    gauge,
    harmonic,
    hybrid_mesh_error,
    pencil,
    propagate,
    pulse_mesh,
    pulse_profile,
    pulse_transfer,
    spectrum,
)

ROOT = Path("outputs/reference")
CHECKS = []


def check(name, condition, residual=None, tolerance=None):
    record = dict(name=name, passed=bool(condition))
    if residual is not None:
        record["residual"] = float(residual)
    if tolerance is not None:
        record["tolerance"] = float(tolerance)
    CHECKS.append(record)
    if not condition:
        raise AssertionError(record)


def close(name, a, b, tol=2e-10):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    error = np.max(np.abs(a - b)) / max(1.0, float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    check(name, error <= tol, error, tol)


def write_csv(name, fields, rows):
    path = ROOT / "data" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def invert_fraction(A):
    n = len(A)
    M = [list(row) + [F(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for j in range(n):
        pivot = next(k for k in range(j, n) if M[k][j])
        M[j], M[pivot] = M[pivot], M[j]
        v = M[j][j]
        M[j] = [x / v for x in M[j]]
        for k in range(n):
            if k != j:
                v = M[k][j]
                M[k] = [x - v * y for x, y in zip(M[k], M[j])]
    return [row[n:] for row in M]


def decimal_pulse_threshold():
    """80-digit bisection, independent of the floating-point root solver."""
    with localcontext() as context:
        context.prec = 80
        one, s = Decimal(1), Decimal(".45")
        root = s.sqrt()

        def trig(z):
            ss, cc = z, one
            ts, tc = z, one
            k = 1
            while max(abs(ts), abs(tc)) > Decimal("1e-82"):
                ts *= -z * z / Decimal((2 * k) * (2 * k + 1))
                tc *= -z * z / Decimal((2 * k - 1) * (2 * k))
                ss += ts
                cc += tc
                k += 1
            return ss, cc

        def mm(A, B):
            return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

        def terminal(a):
            ss, cc = trig(root / a)
            D = [[cc, ss / root], [-root * ss, cc]]
            g = Decimal(".5")
            nu = a - g * (one - a)
            J = [[one, g / nu], [-g * s / a, one - g * g * s / (a * nu)]]
            return mm(mm(mm(mm(D, J), D), J), D)[0][1]

        lo, hi = Decimal(".8"), Decimal(".95")
        for _ in range(245):
            mid = (lo + hi) / 2
            if terminal(mid) < 0:
                lo = mid
            else:
                hi = mid
        value = (lo + hi) / 2
        return str(value)


def exact_examples():
    a = F(3, 4)
    hs = [F(1, k) for k in range(2, 8)]
    y, u, E = F(0), F(1), F(1)
    x, v = F(0), F(1)
    zz, vz = F(1), F(0)
    rows = []
    expected_y = [
        F(0),
        F(4, 5),
        F(7, 30),
        -F(89, 270),
        -F(733, 1134),
        -F(43087, 57834),
        -F(1812397, 2602530),
    ]
    expected_E = [F(1), F(5, 6), F(20, 27), F(55, 81), F(154, 243), F(1309, 2187), F(3740, 6561)]
    for j in range(7):
        check(f"exact/gauge-solution-{j}", x == E * y and v == E * u)
        check(f"exact/harmonic-y-{j}", y == expected_y[j])
        check(f"exact/harmonic-gauge-{j}", E == expected_E[j])
        check(f"exact/Abel-{j}", x * vz - zz * v == -E * E)
        w = u / y if y else None
        rows.append(
            dict(m=j + 1, E0=str(E), y=str(y), x=str(x), w=str(w) if w is not None else "undefined")
        )
        if j == 6:
            break
        h = hs[j]
        nu = a - h * (1 - a)
        yy = y + h * u / nu
        uu = u - h * F(16, 3) * yy
        xx = x + h * (v - (1 - a) * x) / a
        vv = v + h * (-4 * xx - (1 - a) * v) / a
        znext = zz + h * (vz - (1 - a) * zz) / a
        vnext = vz + h * (-4 * znext - (1 - a) * vz) / a
        y, u, E, x, v, zz, vz = yy, uu, E * nu / a, xx, vv, znext, vnext
    write_csv("harmonic_exact.csv", list(rows[0]), rows)
    c = [F(5, 4), F(2), F(11, 4), F(7, 2)]
    A = [[c[0] + c[1], -c[1], F(0)], [-c[1], c[1] + c[2], -c[2]], [F(0), -c[2], c[2] + c[3]]]
    Ai = invert_fraction(A)
    determinant = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] ** 2 * A[2][2]
    check("exact/Green-determinant", determinant == F(1501, 32))
    EE = expected_E[:5]
    hh = hs[:4]
    nu = [a - h * (1 - a) for h in hh]
    G = [[-Ai[j][k] * EE[j + 1] / (nu[k] * EE[k]) for k in range(3)] for j in range(3)]
    expected = [
        [-F(944, 1501), -F(600, 1501), -F(288, 1501)],
        [-F(12800, 27 * 1501), -F(2600, 3 * 1501), -F(416, 1501)],
        [-F(15488, 81 * 1501), -F(3146, 9 * 1501), -F(488, 1501)],
    ]
    check("exact/Green-matrix", G == expected)
    xx = [F(0)] + [sum(G[j][k] * hh[k] for k in range(3)) for j in range(3)] + [F(0)]
    check(
        "exact/Green-forcing-solution",
        xx[1:-1] == [-F(744, 1501), -F(17008, 40527), -F(27064, 121581)],
    )
    xd = [(1 - a) * xx[j] + a * (xx[j + 1] - xx[j]) / hh[j] for j in range(4)]
    check(
        "exact/Green-original-equation",
        all((1 - a) * xd[j] + a * (xd[j + 1] - xd[j]) / hh[j] == 1 for j in range(3)),
    )

    # Full bivariate characteristic polynomial, with exact rational coefficients.
    def add(*polys):
        out = {}
        for p in polys:
            for k, v in p.items():
                out[k] = out.get(k, F(0)) + v
        return {k: v for k, v in out.items() if v}

    def mul(p, q):
        out = {}
        for (ia, il), v in p.items():
            for (ja, jl), w in q.items():
                k = (ia + ja, il + jl)
                out[k] = out.get(k, F(0)) + v * w
        return out

    def scale(p, c):
        return {k: c * v for k, v in p.items()}

    d = [
        {(0, 1): F(1, 2), (2, 0): F(-7), (1, 0): F(2)},
        {(0, 1): F(1, 3), (2, 0): F(-9), (1, 0): F(2)},
        {(0, 1): F(1, 4), (2, 0): F(-11), (1, 0): F(2)},
    ]
    off = [{(1, 0): F(-1), (2, 0): F(4)}, {(1, 0): F(-1), (2, 0): F(5)}]
    det = scale(
        add(
            mul(mul(d[0], d[1]), d[2]),
            scale(mul(d[0], mul(off[1], off[1])), -1),
            scale(mul(d[2], mul(off[0], off[0])), -1),
        ),
        24,
    )
    target = {
        (0, 3): F(1),
        (2, 2): F(-85),
        (1, 2): F(18),
        (4, 1): F(1786),
        (3, 1): F(-792),
        (2, 1): F(86),
        (6, 0): F(-8208),
        (5, 0): F(5712),
        (4, 0): F(-1296),
        (3, 0): F(96),
    }
    check("exact/characteristic-polynomial", det == target)
    # Four-point threshold test and cutoff energy.
    h = hs[:3]
    cval = F(16, 5)
    y = [F(0)]
    u = F(1)
    for hi in h:
        y.append(y[-1] + hi * u / (a - hi * (1 - a)))
        u -= hi * cval * y[-1] / a
    check("exact/threshold-principal", y == [F(0), F(4, 5), F(67, 150), -F(3067, 74250)])
    eta = y[:3] + [F(0)]
    energy = sum(
        (a - hi * (1 - a)) * (eta[j + 1] - eta[j]) ** 2 / hi - hi * cval * eta[j + 1] ** 2 / a
        for j, hi in enumerate(h)
    )
    check("exact/threshold-cutoff", energy == -F(205489, 4050000))
    # Genuine conformable integral and two first-order recurrences.
    integral = [F(0)]
    for j in range(2):
        integral.append(F(2, 3) * integral[-1] + F(4, 3))
    check("exact/conformable-integral", integral == [F(0), F(4, 3), F(20, 9)])
    xs, xu = F(0), F(0)
    vals_s = []
    vals_u = []
    for hi in hs[:3]:
        xs = (a - hi * (1 - a)) / (a + hi) * xs + hi / (a + hi)
        xu = (1 - hi * F(5, 3)) * xu + hi * F(4, 3)
        vals_s.append(xs)
        vals_u.append(xu)
    check("exact/first-order-shifted", vals_s == [F(2, 5), F(36, 65), F(41, 65)])
    check("exact/first-order-unshifted", vals_u == [F(2, 3), F(20, 27), F(62, 81)])


def random_identities():
    rng = np.random.default_rng(20261001)
    for trial in range(16):
        h = rng.uniform(0.025, 0.28, 12)
        t = np.r_[0.0, np.cumsum(h)]
        k0 = rng.uniform(0.65, 1.2, 12)
        k1 = rng.uniform(0.02, 0.35, 12)
        nu, E = gauge(t, k0, k1)
        p = rng.uniform(0.7, 1.6, 12)
        q = rng.uniform(-0.5, 0.5, 12)
        y, z = rng.normal(size=(2, 13))
        x = E * y
        xx = E * z
        dy = np.diff(y) / h
        dx = k1 * x[:-1] + k0 * np.diff(x) / h
        dxx = k1 * xx[:-1] + k0 * np.diff(xx) / h
        name = f"variable-grid-{trial:02d}"
        close(name + "/gauge", dx, nu * E[:-1] * dy)
        v = p * dx
        vz = p * dxx
        P = p * nu
        Q = q / k0
        Lx = k1[:-1] * v[:-1] + k0[:-1] * np.diff(v) / h[:-1] + q[:-1] * x[1:-1]
        Ly = np.diff(P * dy) / h[:-1] + Q[:-1] * y[1:-1]
        close(name + "/reduction", Lx, nu[:-1] * E[:-2] * Ly)
        Lz = k1[:-1] * vz[:-1] + k0[:-1] * np.diff(vz) / h[:-1] + q[:-1] * xx[1:-1]
        wr = p * (xx[:-1] * dx - x[:-1] * dxx) / E[:-1] ** 2
        close(
            name + "/Lagrange",
            xx[1:-1] * Lx - x[1:-1] * Lz,
            k0[:-1] * E[1:-1] ** 2 * np.diff(wr) / h[:-1],
        )
        close(
            name + "/energy",
            np.sum(h * (p * dx**2 / (k0 * E[:-1] * E[1:]) - q * x[1:] ** 2 / (k0 * E[1:] ** 2))),
            np.sum(h * (P * dy**2 - Q * y[1:] ** 2)),
        )
        # Completion of the square with arbitrary w, not a Riccati solution.
        w = rng.uniform(-0.15, 0.15, 13)
        D = P + h * w[:-1]
        dw = k1 * w[:-1] + k0 * np.diff(w) / h
        lhs = P * dy**2 - Q * y[1:] ** 2
        rhs = (
            np.diff(w * y * y) / h
            + (p * dx - w[:-1] * x[:-1]) ** 2 / (E[:-1] ** 2 * D)
            - x[1:] ** 2 / (k0 * E[1:] ** 2) * (dw - k1 * w[:-1] + q + k0 * w[:-1] ** 2 / D)
        )
        close(name + "/completion-square", lhs, rhs)
        # Positive solution for a second operator and its Picone/Riccati identities.
        p2 = rng.uniform(0.8, 1.4, 12)
        q2 = -rng.uniform(0.01, 0.3, 12)
        P2 = p2 * nu
        Q2 = q2 / k0
        zz, uu = [1.0], [0.2]
        for j in range(12):
            zz.append(zz[-1] + h[j] * uu[-1] / P2[j])
            uu.append(uu[-1] - h[j] * Q2[j] * zz[-1])
        zz, uu = np.array(zz), np.array(uu)
        Z = E * zz
        zd = k1 * Z[:-1] + k0 * np.diff(Z) / h
        W = uu / zz
        Dp = P2 + h * W[:-1]
        close(name + "/Riccati", np.diff(W) / h + Q2 + W[:-1] ** 2 / Dp, np.zeros(12))
        Wx = Z[:-1] * dx - x[:-1] * zd
        boundary = W * y * y
        rhs = (
            np.diff(boundary) / h
            + (p - p2) * dx**2 / (k0 * E[:-1] * E[1:])
            + (q2 - q) * x[1:] ** 2 / (k0 * E[1:] ** 2)
            + p2 * Wx**2 / (k0 * Z[:-1] * Z[1:] * E[:-1] ** 2)
        )
        close(name + "/Picone", lhs, rhs)
        # The conformable integral is computed independently as a kernel sum.
        f = rng.normal(size=12)

        def integral(values):
            return E * np.r_[0.0, np.cumsum(h * values / (k0 * E[1:]))]

        If = integral(f)
        close(name + "/integral-inverse", k1 * If[:-1] + k0 * np.diff(If) / h, f)
        close(name + "/integral-of-derivative", integral(dx), x - E * x[0])
        ibp = integral(dx * xx[1:]) + integral(x[:-1] * dxx - k1 * x[:-1] * xx[1:])
        close(name + "/integration-by-parts-1", ibp, x * xx - E * x[0] * xx[0])
        ibp = integral(dx * xx[:-1]) + integral(x[1:] * dxx - k1 * x[1:] * xx[:-1])
        close(name + "/integration-by-parts-2", ibp, x * xx - E * x[0] * xx[0])
        reverse = -np.sum(h * f * E[0] / (k0 * E[1:]))
        close(name + "/integral-reversal", If[-1], -E[-1] * reverse)
        Ft = (1 + t[:, None]) ** 2 + np.cos(t[:-1][None, :]) + 0.1 * t[:, None] * t[:-1][None, :]
        It = np.array(
            [E[j] * np.sum(h[:j] * Ft[j, :j] / (k0[:j] * E[1 : j + 1])) for j in range(13)]
        )
        lhs = k1 * It[:-1] + k0 * np.diff(It) / h
        rhs = []
        for j in range(12):
            dF = (Ft[j + 1, :j] - Ft[j, :j]) / h[j]
            rhs.append(nu[j] * E[j] * np.sum(h[:j] * dF / (k0[:j] * E[1 : j + 1])) + Ft[j + 1, j])
        close(name + "/Leibniz", lhs, rhs)
        # Dirichlet Green matrix against the inverse of the ORIGINAL operator.
        qg = -np.ones(12) * 0.1
        da, ea, b, Eg = pencil(t, 0.75, p, qg, k0=k0, k1=k1)
        A = np.diag(da) + np.diag(ea, 1) + np.diag(ea, -1)
        G = -np.linalg.inv(A)
        Ga = E[1:-1, None] * G / (nu[:-1] * E[:-2])[None, :]
        original = []
        for j in range(11):
            bx = np.zeros(13)
            bx[j + 1] = 1.0
            bd = k1 * bx[:-1] + k0 * np.diff(bx) / h
            bv = p * bd
            original.append(k1[:-1] * bv[:-1] + k0[:-1] * np.diff(bv) / h[:-1] + qg[:-1] * bx[1:-1])
        original = np.array(original).T
        close(name + "/Green", np.linalg.inv(original), Ga * h[:-1][None, :])
        check(name + "/Green-sign", np.max(Ga) <= 0.0)


def numerical_tables():
    t5 = harmonic(5)
    base = spectrum(t5, 1.0)
    expected = np.array(
        [
            [4.6927, 17.6282, 44.6791],
            [3.6450, 13.7829, 35.2221],
            [2.7250, 10.4004, 26.8746],
            [1.2674, 5.0238, 13.5088],
            [0.3143, 1.5007, 4.5850],
        ]
    )
    rows = []
    for i, a in enumerate([1.0, 0.9, 0.8, 0.6, 0.4]):
        lam = spectrum(t5, a)
        close(f"tables/five-point-eigenvalues-{a}", lam, expected[i], 5.1e-5 / max(lam))
        ratio = lam / base
        h = np.diff(t5)
        lo = a * (a - h.max() * (1 - a))
        hi = a * (a - h.min() * (1 - a))
        check(f"tables/five-point-ratio-{a}", np.all((ratio >= lo - 1e-12) & (ratio <= hi + 1e-12)))
        rows.append(
            dict(
                alpha=a,
                lambda1=lam[0],
                lambda2=lam[1],
                lambda3=lam[2],
                ratio1=ratio[0],
                ratio2=ratio[1],
                ratio3=ratio[2],
                lower=lo,
                upper=hi,
            )
        )
    write_csv("five_harmonic_spectrum.csv", list(rows[0]), rows)
    scales = [
        ("harmonic15", harmonic(15)),
        ("sqrt25", np.sqrt(np.arange(1, 26))),
        ("alternating17", alternating(17)),
        ("half_uniform", np.arange(13) * 0.5),
        ("integer", np.arange(11)),
    ]
    target = [
        [8.4170, 8.5593, 8.6398],
        [9.3045, 9.4283, 9.4916],
        [9.8050, 9.8106, 9.8129],
        [9.8134] * 3,
        [9.7887] * 3,
        [10.1322, 9.8891, 9.7765],
    ]
    rows = []
    for i, (name, t) in enumerate(scales):
        h = np.diff(t)
        for j, a in enumerate([0.6, 0.8, 1.0]):
            lam = spectrum(t, a)[0]
            product = lam * np.sum(h / (a - h * (1 - a))) * (t[-1] - t[0]) / a
            check(
                f"tables/Lyapunov-{name}-{a}",
                abs(product - target[i][j]) <= 5.1e-5,
                abs(product - target[i][j]),
                5.1e-5,
            )
            rows.append(dict(scale=name, alpha=a, lambda1=lam, product=product))
    for j, a in enumerate([0.6, 0.8, 1.0]):
        root = first_pulse_root(a)
        product = 4 * root / a * (3 / a + 1 / (a - 0.5 * (1 - a)))
        close(f"tables/Lyapunov-pulse-{a}", product, target[-1][j], 5.1e-5 / product)
        blocks = pulse_profile(a, root)
        check(
            f"tables/pulse-principal-positivity-{a}",
            all(np.min(yy[1:-1]) > 0 for _, yy, _ in blocks),
        )
        rows.append(dict(scale="pulse_exact", alpha=a, lambda1=root, product=product))
    write_csv("lyapunov_products.csv", list(rows[0]), rows)
    # Nodal theorem, shifted orthogonality, residuals and order monotonicity.
    rows = []
    for name, t in scales[:3]:
        previous = None
        for a in [0.6, 0.8, 1.0]:
            pot = 4 * np.sin(2 * np.pi * t[:-1])
            lam, Y, X = spectrum(t, a, q=pot, vectors=True)
            d, e, b, E = pencil(t, a, q=pot)
            A = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
            close(f"spectral/{name}/{a}/residual", A @ Y, b[:, None] * Y * lam[None, :])
            close(f"spectral/{name}/{a}/orthogonality", Y.T @ (b[:, None] * Y), np.eye(len(lam)))
            for n in range(len(lam)):
                values = Y[:, n]
                near = abs(values) < 1e-10 * np.max(abs(values))
                nodal = int(np.sum(near))
                nodal += int(np.sum((values[:-1] * values[1:] < 0) & ~near[:-1] & ~near[1:]))
                check(f"spectral/{name}/{a}/nodes-{n + 1}", nodal == n)
            if previous is not None:
                check(f"spectral/{name}/{a}/monotonicity", np.all(lam > previous))
            previous = lam
            for n in range(3):
                rows.append(dict(scale=name, alpha=a, mode=n + 1, eigenvalue=lam[n]))
    write_csv("sign_changing_potential.csv", list(rows[0]), rows)
    # Exact transfer and pulse-mesh threshold.
    ac = brentq(lambda a: pulse_transfer(0.45, a)[0, 1], 0.8, 0.95, xtol=5e-15)
    close("pulse/exact-threshold", ac, 0.8642536653, 5e-11)
    decimal_root = decimal_pulse_threshold()
    close("pulse/independent-80-digit-threshold", ac, float(decimal_root), 3e-15)
    write_csv(
        "pulse_threshold_high_precision.csv",
        ["alpha_critical", "decimal_precision"],
        [dict(alpha_critical=decimal_root, decimal_precision=80)],
    )
    expected = [0.8661312981, 0.8651758663, 0.8647105555, 0.8644810507, 0.8643670922]
    rows = []
    for m, target in zip([20, 40, 80, 160, 320], expected):
        threshold = brentq(lambda a: spectrum(pulse_mesh(m), a, q=0.45)[0], 0.8, 0.95, xtol=1e-13)
        close(f"pulse/mesh-threshold-{m}", threshold, target, 5e-11)
        rows.append(dict(m=m, threshold=threshold, error=threshold - ac, exact=ac))
    write_csv("pulse_thresholds.csv", list(rows[0]), rows)
    # Continuum refinement errors, both discretizations.
    rows = []
    prev = None
    for N in [25, 50, 100, 200, 400, 800]:
        h = 1 / N
        a = 0.7
        exact = a * a * np.pi**2
        ts = 4 * a * (a - h * (1 - a)) / h**2 * np.sin(np.pi * h / 2) ** 2
        cd = 4 * a * a / h**2 * np.sin(np.pi * h / 2) ** 2
        e1, e2 = abs(ts - exact), abs(cd - exact)
        if prev is not None:
            o1, o2 = np.log2(prev[0] / e1), np.log2(prev[1] / e2)
        else:
            o1 = o2 = ""
        close(f"refinement/spectrum-{N}", spectrum(np.linspace(0, 1, N + 1), a)[0], ts, 5e-10)
        rows.append(dict(N=N, error_TS=e1, order_TS=o1, error_CD=e2, order_CD=o2))
        prev = (e1, e2)
    write_csv("refinement.csv", list(rows[0]), rows)
    rows = []
    for a, targets in zip(
        [0.7, 0.8, 0.9], [(0.0118, 0.000724), (0.00650, 0.000403), (0.00392, 0.000244)]
    ):
        values = []
        for m in [25, 50, 100, 200, 400]:
            error = hybrid_mesh_error(a, m)
            values.append(error)
            rows.append(
                dict(
                    alpha=a,
                    m=m,
                    error_max_state=error,
                    order="" if len(values) == 1 else np.log2(values[-2] / error),
                )
            )
        check(f"hybrid/refinement-{a}", values[-1] < values[0] / 15)
        close(f"hybrid/reported-start-{a}", values[0], targets[0], 5.1e-5 if a == 0.7 else 5.1e-6)
        close(f"hybrid/reported-end-{a}", values[-1], targets[1], 6e-7)
    write_csv("hybrid_refinement.csv", list(rows[0]), rows)

    # Variable gains, explicit initial data and output grid.
    def gains(t):
        return (
            0.75 + 0.05 * np.sin(t),
            0.2 + 0.03 * np.cos(t),
            1 + 0.2 * np.cos(t),
            0.4 + 0.1 * np.sin(t),
        )

    def direct(t, S):
        k0, k1, p, q = gains(t)
        x, v = S
        return [(v / p - k1 * x) / k0, (-q * x - k1 * v) / k0]

    def reduced(t, S):
        k0, k1, p, q = gains(t)
        y, u, logE = S
        return [u / (p * k0), -q * y / k0, -k1 / k0]

    ts = np.linspace(0, 4, 801)
    so = solve_ivp(direct, (0, 4), [0, 1], method="DOP853", rtol=1e-12, atol=1e-12, t_eval=ts)
    sr = solve_ivp(reduced, (0, 4), [0, 1, 0], method="DOP853", rtol=1e-12, atol=1e-12, t_eval=ts)
    error = float(np.max(np.abs(so.y - np.exp(sr.y[2]) * sr.y[:2])))
    check("variable-continuum/state-agreement", error < 1e-10, error, 1e-10)
    write_csv(
        "variable_gain_agreement.csv",
        ["max_state_difference", "rtol", "atol", "x0", "v0", "grid_points"],
        [dict(max_state_difference=error, rtol=1e-12, atol=1e-12, x0=0, v0=1, grid_points=801)],
    )
    print(f"Variable-gain max state difference: {error:.6g}", flush=True)


def zero_counts():
    rows = []
    target = {0.3: [6, 21, 67], 0.5: [3, 12, 39], 0.75: [1, 7, 26], 1.0: [1, 5, 19]}
    grid = np.sqrt(np.arange(1, 1000001))
    h = np.diff(grid)
    for a in [0.3, 0.5, 0.75, 1.0]:
        crossings, actual = count_zeros(grid, a, 1 / grid[:-1])
        for j, T in enumerate([10, 100, 1000]):
            # Count contained gaps rather than a crossing that leaves the horizon.
            n = sum(grid[k + 1] <= T for k in crossings) + sum(grid[k] <= T for k in actual)
            check(f"zeros/sqrt-{a}-{T}", n == target[a][j])
            integral = float(np.sum(h[: T * T - 1] / grid[: T * T - 1]) / a)
            rows.append(dict(alpha=a, T=T, zeros=n, weighted_integral=integral))
        print(f"Square-root counts alpha={a}: {target[a]}", flush=True)
    write_csv("sqrt_zero_counts.csv", list(rows[0]), rows)
    rows = []
    mismatch = 0
    t = np.sqrt(np.arange(1, 90001))
    th = np.arange(1.0, 300.25, 0.5)
    for a in np.linspace(0.34, 1.0, 23):
        cr, az = count_zeros(t, a, 0.2 / t[:-1] ** 2)
        ns = len(cr) + len(az)
        cr, az = count_zeros(th, a, 0.2 / th[:-1] ** 2)
        nh = len(cr) + len(az)
        nr = int(np.floor(np.sqrt(max(0.0, 0.2 / a**2 - 0.25)) * np.log(300) / np.pi))
        mismatch += int(ns != nr)
        rows.append(dict(alpha=a, zeros_sqrt=ns, zeros_half_grid=nh, zeros_continuum=nr))
    check("zeros/Euler-six-differences", mismatch == 6)
    write_csv("euler_zero_counts.csv", list(rows[0]), rows)
    rows = []
    with localcontext() as context:
        context.prec = 80
        for aa, target in zip(["0.34", "0.3334", "0.333334", "0.33333334"], [15, 108, 514, 597]):
            a = Decimal(aa)
            h = Decimal(".5")
            P = a - h * (1 - a)
            y = Decimal(0)
            u = Decimal(1)
            n = 0
            for k in range(598):
                t = 1 + h * k
                z = y + h * u / P
                u -= h * Decimal(".2") / (a * t * t) * z
                n += int(z == 0 or y * z < 0)
                scale = max(abs(z), abs(u))
                y, u = z / scale, u / scale
            check("zeros/Euler-80-digit-" + aa, n == target)
            rows.append(dict(alpha=aa, zeros=n, decimal_precision=80))
    write_csv("euler_near_threshold.csv", list(rows[0]), rows)
    # Riccati loss of positivity on the alternating grid.
    t = alternating(19)
    h = np.diff(t)
    rows = []
    for a, tzero, Dzero in [(1.0, 2.65, -9.20), (0.85, 1.85, -0.325), (0.7, 1.6, -0.225)]:
        y, u, cross, actual = propagate(t, a, q=0.4, initial=(1.0, 0.0))
        j = cross[0]
        D = a - h[j] * (1 - a) + h[j] * u[j] / y[j]
        close(f"Riccati/crossing-time-{a}", t[j], tzero)
        check(f"Riccati/first-sign-loss-{a}", D < 0)
        check(f"Riccati/reported-denominator-{a}", abs(D - Dzero) < 0.005)
        rows.append(dict(alpha=a, left_endpoint=t[j], right_endpoint=t[j + 1], w=u[j] / y[j], D=D))
    write_csv("riccati_sign_loss.csv", list(rows[0]), rows)


def run(output):
    """Run all finite checks and write evidence under an explicit output directory."""
    global ROOT
    ROOT = Path(output).resolve()
    ROOT.mkdir(parents=True, exist_ok=True)
    CHECKS.clear()
    exact_examples()
    print("Exact arithmetic checks passed", flush=True)
    random_identities()
    print("Variable-gain finite-grid identities passed", flush=True)
    numerical_tables()
    print("Spectra, transfer roots, convergence passed", flush=True)
    zero_counts()
    summary = dict(
        date=datetime.now(timezone.utc).date().isoformat(),
        python=platform.python_version(),
        dependencies={
            name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")
        },
        passed=True,
        check_count=len(CHECKS),
        categories=dict(Counter(c["name"].split("/")[0] for c in CHECKS)),
        checks=CHECKS,
        scope="Exact finite arithmetic and finite numerical assertions; not a proof of infinite-horizon theorems.",
    )
    (ROOT / "reports").mkdir(parents=True, exist_ok=True)
    (ROOT / "reports" / "verification_results.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(f"{len(CHECKS)} checks passed; data and verification_results.json written.", flush=True)

    return summary
