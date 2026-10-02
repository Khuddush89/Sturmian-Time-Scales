"""Generate the seven research figures and their numerical data on exact domains.

Run from any directory: sturmian figures --output outputs/reference
Solid arcs occur only inside dense blocks. Dotted gap connectors are guides.
Curves as functions of alpha represent a continuous parameter, not time values.
"""

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from .model import (
    alternating,
    count_zeros,
    first_pulse_root,
    gauge,
    harmonic,
    hybrid_exact,
    propagate,
    pulse_mesh,
    pulse_profile,
    pulse_transfer,
    spectrum,
)

ROOT = Path("outputs/reference")
OUT = ROOT / "figures"
COL = ["#374b9c", "#cd3272", "#ee663b", "#e3a21a", "#17949a"]
plt.rcParams.update(
    {
        "font.family": "DejaVu Serif",
        "font.size": 11,
        "axes.titlesize": 11,
        "axes.labelsize": 11,
        "legend.fontsize": 9.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def base(ax):
    ax.grid(alpha=0.2)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(OUT / (name + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def data(name, fields, rows):
    p = ROOT / "data" / name
    with p.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        w.writerows(rows)


def eigenfunctions():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.25), layout="constrained")
    t = harmonic(15)
    rows = []
    for a, c in zip([1.0, 0.85, 0.6], COL):
        lam, Y, X = spectrum(t, a, vectors=True)
        x = np.r_[0.0, X[:, 0], 0.0]
        if x[1] < 0:
            x = -x
        x /= x.max()
        axs[0].vlines(t, 0, x, color=c, alpha=0.7, lw=0.7)
        axs[0].plot(t, x, "o", ms=3, mfc="white", color=c, label=rf"$\alpha={a:g}$")
        rows.extend([["harmonic", a, z, v, lam[0]] for z, v in zip(t, x)])
        blocks = pulse_profile(a)
        scale = max(np.max(x) for _, _, x in blocks)
        for j, (tt, _, xx) in enumerate(blocks):
            axs[1].plot(
                tt, xx / scale, color=c, lw=1.5, label=rf"$\alpha={a:g}$" if j == 0 else None
            )
            axs[1].plot([tt[0], tt[-1]], [xx[0] / scale, xx[-1] / scale], "o", color=c, ms=3)
            if j:
                pt, _, px = blocks[j - 1]
                axs[1].plot([pt[-1], tt[0]], [px[-1] / scale, xx[0] / scale], ":", color=c, lw=1)
            rows.extend([["pulse", a, z, v / scale, first_pulse_root(a)] for z, v in zip(tt, xx)])
    for left in [0, 1.5, 3]:
        axs[1].axvspan(left, left + 1, color="#f3dfae", alpha=0.3)
    for ax, title in zip(
        axs, ["(a) Harmonic scale: isolated nodes", "(b) Pulse scale: exact block transfer"]
    ):
        base(ax)
        ax.set(xlabel="$t$", ylabel="Normalised $x(t)$", title=title, ylim=(-0.04, 1.12))
        ax.legend(loc="upper right")
    save(fig, "fig1_eigenfunctions")
    data("figure1_profiles.csv", ["scale", "alpha", "t", "x_normalised", "lambda"], rows)


def spectral_ratios():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.25), layout="constrained")
    alphas = np.linspace(0.4, 1, 13)
    rows = []
    scales = [
        ("harmonic", harmonic(15)),
        ("square root", np.sqrt(np.arange(1, 26))),
        ("alternating", alternating(17)),
        ("uniform $h=1/2$", np.arange(13) * 0.5),
    ]
    h = np.diff(harmonic(15))
    lo = alphas * (alphas - h.max() * (1 - alphas))
    hi = alphas * (alphas - h.min() * (1 - alphas))
    axs[0].fill_between(alphas, lo, hi, color=COL[1], alpha=0.15, label="Harmonic ratio bounds")
    for (name, t), c in zip(scales, [COL[1], COL[2], COL[3], COL[0]]):
        ref = spectrum(t, 1.0)[:3]
        ratios = np.array([spectrum(t, a)[:3] / ref for a in alphas])
        axs[0].plot(alphas, ratios[:, 0], "-o", color=c, lw=1.3, ms=3, label=name)
        dev = np.max(abs(ratios - ratios[:, :1]), axis=1)
        mask = alphas < 1
        axs[1].semilogy(
            alphas[mask], np.maximum(dev[mask], 1e-17), "-o", color=c, lw=1.3, ms=3, label=name
        )
        rows.extend([[name, a, *r, d] for a, r, d in zip(alphas, ratios, dev)])
    axs[0].set(title="(a) First-mode spectral ratios", ylabel=r"$\lambda_1(\alpha)/\lambda_1(1)$")
    axs[1].set(
        title="(b) Difference across the first three modes", ylabel="Maximum ratio difference"
    )
    for ax in axs:
        base(ax)
        ax.set_xlabel(r"Order $\alpha$")
        ax.legend(loc="upper left")
    save(fig, "fig2_spectrum")
    data(
        "figure2_ratios.csv",
        ["scale", "alpha", "ratio1", "ratio2", "ratio3", "max_deviation"],
        rows,
    )


def euler():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.25), layout="constrained")
    inds = np.unique(np.r_[np.geomspace(1, 89999, 180).astype(int), 89999])
    t = np.sqrt(inds)
    mu = np.sqrt(inds + 1) - t
    axs[0].loglog(t, mu, "o", mfc="white", ms=3, color=COL[2], label=r"$\{\sqrt{k}\}$")
    t2 = np.arange(1.0, 300.0, 0.5)
    axs[0].loglog(
        t2, np.full_like(t2, 0.5), "s", mfc="white", ms=2.5, color=COL[1], label="$h=1/2$"
    )
    axs[0].axhline(1e-3, color=COL[4], ls="--", label=r"$\mathbb{R}$: $\mu=0$ (axis floor)")
    axs[0].set(title="(a) Graininess on the three domains", xlabel="$t$", ylabel=r"$\mu(t)$")
    with (ROOT / "data" / "euler_zero_counts.csv").open() as f:
        rows = list(csv.DictReader(f))
    a = np.array([float(r["alpha"]) for r in rows])
    for key, c, marker, label in [
        ("zeros_half_grid", COL[1], "s", "$h=1/2$"),
        ("zeros_sqrt", COL[2], "o", r"$\{\sqrt{k}\}$"),
    ]:
        axs[1].plot(a, [int(r[key]) for r in rows], marker, ms=4, mfc="white", color=c, label=label)
    axs[1].step(
        a,
        [int(r["zeros_continuum"]) for r in rows],
        where="post",
        color=COL[4],
        lw=1.4,
        label=r"$\mathbb{R}$, exact count",
    )
    axs[1].axvline(2 * np.sqrt(0.2), color=COL[4], ls="--", lw=1)
    axs[1].axvline(1 / 3, color=COL[1], ls=":", lw=1)
    axs[1].set(
        title="(b) Generalized zeros in $(1,300]$",
        xlabel=r"Order $\alpha$",
        ylabel="Number of zeros",
        xlim=(0.32, 1.02),
    )
    for ax in axs:
        base(ax)
        ax.legend(loc="upper right")
    save(fig, "fig3_euler")


def hybrid():
    fig, axs = plt.subplots(3, 1, figsize=(9, 4.85), sharex=True, layout="constrained")
    rows = []
    for ax, a, c in zip(axs, [0.9, 0.8, 0.7], COL):
        nodes, blocks = hybrid_exact(a)
        scale = max(max(abs(states[0])) for _, states, _ in blocks)
        for k, (tt, states, fun) in enumerate(blocks):
            ax.axvspan(tt[0], tt[-1], color="#f3dfae", alpha=0.3)
            ax.plot(tt, states[0] / scale, color=c, lw=1.4)
            idx = np.flatnonzero(states[0, :-1] * states[0, 1:] < 0)
            for j in idx:
                root = brentq(lambda t: fun(t)[0], tt[j], tt[j + 1])
                ax.plot(root, 0, "+", color=c, ms=8, mew=1.3)
            rows.extend([[a, z, x / scale, v] for z, x, v in zip(tt, states[0], states[1])])
        tn = np.array(sorted(nodes))
        xn = np.array([nodes[t][0] / scale for t in tn])
        ax.plot(tn, xn, "o", color=c, ms=3)
        for j in range(len(tn) - 1):
            if tn[j + 1] - tn[j] > 0.75:
                ax.plot(tn[j : j + 2], xn[j : j + 2], ":", color=c, lw=1)
                if xn[j] * xn[j + 1] < 0:
                    ax.plot(tn[j], 0, "x", color=c, ms=6, mew=1.3)
        ax.axhline(0, color="gray", lw=0.7)
        ax.set(ylabel=rf"$\alpha={a:g}$", ylim=(-1.1, 1.1))
        base(ax)
    axs[0].set_title("Hybrid domain: solid arcs on dense blocks, dotted guides across gaps")
    axs[-1].set_xlabel("$t$")
    save(fig, "fig4_hybrid")
    data("figure4_dense_blocks.csv", ["alpha", "t", "x_normalised", "flux"], rows)


def square_root():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.25), layout="constrained")
    rows = []
    short = np.sqrt(np.arange(1, 101))
    for a, c in [(1.0, COL[0]), (0.3, COL[3])]:
        y, u, cross, actual = propagate(short, a, q=1 / short[:-1])
        profile = y / short**0.25
        profile /= max(abs(profile))
        axs[0].plot(
            short,
            profile,
            "o",
            mfc="white",
            ms=2.6,
            color=c,
            label=rf"$\alpha={a:g}$ ({len(cross) + len(actual)} zeros)",
        )
        axs[0].plot(short[cross], np.zeros(len(cross)), "x", color=c, ms=6, mew=1.2)
        rows.extend([[a, t, v] for t, v in zip(short, profile)])
    grid = np.sqrt(np.arange(1, 1000001))
    horizons = np.sqrt(np.unique(np.geomspace(1, 1000000, 130).astype(int)))
    countrows = []
    for a, c in [(0.3, COL[3]), (0.5, COL[2]), (0.75, COL[1]), (1.0, COL[0])]:
        cr, az = count_zeros(grid, a, 1 / grid[:-1])
        endpoints = np.sort(np.r_[grid[np.array(cr, dtype=int) + 1], grid[np.array(az, dtype=int)]])
        cnt = np.searchsorted(endpoints, horizons, side="right")
        axs[1].semilogx(horizons, cnt, "o", mfc="white", ms=2.7, color=c, label=rf"$\alpha={a:g}$")
        phase = 2 * (np.sqrt(horizons) - 1) / (np.pi * a)
        axs[1].semilogx(horizons, phase, ":", color=c, lw=1)
        countrows.extend([[a, t, n, p] for t, n, p in zip(horizons, cnt, phase)])
    axs[0].axhline(0, color="gray", lw=0.7)
    axs[0].set(
        title="(a) Gauged solution at the isolated nodes",
        xlabel="$t$",
        ylabel="$y(t)/t^{1/4}$, normalised",
    )
    axs[1].set(
        title="(b) Counts and continuum phase heuristic",
        xlabel="Horizon $T$",
        ylabel="Generalized-zero count",
    )
    for ax in axs:
        base(ax)
    axs[0].legend(loc="upper right")
    axs[1].legend(loc="upper left")
    save(fig, "fig5_sqrt")
    data("figure5_short_profile.csv", ["alpha", "t", "y_scaled_normalised"], rows)
    data("figure5_cumulative_counts.csv", ["alpha", "T", "zeros", "phase_heuristic"], countrows)


def riccati():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.25), layout="constrained")
    t = alternating(19)
    h = np.diff(t)
    rows = []
    for a, c in zip([1.0, 0.85, 0.7], COL):
        y, u, cross, actual = propagate(t, a, q=0.4, initial=(1.0, 0.0))
        nu, E = gauge(t, a, 1 - a)
        x = E * y
        j = cross[0]
        axs[0].plot(t, x, "o", mfc="white", color=c, ms=3, label=rf"$\alpha={a:g}$")
        axs[0].plot(t, x, ":", color=c, lw=0.8)
        axs[0].plot(t[j], 0, "x", color=c, ms=7, mew=1.3)
        D = nu + h * u[:-1] / y[:-1]
        axs[1].plot(
            t[: j + 1], D[: j + 1], ":o", mfc="white", color=c, ms=3, label=rf"$\alpha={a:g}$"
        )
        axs[1].plot(t[j], D[j], "v", color=c, ms=6)
        rows.extend([[a, t[k], x[k], u[k] / y[k], D[k]] for k in range(len(h))])
    axs[0].set(title="(a) First generalized zero", xlabel="$t$", ylabel="$x(t)$")
    axs[1].set(
        title="(b) Riccati positivity fails at that gap",
        xlabel="$t$",
        ylabel=r"$\mathcal{D}=P+\mu w$",
    )
    axs[1].set_yscale("symlog", linthresh=0.35)
    for ax in axs:
        base(ax)
        ax.axhline(0, color="gray", lw=0.8)
        ax.legend(loc="lower left")
    save(fig, "fig6_riccati")
    data("figure6_riccati.csv", ["alpha", "t", "x", "w", "D"], rows)


def disconjugacy():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.25), layout="constrained")
    rows = []
    aa = np.linspace(0.36, 1, 24)
    mesh = pulse_mesh(40)
    exact = np.array([first_pulse_root(a) - 0.45 for a in aa])
    grid = np.array([spectrum(mesh, a, q=0.45)[0] for a in aa])
    axs[0].plot(aa, exact, "--", color=COL[0], lw=1.3, label="Exact first transfer parameter")
    axs[0].plot(aa, grid, "o-", mfc="white", color=COL[2], ms=3, lw=0.9, label="Pulse mesh, $m=40$")
    ac = brentq(lambda a: pulse_transfer(0.45, a)[0, 1], 0.8, 0.95)
    axs[0].axvline(ac, color=COL[1], ls=":", lw=1)
    axs[0].axhline(0, color="gray", lw=0.7)
    axs[0].axvspan(0.33334, ac, color=COL[1], alpha=0.09)
    axs[0].axvspan(ac, 1.0, color=COL[3], alpha=0.09)
    axs[0].set(
        title="(a) Closed-interval disconjugacy threshold",
        xlabel=r"Order $\alpha$",
        ylabel="First Dirichlet parameter",
        xlim=(0.34, 1.015),
    )
    axs[0].legend(loc="upper left")
    rows.extend([[a, g, e] for a, g, e in zip(aa, grid, exact)])
    tt = np.arange(11)
    axs[1].axvspan(0, 1, color=COL[1], alpha=0.13, label="Shared sign-changing gap")
    for pair, c, offset, label in [
        ((1.0, -1.0), COL[0], 0.16, "Initial pair $(1,-1)$"),
        ((1.0, -2.0), COL[2], -0.16, "Initial pair $(1,-2)$"),
    ]:
        vals = [*pair]
        for k in range(9):
            vals.append(-vals[-1] - vals[-2])
        vals = np.array(vals)
        axs[1].plot(tt, vals, ":o", color=c, mfc="white", ms=3, lw=0.8, label=label)
        jj = np.flatnonzero(vals[:-1] * vals[1:] < 0)
        axs[1].plot(tt[jj], np.full(len(jj), offset), "x", color=c, ms=6, mew=1.2)
    axs[1].axhline(0, color="gray", lw=0.7)
    axs[1].set(
        title="(b) Independent solutions may share a gap",
        xlabel=r"$k\in\mathbb{Z}$",
        ylabel="$y_k$",
        ylim=(-2.3, 2.25),
    )
    axs[1].legend(loc="upper right")
    for ax in axs:
        base(ax)
    save(fig, "fig7_disconjugacy")
    data("figure7_threshold_curve.csv", ["alpha", "mesh_lambda1", "exact_transfer_parameter"], rows)


def run(output):
    """Create PNG figures and CSV data; verification supplies the Euler count table."""
    global ROOT, OUT
    ROOT = Path(output).resolve()
    OUT = ROOT / "figures"
    OUT.mkdir(parents=True, exist_ok=True)
    (ROOT / "data").mkdir(parents=True, exist_ok=True)
    if not (ROOT / "data" / "euler_zero_counts.csv").exists():
        raise FileNotFoundError("Run sturmian verify with the same --output first.")
    for fun in [eigenfunctions, spectral_ratios, euler, hybrid, square_root, riccati, disconjugacy]:
        fun()
        print(fun.__name__ + " complete", flush=True)
