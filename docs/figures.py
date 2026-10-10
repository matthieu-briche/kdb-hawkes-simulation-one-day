"""Figures for the README (Python port of the q simulation, same algorithm).

    python docs/figures.py      # writes assets/fano.png and assets/residuals_cdf.png

Requires numpy and matplotlib. The q scripts do not depend on this file.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent.parent / "assets"
BETA = 50.0
rng = np.random.default_rng(42)


def hawkes_branch(mu, alpha, beta, T):
    """Exponential Hawkes on [0, T) through the cluster representation."""
    n = alpha / beta
    cur = T * rng.random(rng.poisson(mu * T))
    out = [cur]
    while cur.size:
        parents = np.repeat(cur, rng.poisson(n, cur.size))
        ch = parents + rng.exponential(1 / beta, parents.size)
        cur = ch[ch < T]
        out.append(cur)
    return np.sort(np.concatenate(out))


def compensator_increments(t, mu, alpha, beta):
    """Lambda(t_i) - Lambda(t_{i-1}); i.i.d. Exp(1) under the model."""
    d = np.diff(t, prepend=0.0)
    e = np.exp(-beta * d)
    s = np.empty_like(t)  # s_i = sum_{j<=i} exp(-beta (t_i - t_j))
    acc = 0.0
    for i, ei in enumerate(e):
        acc = 1.0 + ei * acc
        s[i] = acc
    prev = np.concatenate(([0.0], s[:-1]))
    return mu * d + alpha / beta * (1 - e) * prev


def fano_theory(n, beta, w):
    kappa = beta * (1 - n)
    return 1 + (1 / (1 - n) ** 2 - 1) * (1 - (1 - np.exp(-kappa * w)) / (kappa * w))


def fano_empirical(t, T, w):
    counts = np.bincount((t // w).astype(int), minlength=int(T // w))[: int(T // w)]
    return counts.var() / counts.mean()


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.25, linewidth=0.6)


def fig_fano():
    T = 2e5
    windows = np.logspace(-3, 1.3, 30)
    fig, ax = plt.subplots(figsize=(7, 4.2), dpi=150)
    colors = {0.6: "#4C78A8", 0.8: "#F58518", 0.9: "#54A24B"}
    for n, col in colors.items():
        t = hawkes_branch(1.0, n * BETA, BETA, T)
        wg = np.logspace(-3.3, 1.5, 200)
        ax.plot(wg, fano_theory(n, BETA, wg), color=col, lw=1.8, label=f"n = {n}: theory F(T)")
        emp = [fano_empirical(t, T, w) for w in windows]
        ax.plot(windows, emp, "o", color=col, ms=3.5, alpha=0.8)
        ax.axhline(1 / (1 - n) ** 2, color=col, lw=0.8, ls=":")
    ax.axhline(1, color="grey", lw=1, ls="--", label="Poisson: F = 1")
    ax.axvline(10, color="black", lw=0.8, alpha=0.5)
    ax.text(10, 1.25, " 10 s window\n used in the checks", fontsize=8, va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("window length T (s)")
    ax.set_ylabel("Fano factor  Var N(T) / E N(T)")
    ax.set_title("Fano factor rises from 1 to the plateau 1/(1-n)²  (dots: simulation)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fano.png")
    plt.close(fig)


def fig_residuals():
    mu, alpha, T = 1.2, 45.0, 23400.0  # TSLA, n = 0.9
    t = hawkes_branch(mu, alpha, BETA, T)
    t_ms = np.round(t * 1e3) / 1e3
    t_jit = np.sort(np.clip(t_ms + 1e-3 * (rng.random(t.size) - 0.5), 0, T))
    cases = [
        ("continuous times", t, "#4C78A8"),
        ("rounded to 1 ms", t_ms, "#E45756"),
        ("1 ms, ties spread within their ms", t_jit, "#54A24B"),
    ]
    fig, ax1 = plt.subplots(figsize=(7, 4.4), dpi=150)
    x = np.linspace(0, 0.15, 300)
    ax1.plot(x, 1 - np.exp(-x), color="black", lw=1.2, label="Exp(1)")
    for label, tt, col in cases:
        r = np.sort(compensator_increments(tt, mu, alpha, BETA))
        ecdf = np.arange(1, r.size + 1) / r.size
        ks = np.max(np.maximum(ecdf - (1 - np.exp(-r)), (1 - np.exp(-r)) - (ecdf - 1 / r.size)))
        print(f"{label:<36} KS = {ks:.4f}")
        lab = f"{label}  (KS = {ks:.4f})"
        ax1.step(r, ecdf, where="post", color=col, lw=1.3, label=lab)
    ax1.set_xlim(0, 0.15)
    ax1.set_ylim(0, 0.2)
    ax1.set_xlabel("compensator residual")
    ax1.set_ylabel("empirical CDF")
    ax1.legend(fontsize=7.5, frameon=False, loc="lower right")
    style(ax1)
    ax1.set_title(
        f"Residual CDF near 0, n = 0.9, {t.size:,} events (KS 5% value {1.358 / np.sqrt(t.size):.4f})\n"
        "1 ms ties put a mass at 0 that Exp(1) does not have",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(OUT / "residuals_cdf.png")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig_fano()
    fig_residuals()
