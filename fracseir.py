#!/usr/bin/env python3
"""
Core solvers: Caputo fractional SEIR, deterministic (product integration)
and stochastic (fractional Euler-Maruyama, vectorized over sample paths).

Scheme (reduces exactly to Euler-Maruyama at q=1):
  y_n = y_0 + (1/G(q+1)) * sum_j w_{n,j} [ h^q f_j + h^{q-1/2} sigma y_j xi_j ]
  w_{n,j} = (n-j)^q - (n-1-j)^q
Short-memory truncation of length M is used for efficiency.
"""
import numpy as np
from math import gamma as G

# ---------------- parameters (R0 ~ 18, nu_c ~ 0.085) ----------------
BETA0 = 2.8
BETA1 = 0.3
LAM   = 5000.0
MU    = 0.005
DEL   = 0.125
GAM   = 0.143
ALP   = 0.001
N0    = int(LAM / MU)          # 1,000,000

R0 = BETA0 * DEL / ((DEL + MU) * (GAM + MU + ALP))
NU_C = MU * (R0 - 1.0)

Y0 = np.array([0.15 * N0, 2000.0, 1000.0, N0 - 0.15 * N0 - 3000.0])


def beta_t(t):
    return BETA0 * (1.0 + BETA1 * np.sin(2.0 * np.pi * t / 52.0))


def drift(Y, t, nu, b0=BETA0, b1=BETA1, lam=LAM, mu=MU,
          de=DEL, ga=GAM, al=ALP):
    """Y: (...,4). Returns drift of same shape."""
    S, E, I, R = Y[..., 0], Y[..., 1], Y[..., 2], Y[..., 3]
    N = np.maximum(S + E + I + R, 1.0)
    b = b0 * (1.0 + b1 * np.sin(2.0 * np.pi * t / 52.0))
    inf = b * S * I / N
    return np.stack([lam - inf - mu * S - nu * S,
                     inf - (de + mu) * E,
                     de * E - (ga + mu + al) * I,
                     ga * I + nu * S - mu * R], axis=-1)


def _weights(q, M):
    """w[k] = k^q - (k-1)^q  for k=1..M  (lag weights)."""
    k = np.arange(1, M + 1, dtype=float)
    return k**q - (k - 1.0)**q


def frac_det(q, nu_fun, T, h, y0=None, M=400, params=None):
    """Deterministic Caputo SEIR via product integration (short memory M)."""
    n = int(round(T / h))
    t = np.linspace(0.0, T, n + 1)
    y0 = Y0.copy() if y0 is None else np.asarray(y0, float)
    Y = np.zeros((n + 1, 4))
    F = np.zeros((n + 1, 4))
    Y[0] = y0
    kw = params or {}
    nu0 = nu_fun(0.0) if callable(nu_fun) else nu_fun
    F[0] = drift(Y[0], 0.0, nu0, **kw)
    w = _weights(q, min(M, n))
    c = h**q / G(q + 1.0)
    for k in range(1, n + 1):
        m = min(k, M)
        # lags 1..m correspond to F[k-1], F[k-2], ...
        idx = np.arange(k - 1, k - 1 - m, -1)
        s = (w[:m, None] * F[idx]).sum(axis=0)
        Y[k] = np.maximum(y0 + c * s, 0.0)
        nu = nu_fun(t[k]) if callable(nu_fun) else nu_fun
        F[k] = drift(Y[k], t[k], nu, **kw)
    return t, Y


def frac_sde(q, nu_fun, T, h, sigma, n_paths, seed=0, y0=None, M=300):
    """
    Vectorized fractional Euler-Maruyama.
    Returns t, Y with Y shape (n+1, n_paths, 4).
    """
    rng = np.random.default_rng(seed)
    n = int(round(T / h))
    t = np.linspace(0.0, T, n + 1)
    y0 = Y0.copy() if y0 is None else np.asarray(y0, float)

    Y = np.zeros((n + 1, n_paths, 4))
    Ginc = np.zeros((n + 1, n_paths, 4))   # increment integrand
    Y[0] = y0

    w = _weights(q, min(M, n))
    cd = h**q / G(q + 1.0)                 # drift scale
    cs = h**(q - 0.5) / G(q + 1.0)         # noise scale

    nu0 = nu_fun(0.0) if callable(nu_fun) else nu_fun
    xi = rng.standard_normal((n_paths, 4))
    Ginc[0] = cd * drift(Y[0], 0.0, nu0) + cs * sigma * Y[0] * xi

    for k in range(1, n + 1):
        m = min(k, M)
        idx = np.arange(k - 1, k - 1 - m, -1)
        s = np.tensordot(w[:m], Ginc[idx], axes=(0, 0))
        Y[k] = np.maximum(y0 + s, 0.0)
        nu = nu_fun(t[k]) if callable(nu_fun) else nu_fun
        xi = rng.standard_normal((n_paths, 4))
        Ginc[k] = cd * drift(Y[k], t[k], nu) + cs * sigma * Y[k] * xi
    return t, Y


if __name__ == "__main__":
    print(f"R0 = {R0:.3f}   nu_c = {NU_C:.4f}   N0 = {N0}")
    for q in (1.0, 0.95, 0.90):
        t, Y = frac_det(q, 0.05, 260, 0.25)
        print(f"  q={q:.2f}  det peak I = {Y[200:,2].max():9.0f}")
    t, Y = frac_sde(0.95, 0.05, 260, 0.25, 0.04, 20, seed=1)
    pk = Y[200:, :, 2].max(axis=0)
    print(f"  q=0.95 stoch: mean peak {pk.mean():.0f}, min {pk.min():.0f}, max {pk.max():.0f}")
