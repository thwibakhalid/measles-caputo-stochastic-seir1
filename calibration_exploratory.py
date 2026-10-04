#!/usr/bin/env python3
"""
Calibration v2 — corrected demography.

Fixes relative to v1:
  * birth rate decoupled from death rate (Sudan CBR ~33/1000/yr, not 1/64/yr)
  * routine MCV coverage removes a fraction pv of newborns from S
  * reporting fraction rho profiled analytically AND constrained to (0,1]
  * beta0 bounds widened so R0 can reach measles-realistic values
"""
import numpy as np, pandas as pd, time, sys
from math import gamma as G
from scipy.optimize import differential_evolution

from pathlib import Path
CSV = Path(__file__).resolve().parent.parent / "data" / "sudan_measles_monthly_2012_2022.csv"
COL, WPM = "measles_total", 52.0/12.0

DELTA = 0.125            # 1/8 d
GAMMA = 0.143            # 1/7 d
MU    = 1.0/(64*52)      # death rate, per week
BIRTH = 0.033/52         # Sudan crude birth rate, per week
ALPHA = 0.002
NPOP  = 4.0e7

def load_sudan(y0=2012, y1=2022):
    s = pd.read_csv(CSV, parse_dates=["date"]).rename(columns={"date": "t"})
    s = s[(s.Year >= y0) & (s.Year <= y1)].sort_values("t")
    full = pd.date_range(s.t.min(), s.t.max(), freq="MS")
    s = s.set_index("t").reindex(full)
    return full, s[COL].interpolate(limit_direction="both").to_numpy(float)

DATES, YOBS = load_sudan()
NM   = len(YOBS)
SYY  = float(np.sum(YOBS**2))
SST  = float(np.sum((YOBS-YOBS.mean())**2))
H    = 1.0
NSTEP= int(round(NM*WPM/H))
TGRID= np.linspace(0.0, NM*WPM, NSTEP+1)
MIDX = np.minimum((TGRID/WPM).astype(int), NM-1)

def simulate_E(q, b0, b1, ph, pv, S0f, E0f, M=160):
    """Unit (rho=1) monthly incidence delta*E."""
    n = NSTEP
    k = np.arange(1, min(M, n)+1, dtype=float)
    w = k**q - (k-1.0)**q
    c = H**q / G(q+1.0)
    S0, E0 = S0f*NPOP, E0f*NPOP
    I0 = 0.3*E0
    Y0 = np.array([S0, E0, I0, NPOP-S0-E0-I0])
    Y = np.empty((n+1,4)); F = np.empty((n+1,4))
    Y[0] = Y0
    tp = 2*np.pi/52.0
    recruit = BIRTH*NPOP*(1.0-pv)          # susceptible newborns

    def f(y, t):
        S,E,I,R = y
        N = S+E+I+R
        if N < 1.0: N = 1.0
        b = b0*(1.0 + b1*np.sin(tp*(t-ph)))
        inf = b*S*I/N
        return np.array([recruit - inf - MU*S,
                         inf - (DELTA+MU)*E,
                         DELTA*E - (GAMMA+MU+ALPHA)*I,
                         GAMMA*I + BIRTH*NPOP*pv - MU*R])

    F[0] = f(Y0, 0.0)
    for j in range(1, n+1):
        m = min(j, M)
        sl = slice(j-1, j-1-m if j-1-m >= 0 else None, -1)
        Y[j] = Y0 + c*(w[:m] @ F[sl])
        np.maximum(Y[j], 0.0, out=Y[j])
        F[j] = f(Y[j], TGRID[j])

    out = np.zeros(NM)
    np.add.at(out, MIDX, DELTA*Y[:,1]*H)
    return out

def sse_profiled(p, q):
    try:
        pred = simulate_E(q, *p)
        if not np.all(np.isfinite(pred)): return 1e18
        spp = float(pred @ pred)
        if spp <= 0: return 1e18
        syp = float(YOBS @ pred)
        if syp <= 0: return SYY
        rho = syp/spp
        if rho > 1.0:                      # infeasible: >100% reporting
            rho = 1.0                      # clamp and evaluate honestly
            return float(np.sum((pred-YOBS)**2))
        return SYY - syp*syp/spp
    except Exception:
        return 1e18

def rho_of(p, q):
    pred = simulate_E(q, *p)
    return min(float(YOBS @ pred)/float(pred @ pred), 1.0), pred

BOUNDS = [(1.0, 9.0),      # beta0  -> R0 up to ~60
          (0.05, 0.95),    # beta1
          (0.0, 52.0),     # phase
          (0.30, 0.95),    # pv  routine MCV coverage
          (0.01, 0.25),    # S0 fraction
          (2e-6, 3e-3)]    # E0 fraction

def fit(q, seed=5, maxiter=70, popsize=14):
    t0 = time.time()
    r = differential_evolution(sse_profiled, BOUNDS, args=(q,), seed=seed,
                               maxiter=maxiter, popsize=popsize, tol=1e-10,
                               mutation=(0.4,1.0), recombination=0.9,
                               polish=True, init='sobol')
    return r.x, r.fun, time.time()-t0

def report(q, p, S, dt):
    rho, pred = rho_of(p, q)
    rmse = np.sqrt(S/NM); r2 = 1.0 - S/SST
    k = 7 if q == 1.0 else 8
    aicc = NM*np.log(S/NM) + 2*k + (2*k*(k+1))/(NM-k-1)
    R0 = p[0]*DELTA/((DELTA+MU)*(GAMMA+MU+ALPHA))
    nm = ["beta0","beta1","phase","pv","S0f","E0f"]
    at = [n for n,v,(lo,hi) in zip(nm,p,BOUNDS)
          if abs(v-lo) < 1e-3*(hi-lo) or abs(v-hi) < 1e-3*(hi-lo)]
    print(f"q={q:.2f} SSE={S:.6g} RMSE={rmse:.1f} R2={r2:.4f} AICc={aicc:.2f} "
          f"rho={rho:.4f} R0={R0:.1f} t={dt:.0f}s")
    print("   " + "  ".join(f"{n}={v:.4g}" for n,v in zip(nm,p)))
    print("   at-bound: " + (", ".join(at) if at else "none"), flush=True)
    np.savez(f"v2_q{q:.2f}.npz", p=p, sse=S, q=q, rho=rho, pred=pred,
             yobs=YOBS, dates=np.array([str(d) for d in DATES]), r2=r2, aicc=aicc)
    return r2, rho, aicc

if __name__ == "__main__":
    for q in [float(x) for x in sys.argv[1:]]:
        p, S, dt = fit(q)
        report(q, p, S, dt)
