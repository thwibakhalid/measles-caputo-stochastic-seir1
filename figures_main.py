#!/usr/bin/env python3
"""Six essential figures, all at Caputo fractional order."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import rankdata, spearmanr
from math import gamma as G
import warnings; warnings.filterwarnings('ignore')

from fracseir import (frac_det, frac_sde, drift, BETA0, BETA1, LAM, MU,
                      DEL, GAM, ALP, N0, R0, NU_C, Y0)

plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                     'axes.labelsize': 12, 'axes.titlesize': 11,
                     'figure.dpi': 300, 'savefig.dpi': 300,
                     'savefig.bbox': 'tight', 'axes.grid': True,
                     'grid.alpha': 0.3})

QREF = 0.95          # representative fractional order used throughout
print(f"R0={R0:.2f}  nu_c={NU_C:.4f}  q_ref={QREF}")

# =====================================================================
# FIGURE A — memory effect of the fractional order
# =====================================================================
print("\n[A] memory effect ...")
figA, ax = plt.subplots(1, 2, figsize=(13, 4.8))

for q, c in zip([0.85, 0.90, 0.95, 1.00],
                ['#e74c3c', '#e67e22', '#2980b9', '#2c3e50']):
    t, Y = frac_det(q, 0.05, 300, 0.25)
    ax[0].plot(t, Y[:, 2], color=c, lw=1.8, label=f'$q={q:.2f}$')
ax[0].set_xlabel('Time (weeks)'); ax[0].set_ylabel('Infectious $I(t)$')
ax[0].set_title('(a) Caputo memory reshapes the outbreak')
ax[0].legend(framealpha=0.9)

from scipy.signal import find_peaks
qs = np.linspace(0.80, 1.00, 11)
pk, tp = [], []
for q in qs:
    t, Y = frac_det(q, 0.05, 300, 0.25)
    I = Y[:, 2]
    loc, _ = find_peaks(I, prominence=I.max() * 0.02)
    j = int(loc[0]) if len(loc) else int(np.argmax(I))
    pk.append(I[j]); tp.append(t[j])   # first epidemic peak
l1 = ax[1].plot(qs, np.array(pk) / 1e3, 'o-', color='#8e44ad', lw=2,
                label='First peak $I_{\\max}$ (thousands)')
ax[1].set_xlabel('Fractional order $q$')
ax[1].set_ylabel('First peak $I_{\\max}$ (thousands)', color='#8e44ad')
ax[1].tick_params(axis='y', labelcolor='#8e44ad')
axb = ax[1].twinx()
l2 = axb.plot(qs, tp, 's--', color='#16a085', lw=2, label='Time to first peak (wk)')
axb.set_ylabel('Time to first peak (weeks)', color='#16a085')
axb.tick_params(axis='y', labelcolor='#16a085')
axb.grid(False)
ax[1].set_title('(b) First-peak magnitude and timing versus $q$')
ax[1].legend(l1 + l2, [x.get_label() for x in l1 + l2],
             framealpha=0.9, loc='center right', fontsize=9)
figA.tight_layout(); figA.savefig('FIGA_memory.pdf'); plt.close(figA)
print(f"   peaks {pk[0]:.0f} -> {pk[-1]:.0f}")

# =====================================================================
# FIGURE B — bifurcation / critical threshold (fractional)
# =====================================================================
print("[B] bifurcation ...")
figB, ax = plt.subplots(1, 2, figsize=(13, 4.8))

nus = np.linspace(0.005, 0.16, 32)
det_pk, sto_pk = [], []
for nu in nus:
    t, Y = frac_det(QREF, nu, 300, 0.5)
    det_pk.append(Y[150:, 2].max())
    t, Ys = frac_sde(QREF, nu, 300, 0.5, 0.05, 25, seed=int(nu * 1e4))
    sto_pk.append(Ys[150:, :, 2].max(axis=0))
det_pk = np.array(det_pk)

for i, nu in enumerate(nus):
    ax[0].scatter([nu] * len(sto_pk[i]), np.array(sto_pk[i]) / 1e3,
                  c='gray', alpha=0.18, s=7, zorder=1)
ax[0].plot(nus, det_pk / 1e3, 'b-', lw=2.5, zorder=3,
           label=f'Deterministic ($q={QREF}$)')
ax[0].axvline(NU_C, color='red', ls='--', lw=1.6,
              label=f'$\\nu_c={NU_C:.3f}$')
ax[0].set_xlabel('Vaccination rate $\\nu$ (wk$^{-1}$)')
ax[0].set_ylabel('Peak infections (thousands)')
ax[0].set_title('(a) Bifurcation diagram')
ax[0].legend(fontsize=9, framealpha=0.9)

for nu, c, lab in [(0.06, '#e74c3c', '$\\nu=0.06<\\nu_c$'),
                   (0.085, '#f39c12', '$\\nu\\approx\\nu_c$'),
                   (0.12, '#27ae60', '$\\nu=0.12>\\nu_c$')]:
    t, Y = frac_det(QREF, nu, 300, 0.25)
    ax[1].plot(t, Y[:, 2] / 1e3, color=c, lw=1.8, label=lab)
ax[1].set_xlabel('Time (weeks)'); ax[1].set_ylabel('Infectious (thousands)')
ax[1].set_title(f'(b) Regimes across $\\nu_c$ ($q={QREF}$)')
ax[1].legend(fontsize=9, framealpha=0.9)
figB.tight_layout(); figB.savefig('FIGB_bifurcation.pdf'); plt.close(figB)

# =====================================================================
# FIGURE C — vaccination strategies (fractional)
# =====================================================================
print("[C] strategies ...")
figC, ax = plt.subplots(1, 2, figsize=(13, 4.8))

def nu_const(t):  return 0.06
def nu_pw(t):     return 0.15 if (t % 52) < 6 else 0.045
def nu_pulse(t):  return 0.78 if (t % 52) < 4 else 0.0
def nu_zero(t):   return 0.0

strat = [('Constant ($\\nu=0.06$)', nu_const, '#2980b9'),
         ('Piecewise constant',     nu_pw,    '#27ae60'),
         ('Annual pulse',           nu_pulse, '#8e44ad'),
         ('No vaccination',         nu_zero,  '#c0392b')]

tt = np.linspace(0, 104, 800)
for nm, f, c in strat[:3]:
    ax[0].plot(tt, [f(x) for x in tt], color=c, lw=1.6, label=nm, alpha=.85)
ax[0].set_xlabel('Time (weeks)'); ax[0].set_ylabel('$\\nu(t)$ (wk$^{-1}$)')
ax[0].set_title('(a) Vaccination schedules')
ax[0].legend(fontsize=9, framealpha=0.9)

rows = []
for nm, f, c in strat:
    t, Y = frac_det(QREF, f, 300, 0.25)
    I = Y[:, 2]
    ax[1].plot(t, I / 1e3, color=c, lw=1.7, label=nm)
    doses = np.trapezoid(np.array([f(x) for x in t]) * Y[:, 0], t)
    rows.append((nm, I[150:].max(), doses))
ax[1].set_xlabel('Time (weeks)'); ax[1].set_ylabel('Infectious (thousands)')
ax[1].set_title(f'(b) Response under each schedule ($q={QREF}$)')
ax[1].legend(fontsize=8, framealpha=0.9)
figC.tight_layout(); figC.savefig('FIGC_strategies.pdf'); plt.close(figC)

base_pk, base_ds = rows[0][1], rows[0][2]
print("   strategy table:")
for nm, p, d in rows:
    eff = (base_pk / p) * (base_ds / d) if d > 0 else float('nan')
    print(f"     {nm:26s} peak={p:9.0f} doses={d:.3e} eff={eff:.2f}")
np.save('strategy_rows.npy', np.array([[r[1], r[2]] for r in rows]))

# =====================================================================
# FIGURE D — stochastic vs deterministic (fractional)
# =====================================================================
print("[D] stochastic ...")
figD, ax = plt.subplots(1, 3, figsize=(16, 4.8))
NU_D, SIG, T_D, H_D = 0.05, 0.04, 300, 0.25

t, Yd = frac_det(QREF, NU_D, T_D, H_D)
ax[0].plot(t, Yd[:, 2] / 1e3, 'b-', lw=2.2, zorder=5, label='Deterministic')
t, Ys = frac_sde(QREF, NU_D, T_D, H_D, SIG, 20, seed=11)
for p in range(20):
    ax[0].plot(t, Ys[:, p, 2] / 1e3, color='gray', alpha=.18, lw=.6)
ax[0].plot(t, Ys[:, 3, 2] / 1e3, color='red', alpha=.8, lw=1.0,
           label='Stochastic path')
ax[0].set_xlabel('Time (weeks)'); ax[0].set_ylabel('Infectious (thousands)')
ax[0].set_title(f'(a) Sample paths ($q={QREF}$, $\\sigma={SIG}$)')
ax[0].legend(fontsize=9, framealpha=0.9)

t, Ys = frac_sde(QREF, NU_D, T_D, H_D, SIG, 400, seed=202)
peaks = Ys[150:, :, 2].max(axis=0)
dpk = Yd[150:, 2].max()
fades = int((Ys[-40:, :, 2].max(axis=0) < 50).sum())
ax[1].hist(peaks / 1e3, bins=28, color='steelblue', alpha=.75,
           edgecolor='white')
ax[1].axvline(dpk / 1e3, color='red', ls='--', lw=2,
              label=f'Det. {dpk/1e3:.1f}k')
ax[1].axvline(peaks.mean() / 1e3, color='green', lw=2,
              label=f'Mean {peaks.mean()/1e3:.1f}k')
ax[1].set_xlabel('Peak infections (thousands)'); ax[1].set_ylabel('Count')
ax[1].set_title(f'(b) Peak distribution ($n=400$, fade-outs {fades})')
ax[1].legend(fontsize=9, framealpha=0.9)

sigs = [0.01, 0.03, 0.05, 0.08, 0.12]
box = []
for s in sigs:
    _, Yb = frac_sde(QREF, NU_D, T_D, H_D, s, 80, seed=int(s * 1e4))
    box.append(Yb[150:, :, 2].max(axis=0) / 1e3)
bp = ax[2].boxplot(box, labels=[str(s) for s in sigs], patch_artist=True,
                   medianprops=dict(color='red', lw=2))
for pt, c in zip(bp['boxes'], ['#27ae60', '#2980b9', '#f39c12',
                               '#e67e22', '#e74c3c']):
    pt.set_facecolor(c); pt.set_alpha(.6)
ax[2].axhline(dpk / 1e3, color='blue', ls='--', lw=1.4, label='Deterministic')
ax[2].set_xlabel('Noise intensity $\\sigma$')
ax[2].set_ylabel('Peak infections (thousands)')
ax[2].set_title('(c) Variability versus noise')
ax[2].legend(fontsize=9, framealpha=0.9)
figD.tight_layout(); figD.savefig('FIGD_stochastic.pdf'); plt.close(figD)
print(f"   det={dpk:.0f} mean={peaks.mean():.0f} fades={fades}/400")

# =====================================================================
# FIGURE E — extinction probability + corrected threshold
# =====================================================================
print("[E] extinction ...")
figE, ax = plt.subplots(1, 2, figsize=(13, 4.8))

NU_E = 0.083                       # just below nu_c
S_dfe = LAM / (MU + NU_E)
y0e = np.array([S_dfe, 300.0, 150.0, N0 - S_dfe - 450.0])
sig_grid = np.arange(0.02, 0.21, 0.02)
pext = []
for s in sig_grid:
    _, Ye = frac_sde(QREF, NU_E, 300, 0.5, s, 150, seed=int(s * 1e5),
                     y0=y0e)
    pext.append(float((Ye[-60:, :, 2].max(axis=0) < 20).mean()))
    print(f"     sigma={s:.2f}  P(ext)={pext[-1]:.2f}")

ax[0].plot(sig_grid, pext, 'o-', color='#8e44ad', lw=2.4, ms=7,
           mfc='white', mew=2)
ax[0].fill_between(sig_grid, pext, alpha=.12, color='purple')
ax[0].axhline(0.5, color='gray', ls='--', alpha=.5, label='50%')
ax[0].set_ylim(-0.04, 1.04)
ax[0].set_xlabel('Noise intensity $\\sigma$')
ax[0].set_ylabel('Extinction probability')
ax[0].set_title(f'(a) Fade-out probability ($\\nu={NU_E}$, $q={QREF}$)')
ax[0].legend(fontsize=9, framealpha=0.9)

sg = np.linspace(0.001, 0.40, 300)
R_old = R0 * np.exp(-0.5 * 3 * sg**2)
R_new = np.array([BETA0 * DEL /
                  ((DEL + MU + .5 * s**2) * (GAM + MU + ALP + .5 * s**2)) *
                  MU / (MU + NU_E + .5 * s**2) for s in sg])
R_nu = R0 * MU / (MU + NU_E)
ax[1].axhline(R_nu, color='gray', ls='-.', lw=1.5,
              label=f'$\\mathcal{{R}}_\\nu={R_nu:.3f}$')
ax[1].plot(sg, R_old, 'r-', lw=2.3, alpha=.75,
           label='Superseded: $R_0e^{-\\sigma^2/2}$')
ax[1].plot(sg, R_new, 'b-', lw=2.5, label='Corrected $\\mathcal{R}^s_0$')
ax[1].axhline(1.0, color='green', ls=':', lw=2.2, label='Threshold $=1$')
ax[1].axhline(R0, color='k', ls='--', lw=1.0, alpha=.45,
              label=f'$\\mathcal{{R}}_0={R0:.1f}$ (no vaccination)')
ax[1].set_yscale('log')
ax[1].set_ylim(0.1, 40); ax[1].set_xlim(0, 0.40)
ax[1].set_xlabel('Noise intensity $\\sigma$')
ax[1].set_ylabel('Reproduction number')
ax[1].set_title(f'(b) Stochastic thresholds ($\\nu={NU_E}$)')
ax[1].legend(fontsize=8, loc='center right', framealpha=0.9)
figE.tight_layout(); figE.savefig('FIGE_extinction.pdf'); plt.close(figE)

# =====================================================================
# FIGURE F — PRCC global sensitivity (fractional)
# =====================================================================
print("[F] PRCC ...")
names = ['$\\beta_0$', '$\\beta_1$', '$\\delta$', '$\\gamma$',
         '$\\mu$', '$\\alpha$', '$\\nu$', '$q$']
base = np.array([2.8, 0.3, 0.125, 0.143, 0.005, 0.001, 0.05, 0.95])
NL, NP = 300, len(base)
rng = np.random.default_rng(7)
X = np.zeros((NL, NP))
for j in range(NP):
    e = np.linspace(0, 1, NL + 1)
    p = rng.uniform(e[:-1], e[1:]); rng.shuffle(p)
    lo, hi = base[j] * 0.85, min(base[j] * 1.15, 1.0 if j == 7 else 1e9)
    X[:, j] = lo + p * (hi - lo)

Ypk = np.full(NL, np.nan); Ytp = np.full(NL, np.nan)
for i in range(NL):
    b0, b1, de, ga, mu, al, nu, q = X[i]
    lam = mu * N0
    pr = dict(b0=b0, b1=b1, lam=lam, mu=mu, de=de, ga=ga, al=al)
    try:
        t, Y = frac_det(q, nu, 200, 0.5, params=pr)
        I = Y[:, 2]
        if np.all(np.isfinite(I)) and I.max() > 100:
            j = 100 + int(np.argmax(I[100:]))
            Ypk[i] = I[j]; Ytp[i] = t[j]
    except Exception:
        pass
ok = np.isfinite(Ypk) & np.isfinite(Ytp)
print(f"   valid {ok.sum()}/{NL}")

def prcc(Xm, y):
    n, p = Xm.shape
    Xr = np.apply_along_axis(rankdata, 0, Xm); yr = rankdata(y)
    v = np.zeros(p); pv = np.ones(p)
    for j in range(p):
        A = np.column_stack([np.delete(Xr, j, 1), np.ones(n)])
        rx = Xr[:, j] - A @ np.linalg.lstsq(A, Xr[:, j], rcond=None)[0]
        ry = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
        if rx.std() > 1e-12 and ry.std() > 1e-12:
            c, p_ = spearmanr(rx, ry)
            v[j] = c if np.isfinite(c) else 0
            pv[j] = p_ if np.isfinite(p_) else 1
    return v, pv

figF, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, (yv, ttl) in enumerate([(Ypk[ok], '(a) PRCC: peak infections $I_{\\max}$'),
                               (Ytp[ok], '(b) PRCC: time to peak $t^{*}$')]):
    v, pv = prcc(X[ok], yv)
    cols = ['#e74c3c' if z > 0 else '#2980b9' for z in v]
    ax[k].barh(range(NP), v, height=.66, color=cols, edgecolor='white',
               alpha=.88)
    ax[k].set_yticks(range(NP)); ax[k].set_yticklabels(names, fontsize=12)
    ax[k].set_xlabel('PRCC'); ax[k].set_title(ttl, fontweight='bold')
    ax[k].axvline(0, color='k', lw=.8)
    for x in (0.5, -0.5):
        ax[k].axvline(x, color='gray', ls='--', alpha=.4)
    ax[k].set_xlim(-1.1, 1.1)
    for i in range(NP):
        st = '***' if pv[i] < .001 else ('**' if pv[i] < .01
                                         else ('*' if pv[i] < .05 else ''))
        inside = abs(v[i]) > 0.62
        if inside:                       # keep label inside long bars
            xx = v[i] + (-0.03 if v[i] > 0 else 0.03)
            ha = 'right' if v[i] > 0 else 'left'
            col = 'white'
        else:
            xx = v[i] + (0.04 if v[i] >= 0 else -0.04)
            ha = 'left' if v[i] >= 0 else 'right'
            col = 'black'
        ax[k].text(xx, i, f'{v[i]:.2f}{st}', va='center', ha=ha,
                   fontsize=8.5, fontweight='bold', color=col)
    print(f"   {ttl}: {[f'{z:.2f}' for z in v]}")
figF.suptitle(f'Global sensitivity analysis (LHS/PRCC, $n={ok.sum()}$)',
              fontsize=13, fontweight='bold')
figF.tight_layout(); figF.savefig('FIGF_prcc.pdf'); plt.close(figF)

print("\nDONE: FIGA..FIGF")
