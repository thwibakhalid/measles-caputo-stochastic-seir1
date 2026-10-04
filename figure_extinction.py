#!/usr/bin/env python3
"""Figure E with the rigorously derived threshold R0^s = R0 * mu/(mu+nu+sigma^2/2)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fracseir import frac_sde, LAM, MU, DEL, GAM, ALP, N0, R0, NU_C

plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans',
                     'axes.labelsize': 12, 'axes.titlesize': 11,
                     'figure.dpi': 300, 'savefig.dpi': 300,
                     'savefig.bbox': 'tight', 'axes.grid': True,
                     'grid.alpha': 0.3})

QREF, NU_E = 0.95, 0.083
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))

# ---- (a) simulated fade-out probability -------------------------------
S_dfe = LAM / (MU + NU_E)
y0e = np.array([S_dfe, 300.0, 150.0, N0 - S_dfe - 450.0])
sg = np.arange(0.02, 0.21, 0.02)
pe = []
for s in sg:
    _, Y = frac_sde(QREF, NU_E, 300, 0.5, s, 150, seed=int(s * 1e5), y0=y0e)
    pe.append(float((Y[-60:, :, 2].max(axis=0) < 20).mean()))
    print(f"  sigma={s:.2f}  P(ext)={pe[-1]:.2f}")

ax[0].plot(sg, pe, 'o-', color='#8e44ad', lw=2.4, ms=7, mfc='white', mew=2)
ax[0].fill_between(sg, pe, alpha=.12, color='purple')
# theoretical threshold crossing
s_star = np.sqrt(2 * (R0 * MU - MU - NU_E))
ax[0].axvline(s_star, color='green', ls=':', lw=2.2,
              label=f'$\\mathcal{{R}}^s_0=1$ at $\\sigma^*={s_star:.3f}$')
ax[0].axhline(0.5, color='gray', ls='--', alpha=.5, label='50%')
ax[0].set_ylim(-0.04, 1.04)
ax[0].set_xlabel('Noise intensity $\\sigma$')
ax[0].set_ylabel('Extinction probability')
ax[0].set_title(f'(a) Fade-out probability ($\\nu={NU_E}$, $q={QREF}$)')
ax[0].legend(fontsize=8.5, framealpha=0.9)

# ---- (b) corrected vs superseded threshold ----------------------------
s = np.linspace(0.001, 0.40, 400)
R_old = R0 * np.exp(-0.5 * 3 * s**2)                    # superseded
R_new = R0 * MU / (MU + NU_E + 0.5 * s**2)              # derived here
R_nu = R0 * MU / (MU + NU_E)

ax[1].axhline(R0, color='k', ls='--', lw=1.0, alpha=.45,
              label=f'$\\mathcal{{R}}_0={R0:.1f}$ (no vaccination)')
ax[1].plot(s, R_old, 'r-', lw=2.3, alpha=.8,
           label='Superseded: $\\mathcal{R}_0e^{-\\sigma^2/2}$')
ax[1].plot(s, R_new, 'b-', lw=2.6,
           label='Derived: $\\mathcal{R}_0\\mu/(\\mu+\\nu+\\sigma_S^2/2)$')
ax[1].axhline(1.0, color='green', ls=':', lw=2.2, label='Threshold $=1$')
ax[1].axvline(s_star, color='green', ls=':', lw=1.2, alpha=.6)
ax[1].annotate(f'$\\sigma^*={s_star:.3f}$', xy=(s_star, 1.0),
               xytext=(s_star + 0.05, 2.2), fontsize=9, color='green',
               arrowprops=dict(arrowstyle='->', color='green', lw=1.4))
ax[1].set_yscale('log'); ax[1].set_ylim(0.4, 40); ax[1].set_xlim(0, 0.40)
ax[1].set_xlabel('Noise intensity $\\sigma$')
ax[1].set_ylabel('Reproduction number')
ax[1].set_title(f'(b) Stochastic thresholds ($\\nu={NU_E}$)')
ax[1].legend(fontsize=7.8, loc='lower left', framealpha=0.9)

fig.tight_layout(); fig.savefig('FIGE_extinction.pdf'); plt.close(fig)
print(f"sigma* = {s_star:.4f}   R_nu = {R_nu:.4f}")
print("saved FIGE_extinction.pdf")
