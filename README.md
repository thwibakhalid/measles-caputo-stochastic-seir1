# Caputo fractional-order stochastic SEIR model for measles

Code and data accompanying the manuscript

> T. A. Khalid, *A Caputo fractional-order stochastic SEIR model for measles
> with piecewise constant vaccination and seasonal forcing: Lyapunov
> stability, threshold dynamics and vaccination controllability*,
> submitted to the *Journal of the Nigerian Society of Physical Sciences*.

Everything needed to reproduce every figure and every number reported in the
paper is in this repository.

---

## What the model is

A Susceptible–Exposed–Infectious–Recovered model in which

* the transmission rate is seasonally forced,
  `beta(t) = beta_0 (1 + beta_1 sin(2 pi t / 52))`;
* vaccination is a **piecewise constant** schedule, constant on windows of
  fixed length and then reset;
* environmental variability enters as **multiplicative noise**;
* memory is represented by a **Caputo derivative of order `q`**.

The equation is solved in its stochastic Volterra form. The order is
restricted to `q` in `(1/2, 1]`, which is the range in which both the
well-posedness theory and the convergence of the numerical scheme are
established.

---

## Repository layout

```
src/
  fracseir.py                  core solvers (deterministic and stochastic)
  figures_main.py              Figures 2-5 and 7 of the paper
  figure_extinction.py         Figure 6 (fade-out probability, thresholds)
  figure_sudan.py              Figure 1 (Sudan surveillance)
  calibration_exploratory.py   exploratory fit reported in the Discussion
data/
  sudan_measles_monthly_2012_2022.csv
figures/
  FIGA..FIGG                   the figures as they appear in the paper
```

---

## Numerical method

Integration uses the Euler–Maruyama scheme for Caputo stochastic fractional
differential equations,

```
X_n = X_0 + (1/Gamma(q+1)) * sum_{j<n} w_{n,j} [ h^q f(t_j, X_j)
                                               + h^(q-1/2) G(X_j) xi_j ]
w_{n,j} = (n-j)^q - (n-1-j)^q,     xi_j ~ N(0,1)
```

The drift and noise increments carry the different scalings `h^q` and
`h^(q-1/2)`; at `q = 1` all weights equal one and the scheme reduces to the
classical Euler–Maruyama method. Strong convergence of order `q - 1/2` for
`q` in `(1/2, 1)` is due to Doan, Huong, Kloeden and Vu, *J. Comput. Appl.
Math.* **380** (2020) 112989. A short-memory truncation of length `M` is used
for efficiency; results were checked for sensitivity to `M` and to the step
size by halving each.

---

## Reproducing the figures

Requires Python 3.10 or later.

```bash
pip install -r requirements.txt

cd src
python figure_sudan.py        # Figure 1
python figures_main.py        # Figures 2, 3, 4, 5, 7
python figure_extinction.py   # Figure 6
```

Each script writes its PDF into the working directory. Random seeds are
fixed inside the scripts, so the figures are bit-for-bit reproducible.

Representative run times on a single core: `figure_sudan.py` a few seconds,
`figures_main.py` about twenty minutes, `figure_extinction.py` about ten
minutes.

---

## Data

`data/sudan_measles_monthly_2012_2022.csv` contains the monthly count of
reported measles cases for Sudan from January 2012 to December 2022, 130
months. It is extracted unchanged from the World Health Organization
Immunization Data portal, *Provisional monthly measles and rubella data:
distribution of cases by country, year and month*,
<https://immunizationdata.who.int>.

The series ends in December 2022, before the April 2023 conflict, after
which reporting was disrupted.

**On the use of these data.** They are used in the paper to motivate the
inclusion of a seasonal term and to bound its amplitude. They are *not* used
to calibrate the model and do not validate it: a variance decomposition
shows that the within-year seasonal component accounts for about 15% of the
total variance while differences between years account for roughly 50%,
which an autonomous seasonally forced model cannot reproduce. The
exploratory fit in `calibration_exploratory.py` is reported in the
Discussion precisely because it is inconclusive: the estimated memory order
drifts towards `q = 1` and the fractional model does not improve the
residual sum of squares on these data.

---

## Citation

If you use this code, please cite the paper above. Once the archived release
has a DOI, cite that as well.

## Licence

MIT, see `LICENSE`. The WHO data are redistributed under the terms of the
WHO data policy.
