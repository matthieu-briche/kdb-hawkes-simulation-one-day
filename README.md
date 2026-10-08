# Hawkes Quote Simulator — kdb+/q

End-to-end simulation of an intraday bid/ask quote stream driven by **self-exciting Hawkes processes**, built entirely in **kdb+/q**: simulation, statistical validation, partitioned historical database, as-of analytics, real-time replay through a tickerplant, and **maximum-likelihood recovery of the model parameters**.

> All data is simulated. This is a personal project to explore market-microstructure modelling and the kdb+ tick stack.

---

## Why Hawkes?

On real markets, events are not independent: a quote update tends to trigger further updates within milliseconds. A homogeneous Poisson process cannot reproduce this **clustering**. The Hawkes process can: each event temporarily raises the arrival intensity.

$$
\lambda(t) = \mu + \sum_{t_i < t} \alpha\, e^{-\beta (t - t_i)}
$$

| Parameter | Meaning | Value used |
|---|---|---|
| $\mu$ | baseline intensity (events/s) | 0.5 – 1.2 |
| $\alpha$ | jump in intensity after each event | 30 – 45 |
| $\beta$ | decay rate (memory $1/\beta$ = 20 ms) | 50 |
| $n = \alpha/\beta$ | branching ratio (must be < 1) | 0.6 – 0.9 |

Stationary mean rate: $\mu / (1-n)$, i.e. up to ~12 quotes/s on average for TSLA, with bursts much higher.

---

## Architecture

```
                       hawkes_quotes.q
   ┌──────────────┐    ┌──────────────────────────────┐      ┌──────────┐
   │  Simulation  │───►│ statistical checks            │─────►│  HDB     │──► mle.q
   │  (branching) │    │ (compensator, KS, ACF, Fano)  │      │  hdbq/   │    parameter
   └──────────────┘    └──────────────────────────────┘      │  date/   │    recovery
                                                              │  quote   │
                                                              └────┬─────┘
                                                                   └──► aj : bid/ask on a 1 ms grid
   feed.q                     tick.q (port 5010)            r.q
   ┌──────────────┐  .u.upd   ┌──────────────┐  publish   ┌──────────┐
   │ accelerated  │──────────►│  tickerplant │───────────►│   RDB    │
   │ replay (x60) │           │  schema sym.q│            │          │
   └──────────────┘           └──────────────┘            └──────────┘
```

| File | Role |
|---|---|
| [`hawkes_quotes.q`](hawkes_quotes.q) | simulation, statistical checks, writes a 10-day partitioned HDB, 1 ms as-of grid |
| [`mle.q`](mle.q) | maximum-likelihood estimation of $(\mu,\alpha,\beta)$ from the HDB, standard errors, goodness of fit |
| [`feed.q`](feed.q) | replays a simulated session into a tickerplant in accelerated time |
| [`sym.q`](sym.q) | kdb+tick schema for the `quote` table |

---

## 1. Simulation — cluster representation

Instead of the classic thinning algorithm (sequential, one event at a time), the process is simulated through its **branching (cluster) representation**, which vectorises naturally in q:

1. **Immigrants**: Poisson($\mu T$) events uniformly distributed on $[0, T)$.
2. **Offspring**: each event generates Poisson($n$) children, each delayed by an Exp($\beta$) time.
3. Repeat on each new generation until extinction (guaranteed since $n<1$).

Each generation is processed as a whole vector: no per-event loop.

Prices follow a geometric random walk on the event clock, in integer cents, with a 1–3 cent spread. Timestamps are stored as `timespan` at **1 ms resolution**, partitioned by date with the `p#` attribute on `sym` (standard kdb+tick layout, written with `.Q.dpft`). The random seed is fixed (`\S 42`) so every run is reproducible.

---

## 2. Statistical validation

A simulator is only useful if it is checked. `hawkes_quotes.q` runs four tests per symbol:

- **Time-rescaling theorem**: the compensator increments $\Lambda(t_i) - \Lambda(t_{i-1})$ must be i.i.d. Exp(1). Checked with the mean, the variance and a **Kolmogorov–Smirnov test**.
- **Independence**: lag-1 autocorrelation of these increments ≈ 0.
- **Clustering**: the Fano factor (variance/mean of counts over 10 s windows) must be ≫ 1 and approach the theoretical value $1/(1-n)^2$, versus 1 for a Poisson process.
- **Data integrity**: strictly non-decreasing timestamps per symbol before writing to disk, and type checks on reload.

---

## 3. Parameter recovery by maximum likelihood (`mle.q`)

The final check: can the true parameters be recovered **from the stored 1 ms data alone**?

The exponential kernel gives a log-likelihood computable in $O(n)$ through a recursion:

$$
\log L = \sum_i \log(\mu + \alpha A_i) - \mu T - \frac{\alpha}{\beta}\sum_i \left(1 - e^{-\beta (T - t_i)}\right),
\qquad A_i = e^{-\beta (t_i - t_{i-1})}(1 + A_{i-1})
$$

- **Optimiser**: a Nelder–Mead simplex written in q, run on $(\log\mu,\ \operatorname{logit} n,\ \log\beta)$ so that positivity and stationarity ($n<1$) hold without constraints.
- **Standard errors**: inverse of the numerical Hessian of $-\log L$ (observed Fisher information).
- **Hawkes vs Poisson**: likelihood-ratio test against a homogeneous Poisson process ($\chi^2_2$ at 5 % = 5.99).
- **Goodness of fit**: KS test on the compensator residuals computed with the *estimated* parameters.

### Results

Output of `q mle.q` on the first day of the HDB (seed 42):

| sym | N events | $\mu$ (true → est.) | $\alpha$ (true → est.) | $\beta$ (true → est.) | $n$ (true → est.) | LR vs Poisson |
|---|---|---|---|---|---|---|
| AAPL | | 1.0 → | 40 → | 50 → | 0.80 → | |
| MSFT | | 0.8 → | 35 → | 50 → | 0.70 → | |
| GOOG | | 0.5 → | 30 → | 50 → | 0.60 → | |
| AMZN | | 0.6 → | 35 → | 50 → | 0.70 → | |
| TSLA | | 1.2 → | 45 → | 50 → | 0.90 → | |

<!-- Fill in from the `res` and `z` tables printed by `q mle.q`. -->

---

## How to run

Requires [kdb+](https://kx.com/kdb-personal-edition-download/) (free personal / community edition).

```bash
# 1. simulate, validate, write the HDB (hdbq/) and build the 1 ms grid
q hawkes_quotes.q

# 2. recover the parameters from the HDB
q mle.q
```

Real-time replay with the standard [kdb+tick](https://github.com/KxSystems/kdb-tick) scripts (place `sym.q` in `tick/`):

```bash
q tick.q sym . -p 5010          # tickerplant
q tick/r.q :5010 -p 5011        # real-time database
q feed.q 60                     # replay one session at x60 (6h30 in 6 min 30 s)
```

---

## Possible extensions

- Real-time subscriber computing spread and estimated intensity on the fly
- Multivariate Hawkes (cross-excitation between symbols, or between bid and ask sides)
- Intraday seasonality: time-varying baseline $\mu(t)$ (U-shaped activity)
- Price impact: make price moves depend on order-flow clustering

## Tech

`kdb+/q` · `kdb+tick` · `.Q.dpft` partitioned HDB · `aj` as-of join · point processes · maximum likelihood · Nelder–Mead
