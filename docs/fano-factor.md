# The Fano factor of a Hawkes process

This note derives the Fano factor used as a clustering check in [`hawkes_quotes.q`](../hawkes_quotes.q), and explains why the empirical value sits slightly **below** the asymptotic target $1/(1-n)^2$.

## Definition

The Fano factor over a window of length $T$ measures the dispersion of the event count $N(T)$:

$$
F(T) = \frac{\operatorname{Var}[N(T)]}{\mathbb{E}[N(T)]}
$$

For a Poisson process, $F(T) = 1$ for every $T$. A Hawkes process is **self-exciting**: each event temporarily raises the intensity, so events arrive in clusters. This produces **overdispersion**, i.e. $F > 1$.

## Asymptotic result (any kernel)

Consider a stationary linear Hawkes process with intensity

$$
\lambda(t) = \mu + \int_{-\infty}^{t} \varphi(t-s)\, dN(s),
$$

and **branching ratio** $n = \int_0^\infty \varphi(t)\,dt < 1$ (required for stationarity).

- The mean intensity is $\Lambda = \dfrac{\mu}{1-n}$.
- For large $T$, the variance grows as $\operatorname{Var}[N(T)] \sim \dfrac{\mu}{(1-n)^3}\,T$.

Hence:

$$
F(\infty) = \lim_{T\to\infty} F(T) = \frac{1}{(1-n)^2}
$$

Two remarks:

- It depends **only on $n$**, not on the shape of the kernel $\varphi$.
- It blows up as $n \to 1$, i.e. near criticality.

Inverting it gives a simple estimator of the branching ratio, $\hat n = 1 - 1/\sqrt{\hat F}$, with $\hat F$ computed over long windows.

## Exponential kernel: exact formula in $T$

With $\varphi(t) = \alpha e^{-\beta t}$ (the kernel used in this repo), $n = \alpha/\beta$. Let $\kappa = \beta - \alpha = \beta(1-n)$.

The (off-diagonal) covariance density is

$$
c(\tau) = \Lambda\,\frac{\alpha(2\beta-\alpha)}{2(\beta-\alpha)}\, e^{-\kappa|\tau|}.
$$

Integrating over the window, $\operatorname{Var}[N(T)] = \Lambda T + 2\int_0^T (T-u)\,c(u)\,du$, which gives:

$$
F(T) = 1 + \left(\frac{1}{(1-n)^2} - 1\right)\left(1 - \frac{1 - e^{-\kappa T}}{\kappa T}\right)
$$

How the Fano factor behaves with window size:

- For $T \ll 1/\kappa$, $F(T) \approx 1$: at short scales the process looks locally Poisson.
- For $T \gg 1/\kappa$, $F(T) \to 1/(1-n)^2$.
- The transition happens on the time scale $1/\kappa = 1/(\beta(1-n))$, which becomes very long near criticality.

## Application to this repo

`hawkes_quotes.q` counts events over 10 s windows and compares the result with `fanoTheo` $= 1/(1-n)^2$. With $\beta = 50$, the transition scale $1/\kappa$ is 50 ms to 200 ms, so 10 s windows are well into the asymptotic regime, but not exactly at the limit:

| sym | $n$ | $1/\kappa$ (ms) | $F(\infty)$ | $F(10\text{ s})$ | gap |
|---|---|---|---|---|---|
| AAPL | 0.80 | 100 | 25.00 | 24.76 | −1.0 % |
| MSFT | 0.70 | 67 | 11.11 | 11.04 | −0.6 % |
| GOOG | 0.60 | 50 | 6.25 | 6.22 | −0.4 % |
| AMZN | 0.70 | 67 | 11.11 | 11.04 | −0.6 % |
| TSLA | 0.90 | 200 | 100.00 | 98.02 | −2.0 % |

This finite-window term is the systematic reason the empirical Fano factor comes out "a bit less" than `fanoTheo`. On top of it, a 6h30 session gives about 2,340 windows, so the sampling error on the variance-to-mean ratio is roughly $\sqrt{2/2340} \approx 3\%$: a single run can land a few percent either side.

A quick numerical check (cluster simulation over $2 \times 10^5$ s, Python) agrees with the formula across window sizes:

| $n$ | window | empirical $F$ | theoretical $F(T)$ |
|---|---|---|---|
| 0.8 | 10 ms | 2.16 | 2.16 |
| 0.8 | 100 ms | 9.77 | 9.83 |
| 0.8 | 1 s | 22.2 | 22.6 |
| 0.8 | 10 s | 24.1 | 24.8 |
| 0.9 | 10 ms | 3.47 | 3.43 |
| 0.9 | 100 ms | 22.4 | 22.1 |
| 0.9 | 1 s | 82.6 | 80.3 |
| 0.9 | 10 s | 102 | 98.0 |

## In practice

Plotting $F(T)$ against $T$ on data is a good diagnostic for self-excitation. A curve that rises from 1 to a plateau suggests a subcritical Hawkes process: the plateau gives $n$, and the position of the transition gives the time scale of the kernel.

Caveat: a non-stationary baseline (e.g. intraday seasonality, a U-shaped $\mu(t)$) also inflates $F$. Remove it before reading the plateau as self-excitation.
