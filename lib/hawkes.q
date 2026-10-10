/ =====================================================================
/ Shared Hawkes model code: parameters, simulation, residual tests.
/ Loaded by hawkes_quotes.q, mle.q and feed.q (\l lib/hawkes.q). Single source
/ of truth for the true parameters, so the scripts cannot drift apart.
/ Loading this file draws no random numbers: the caller sets the seed.
/ =====================================================================

/ ---------- true parameters per symbol (time in seconds) ----------
/ mu: baseline (events/s) ; alpha, beta: excitation (memory 1/beta = 20 ms)
/ branching ratio n = alpha/beta ; stationary mean rate = mu/(1-n)
/ p0: initial mid price ; sig: per-event log-price volatility
par:([sym:`AAPL`MSFT`GOOG`AMZN`TSLA]
  mu:   1.0 0.8 0.5 0.6 1.2;
  alpha:40  35  30  35  45f;
  beta: 50  50  50  50  50f;
  p0:   150 300 140 130 250f;
  sig:  0.0001 0.0001 0.0001 0.0001 0.0002)

T:23400f                                                   / session length, 09:30-16:00
sessOpen:0D09:30:00

/ ---------- tools ----------
normal:{sqrt[-2*log 1-x?1f]*cos 6.283185307179586*x?1f}    / Box-Muller, x draws

/ Exponential Hawkes through the branching (cluster) representation, vectorised.
/ immigrants ~ Poisson(mu*T) (normal approximation, negligible for large mu*T);
/ each event has Poisson(alpha/beta) children, delays ~ Exp(beta).
/ Returns the sorted event times in seconds, in [0,T).
hawkesBranch:{[mu;alpha;beta;T]
  if[alpha>=beta; '"non-stationary: alpha must be < beta"];
  n:alpha%beta;
  cdf:sums exp[neg n]*(n xexp til 40)%1f,prds 1f+til 39;   / cdf of Poisson(n)
  lam0:mu*T;
  n0:floor 0|0.5+lam0+sqrt[lam0]*first normal 1;           / number of immigrants
  cur:T*n0?1f;
  res:cur;
  while[count cur;
    k:1+cdf bin count[cur]?1f;                             / children per parent
    prn:cur where k;
    ch:prn+(neg log 1-count[prn]?1f)%beta;
    cur:ch where ch<T;
    res,:cur];
  asc res}

/ ---------- residual tests ----------
/ Compensator increments Lambda(t_i)-Lambda(t_{i-1}): i.i.d. Exp(1) under the model
/ (time-rescaling theorem). O(n) through S_i = 1 + exp(-beta d_i) S_{i-1}.
comp:{[t;mu;al;be]
  d:deltas t;
  e:exp neg[be]*d;
  S:{1+x*y}\[0f;e];
  (mu*d)+(al%be)*(1-e)*0f,-1_S}

/ Kolmogorov-Smirnov statistic against Exp(1); 5% critical value ~ 1.358/sqrt n
ksStat:{s:asc x; n:count s; F:1-exp neg s; i:1+til n;
  max[(i%n)-F] | max F-(i-1)%n}

/ autocorrelation of x at lag k
acf:{[x;k] y:x-avg x; (sum (k _ y)*(neg k) _ y)%sum y*y}

/ ---------- one day of quotes for one symbol ----------
/ time = time of day (timespan, 1 ms resolution); the date is carried by the partition.
/ Integer-cent mid on a geometric random walk (event clock), spread 1 to 3 cents
/ placed around the mid: |(bid+ask)/2 - mid| <= 0.5 cent.
mkQuote:{[s]
  r:par s;
  tk:hawkesBranch[r`mu;r`alpha;r`beta;T];
  n:count tk;
  m:floor 100*r[`p0]*exp sums r[`sig]*normal n;            / mid in integer cents
  spr:1+n?3;                                               / spread in cents
  b:m-spr div 2;                                           / bid in integer cents
  ([] time:`timespan$sessOpen+"n"$1000000*"j"$1e3*tk;      / rounded to the nearest ms
      sym:n#s;
      bid:0.01*b;
      ask:0.01*b+spr;
      bsize:`long$100*1+n?10;
      asize:`long$100*1+n?10)}

/ ---------- markdown table (to paste results into the README) ----------
mdCell:{$[-9h=type x; .Q.f[3;x]; 10h=type x; x; string x]}
mdTable:{[t] t:0!t; c:cols t;
  ("| ",(" | " sv string c)," |";"|",raze (count c)#enlist"---|"),
  {"| ",(" | " sv mdCell each value x)," |"} each t}
