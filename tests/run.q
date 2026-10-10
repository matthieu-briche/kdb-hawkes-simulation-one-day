/ =====================================================================
/ Unit and statistical tests for lib/hawkes.q and lib/mle.q.
/ Usage (from the repository root):  q tests/run.q
/ Exit code 0 if every test passes, 1 otherwise (usable in CI).
/ No HDB is needed: everything is simulated here, with a fixed seed.
/ ---
/ Two kinds of tests:
/  - exact: compared with a brute-force O(n^2) version or a closed form,
/    to a numerical tolerance;
/  - statistical: a random quantity against its theoretical value, with a
/    tolerance of about 5 standard deviations (or twice a 5% critical value),
/    so that a pass does not depend on the seed and a failure means a bug.
/ =====================================================================

\S 42
\l lib/hawkes.q
\l lib/mle.q

/ ---------- minimal harness ----------
.t.pass:0; .t.fail:();
/ test[name;f]: f is a niladic function returning a boolean; an error counts as a failure
test:{[name;f]
  r:@[f;::;{(`err;x)}];
  ok:$[-1h=type r; r; 0b];
  $[ok; [.t.pass+:1; -1"  ok    ",name];
        [.t.fail,:enlist name; -1"  FAIL  ",name,$[`err~first r; "  (error: ",last[r],")"; ""]]];}
near:{[a;b;tol] all (raze/)tol>=abs a-b}                  / works on atoms, vectors, matrices

/ brute-force O(n^2) references
lamB:{[t;mu;al;be;x] mu+al*sum exp neg be*x-t where t<x}                   / intensity just before x
LamB:{[t;mu;al;be;x] (mu*x)+(al%be)*sum 1-exp neg be*x-t where t<x}        / compensator at x
logLB:{[t;T;mu;al;be] sum[log lamB[t;mu;al;be] each t]-LamB[t;mu;al;be;T]}

/ small sample without ties, reused by the exact tests
ts:asc 300?20f

-1"== lib/hawkes.q: simulator ==";
test["hawkesBranch: sorted times in [0,T)";
  {t:hawkesBranch[1f;30f;50f;1000f]; (t~asc t)&(all t>=0)&all t<1000f}];
test["hawkesBranch: alpha>=beta is rejected";
  {r:@[hawkesBranch[1f;60f;50f];1000f;{x}]; (10h=type r)&r like "non-stationary*"}];
test["hawkesBranch: alpha=0 gives a Poisson count, mean mu*T (+/- 5 sd)";
  {N:count hawkesBranch[2f;0f;50f;10000f]; 5>abs(N-20000)%sqrt 20000}];
test["hawkesBranch: mean count mu*T/(1-n) (+/- 5 sd, sd = sqrt(mu*T/(1-n)^3))";
  {N:count hawkesBranch[1f;25f;50f;20000f]; 5>abs(N-40000)%sqrt 20000%0.5 xexp 3}];
test["hawkesBranch: compensator residuals ~ Exp(1): mean, variance, KS, lag-1 acf";
  {t:hawkesBranch[1f;40f;50f;5000f]; x:comp[t;1f;40f;50f]; N:count x; se:1%sqrt N;
   (5>abs(avg[x]-1)%se)&(5>abs(var[x]-1)%sqrt 8%N)&(ksStat[x]<2*1.358*se)&5>abs[acf[x;1]]%se}];

-1"== lib/hawkes.q: residual tests ==";
test["comp: equals the brute-force compensator increments";
  {near[comp[ts;1f;30f;50f]; deltas LamB[ts;1f;30f;50f] each ts; 1e-9]}];
test["comp: alpha=0 reduces to mu * inter-event times";
  {near[comp[ts;2f;0f;50f]; 2*deltas ts; 1e-12]}];
test["comp: a tie (equal timestamps) gives an increment of exactly 0";
  {0f=last comp[ts,last ts;1f;30f;50f]}];
test["ksStat: one point at the Exp(1) median gives 0.5";
  {near[ksStat enlist log 2; 0.5; 1e-12]}];
test["ksStat: a mass p at 0 forces KS >= p (the 1 ms tie artefact)";
  {x:(1000#0f),neg log 1-9000?1f; 0.0999<ksStat x}];
test["ksStat: Exp(1) sample below twice the 5% critical value";
  {N:20000; ksStat[neg log 1-N?1f]<2*1.358%sqrt N}];
test["acf: lag 1 of an alternating series is -(n-1)/n";
  {near[acf[100#1 -1f;1]; -0.99; 1e-12]}];
test["acf: lag 0 is 1";
  {near[acf[ts;0]; 1f; 1e-12]}];

-1"== lib/hawkes.q: quotes ==";
test["mkQuote: schema, ms timestamps, sorted, inside the session";
  {qt:mkQuote`GOOG; tm:exec time from qt;
   (`time`sym`bid`ask`bsize`asize~cols qt)&("nsffjj"~exec t from meta qt)&
   (tm~asc tm)&(all 0=("j"$tm) mod 1000000)&(all tm>=sessOpen)&all tm<=sessOpen+"n"$1e9*T}];
test["mkQuote: spread of 1 to 3 cents, sizes are round lots of 100 to 1000";
  {qt:mkQuote`MSFT; s:"j"$100*qt[`ask]-qt`bid;
   (all s within 1 3)&(all qt[`bsize] in 100*1+til 10)&all qt[`asize] in 100*1+til 10}];
test["mdTable: header, separator, one line per row";
  {m:mdTable ([]a:1 2;b:`x`y); (4=count m)&(m[0]~"| a | b |")&(m[1]~"|---|---|")&m[2]~"| 1 | x |"}];

-1"== lib/mle.q: likelihood ==";
test["logL: recursion equals the brute-force O(n^2) log-likelihood";
  {near[logL[ts;20f;1f;30f;50f]; logLB[ts;20f;1f;30f;50f]; 1e-8]}];
test["logL: alpha=0 is the Poisson log-likelihood N log mu - mu T";
  {near[logL[ts;20f;2f;0f;50f]; (count[ts]*log 2)-2*20; 1e-9]}];
test["toNat: (log mu; logit n; log beta) maps back to (mu; alpha; beta)";
  {near[toNat (log 0.8;log 0.7%0.3;log 50); 0.8 35 50f; 1e-9]}];

-1"== lib/mle.q: optimiser and Hessian ==";
rosen:{a:x[1]-x[0]*x 0; b:1-x 0; (100*a*a)+b*b}
test["nelderMead: Rosenbrock from (-1.2;1) reaches (1;1)";
  {r:nelderMead[rosen;-1.2 1f;1e-12;5000]; near[r`x;1 1f;1e-4]&r[`iter]<5000}];
test["nelderMead: 3-d quadratic, minimum at (1;-2;3)";
  {r:nelderMead[{sum d*d:x-1 -2 3f};0 0 0f;1e-12;5000]; near[r`x;1 -2 3f;1e-4]}];
test["nelderMead: maxit is respected";
  {5=(nelderMead[rosen;-1.2 1f;1e-12;5])`iter}];
A:(4 1 0f;1 3 0.5;0 0.5 2f)
test["hess: exact on a quadratic form";
  {near[hess[{(0.5*x mmu A mmu x)+sum x};1 2 3f]; A; 1e-5]}];

-1"== end to end: simulate, fit, check ==";
/ one fit on continuous times, AAPL parameters, 5000 s (about 25,000 events)
e2e:{[T]
  t:hawkesBranch[1f;40f;50f;T]; N:count t;
  r:nelderMead[nllX[t;T];(log 0.5*N%T;0f;log 10f);1e-9;1000];
  p:toNat r`x;
  se:sqrt {x[y;y]}[inv hess[nllP[t;T];p]] each til 3;
  `p`se`z`iter!(p;se;(p-1 40 50f)%se;r`iter)}[5000f]
-1"  estimates ",(-3!e2e`p),"  SE ",(-3!e2e`se),"  z ",-3!e2e`z;
test["MLE: converges before maxit";             {e2e[`iter]<1000}];
test["MLE: positive, finite standard errors";   {all (e2e[`se]>0)&e2e[`se]<0w}];
test["MLE: |z| < 4 for mu, alpha, beta";        {all 4>abs e2e`z}];
test["MLE: LR against Poisson is large";
  {t:hawkesBranch[1f;40f;50f;2000f]; N:count t; lr:2*logL[t;2000f;1f;40f;50f]-(N*log N%2000)-N; lr>1000}];

-1"== 1 ms rounding and the residual KS test (README, section 3) ==";
/ TSLA parameters, 2000 s (about 24,000 events, ~11% ties after rounding)
rnd:{[T]
  t:hawkesBranch[1.2;45f;50f;T]; N:count t;
  tr:0.001*"j"$1000*t;                                    / rounded to the nearest ms, as in mkQuote
  tj:asc T&0|tr+0.001*-0.5+N?1f;                          / spread uniformly within the ms, as in mle.q
  `N`ties`ks`ksj`crit!(N;sum[0=1_deltas tr]%N;ksStat comp[tr;1.2;45f;50f];ksStat comp[tj;1.2;45f;50f];1.358%sqrt N)}[2000f]
-1"  N ",string[rnd`N],"  ties ",(.Q.f[4;rnd`ties]),"  KS rounded ",(.Q.f[4;rnd`ks]),"  KS de-quantised ",(.Q.f[4;rnd`ksj]),"  5% crit ",.Q.f[4;rnd`crit];
test["rounding: KS on rounded times is at least the share of ties";   {rnd[`ks]>=rnd[`ties]-1e-9}];
test["rounding: KS on rounded times rejects";                         {rnd[`ks]>rnd`crit}];
test["rounding: de-quantised times pass (below twice the critical value)"; {rnd[`ksj]<2*rnd`crit}];

/ ---------- summary ----------
-1"\n",string[.t.pass]," passed, ",string[count .t.fail]," failed";
if[count .t.fail; -1"failed: ",", " sv .t.fail];
exit $[count .t.fail;1;0]
