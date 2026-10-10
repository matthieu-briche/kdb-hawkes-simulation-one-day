/ =====================================================================
/ Maximum-likelihood estimation (MLE) of the Hawkes parameters
/ (mu, alpha, beta) from the quotes stored in the "hdbq" HDB,
/ on every (day, symbol) pair: 10 days x 5 symbols = 50 fits.
/ Usage: q mle.q   (after q hawkes_quotes.q, which creates hdbq)
/ ---
/ Exponential Hawkes log-likelihood on [0,T]:
/   log L = sum_i log(mu + alpha*A_i) - mu*T - (alpha/beta)*sum_i (1-exp(-beta(T-t_i)))
/   A_1 = 0 ;  A_i = exp(-beta(t_i - t_{i-1})) * (1 + A_{i-1})     (O(n) recursion)
/ Optimisation: Nelder-Mead on x = (log mu ; logit n ; log beta), n = alpha/beta
/   -> mu>0, beta>0 and 0<n<1 (stationarity) hold without constraints.
/ Standard errors: inverse of the numerical Hessian of -log L (observed information).
/ =====================================================================

\S 42                                                            / seed for the tie jitter only
\l lib/hawkes.q
-1"Loading mle.q ...";

/ ---------- likelihood ----------
logL:{[t;T;mu;al;be]
  e:exp neg[be]*1_deltas t;
  A:0f,{y*1+x}\[0f;e];
  sum[log mu+al*A]-(mu*T)+(al%be)*sum 1-exp neg[be]*T-t}

toNat:{[x] be:exp x 2; n:1%1+exp neg x 1; (exp x 0;n*be;be)}   / x -> (mu;alpha;beta)
nllX:{[t;T;x] neg logL[t;T] . toNat x}                          / objective to minimise
nllP:{[t;T;p] neg logL[t;T] . p}                                / same, natural parameters

/ ---------- Nelder-Mead (reflection / expansion / contraction / shrink) ----------
/ RELATIVE stopping rule: (worst - best) < tol * (1 + |best|)
nelderMead:{[f;x0;tol;maxit]
  n:count x0;
  S:(enlist x0),x0+/:0.5*{x=/:x}til n;                          / initial simplex
  F:f each S;
  i:iasc F; S:S i; F:F i; it:0;
  while[(it<maxit)&(tol*1+abs first F)<(last F)-first F;
    c:avg -1_S;                                                 / centroid without the worst point
    xr:c+c-last S; fr:f xr;
    $[fr<first F;
        [xe:c+2*xr-c; fe:f xe; $[fe<fr; [S[n]:xe; F[n]:fe]; [S[n]:xr; F[n]:fr]]];
      fr<F n-1;
        [S[n]:xr; F[n]:fr];
      [xc:$[fr<last F; c+0.5*xr-c; c+0.5*(last S)-c]; fc:f xc;
       $[fc<fr&last F;
          [S[n]:xc; F[n]:fc];
          [S:first[S]+/:0.5*S-\:first S; F:f each S]]]];
    i:iasc F; S:S i; F:F i; it+:1];
  `x`f`iter!(first S;first F;it)}

/ ---------- numerical Hessian (central differences) ----------
hess:{[f;p]
  n:count p;
  E:(1e-4*abs p)*{x=/:x}til n;
  g:{[f;p;E;i;j]
    a:f p+E[i]+E j; b:f p+E[i]-E j; c:f p+E[j]-E i; d:f p-E[i]+E j;
    (a+d-b+c)%4*E[i;i]*E[j;j]};
  {[gg;k;i] gg[i] each k}[g[f;p;E];til n] each til n}

/ ---------- estimation on every (day, symbol) of the HDB ----------
system"l hdbq";
-1"== Hawkes MLE on the HDB: ",string[count date]," days x ",string[count par]," symbols ==";

fit:{[d;s]
  t0:.z.p;
  t:1e-9*"j"$(exec time from select time from quote where date=d,sym=s)-sessOpen;   / seconds since the open (select then exec: KDB-X does not support a computed exec on a partitioned table)
  N:count t;
  x0:(log 0.5*N%T;0f;log 10f);                                   / deliberately far start
  r:nelderMead[nllX[t;T];x0;1e-9;1000];
  p:toNat r`x;
  se:sqrt {x[y;y]}[inv hess[nllP[t;T];p]] each til 3;
  lr:2*neg[r`f]-(N*log N%T)-N;                                   / LR vs homogeneous Poisson
  ks:ksStat comp[t] . p;                                         / residuals on the 1 ms times
  tj:asc T&0|t+0.001*-0.5+N?1f;                                  / each event spread uniformly within its ms
  ksj:ksStat comp[tj] . p;                                       / residuals on the de-quantised times
  -1"  ",string[d]," ",string[s],": N=",string[N],", ",string[r`iter]," iterations, ",string[`second$.z.p-t0]," (hh:mm:ss)";
  tr:par s;
  `date`sym`N`iter`mu`muHat`muSE`alpha`alphaHat`alphaSE`beta`betaHat`betaSE`n`nHat`LRvsPoisson`KS`KSjit`KS5pct!
   (d;s;N;r`iter;tr`mu;p 0;se 0;tr`alpha;p 1;se 1;tr`beta;p 2;se 2;tr[`alpha]%tr`beta;p[1]%p 2;lr;ks;ksj;1.358%sqrt N)}

res:raze {[d] fit[d] each exec sym from par} each date
z:update zMu:(muHat-mu)%muSE, zAlpha:(alphaHat-alpha)%alphaSE, zBeta:(betaHat-beta)%betaSE from res

/ ---------- first day ----------
d0:first date
-1"\n== first day (",string[d0],") ==";
show delete date from select from z where date=d0;
-1"\nmarkdown:";
-1 mdTable select sym,N,iter,muHat,muSE,alphaHat,alphaSE,betaHat,betaSE,nHat,LRvsPoisson,KS,KSjit,KS5pct from z where date=d0;

/ ---------- all days: is the estimator unbiased, are the SEs right? ----------
/ if both hold, z = (estimate - truth)/SE ~ N(0,1): mean ~ 0, sd ~ 1, and the
/ 95% interval estimate +/- 1.96 SE contains the truth ~95% of the time.
cover:{avg 1.96>abs x}
-1"\n== all days: standardised errors z = (estimate - truth)/SE, expected mean ~0, sd ~1 ==";
summ:select days:count i,
  zMuMean:avg zMu, zMuSD:dev zMu, zAlphaMean:avg zAlpha, zAlphaSD:dev zAlpha, zBetaMean:avg zBeta, zBetaSD:dev zBeta,
  cov95Mu:cover zMu, cov95Alpha:cover zAlpha, cov95Beta:cover zBeta,
  KSpass:sum KS<KS5pct, KSjitPass:sum KSjit<KS5pct by sym from z
show summ;
pool:select fits:count i, zMuMean:avg zMu, zAlphaMean:avg zAlpha, zBetaMean:avg zBeta,
  cov95Mu:cover zMu, cov95Alpha:cover zAlpha, cov95Beta:cover zBeta,
  KSpass:sum KS<KS5pct, KSjitPass:sum KSjit<KS5pct from z
-1"pooled over all fits:";
show pool;
-1"\nmarkdown:";
-1 mdTable summ;
-1 mdTable pool;

/ ---------- verdict ----------
-1"\nLR vs Poisson, chi2(2) 5% value 5.99 (indicative only: under H0 beta is not identified, see README) -> self-excitation ",$[all res[`LRvsPoisson]>5.99;"significant in every fit";"NOT significant in some fits"];
-1"KS on residuals, estimated parameters: 1 ms times pass ",string[sum res[`KS]<res`KS5pct],"/",string[count res],", de-quantised times pass ",string[sum res[`KSjit]<res`KS5pct],"/",string count res;
