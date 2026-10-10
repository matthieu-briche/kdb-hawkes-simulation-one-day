/ =====================================================================
/ Estimation tools shared by mle.q and tests/run.q: exponential Hawkes
/ log-likelihood, reparametrisation, Nelder-Mead simplex, numerical Hessian.
/ Pure functions, no I/O, no random draws: loading this file has no side effect.
/ =====================================================================

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
