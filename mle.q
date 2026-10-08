/ =====================================================================
/ Estimation par maximum de vraisemblance (MLE) des parametres Hawkes
/ (mu, alpha, beta) a partir des quotes stockees dans la HDB "hdbq"
/ Lancement : q mle.q   (apres q hawkes_quotes.q, qui cree hdbq)
/
/ Vraisemblance du Hawkes exponentiel sur [0,T] :
/   log L = sum_i log(mu + alpha*A_i) - mu*T - (alpha/beta)*sum_i (1-exp(-beta(T-t_i)))
/   A_1 = 0 ;  A_i = exp(-beta(t_i - t_{i-1})) * (1 + A_{i-1})     (recursion O(n))
/ Optimisation : Nelder-Mead sur x = (log mu ; logit n ; log beta), n = alpha/beta
/   -> mu>0, beta>0 et 0<n<1 (stationnarite) garantis sans contrainte.
/ Erreurs-types : inverse de la Hessienne numerique de -log L (information observee).
/ =====================================================================

-1"Chargement de mle.q ...";

/ ---------- vraisemblance ----------
logL:{[t;T;mu;al;be]
  e:exp neg[be]*1_deltas t;
  A:0f,{y*1+x}\[0f;e];
  sum[log mu+al*A]-(mu*T)+(al%be)*sum 1-exp neg[be]*T-t}

toNat:{[x] be:exp x 2; n:1%1+exp neg x 1; (exp x 0;n*be;be)}   / x -> (mu;alpha;beta)
nllX:{[t;T;x] neg logL[t;T] . toNat x}                          / objectif a minimiser
nllP:{[t;T;p] neg logL[t;T] . p}                                / idem, parametres naturels

/ ---------- Nelder-Mead (reflexion / expansion / contraction / retrecissement) ----------
/ critere d'arret RELATIF : (pire - meilleur) < tol * (1 + |meilleur|)
nelderMead:{[f;x0;tol;maxit]
  n:count x0;
  S:(enlist x0),x0+/:0.5*{x=/:x}til n;                          / simplexe initial
  F:f each S;
  i:iasc F; S:S i; F:F i; it:0;
  while[(it<maxit)&(tol*1+abs first F)<(last F)-first F;
    c:avg -1_S;                                                 / centroide hors pire point
    xr:c+c-last S; fr:f xr;
    $[fr<first F;
        [xe:c+2*xr-c; fe:f xe; $[fe<fr; [S[n]:xe; F[n]:fe]; [S[n]:xr; F[n]:fr]]];
      fr<F n-1;
        [S[n]:xr; F[n]:fr];
      [xc:$[fr<last F; c+0.5*xr-c; c+0.5*(last S)-c]; fc:f xc;
       $[fc<fr&last F;
          [S[n]:xc; F[n]:fc];
          [S:first[S]+/:0.5*S-\:first S; F:f each S]]]];
    i:iasc F; S:S i; F:F i; it+:1;
    if[0=it mod 50; -1"      it=",string[it],"  -logL=",string first F]];
  `x`f`iter!(first S;first F;it)}

/ ---------- Hessienne numerique (differences centrees) ----------
hess:{[f;p]
  n:count p;
  E:(1e-4*abs p)*{x=/:x}til n;
  g:{[f;p;E;i;j]
    a:f p+E[i]+E j; b:f p+E[i]-E j; c:f p+E[j]-E i; d:f p-E[i]+E j;
    (a+d-b+c)%4*E[i;i]*E[j;j]};
  {[gg;k;i] gg[i] each k}[g[f;p;E];til n] each til n}

/ ---------- adequation : residus du compensateur (doivent etre Exp(1)) ----------
comp:{[t;mu;al;be] d:deltas t; e:exp neg[be]*d; S:{1+x*y}\[0f;e];
  (mu*d)+(al%be)*(1-e)*0f,-1_S}
ksStat:{s:asc x; n:count s; F:1-exp neg s; i:1+til n;
  max[(i%n)-F] | max F-(i-1)%n}

/ ---------- vrais parametres (identiques a hawkes_quotes.q) ----------
par:([sym:`AAPL`MSFT`GOOG`AMZN`TSLA]
  mu:   1.0 0.8 0.5 0.6 1.2;
  alpha:40  35  30  35  45f;
  beta: 50  50  50  50  50f)
T:23400f

/ ---------- estimation sur un jour de la HDB ----------
system"l hdbq";
d0:first date
-1"== MLE Hawkes sur la HDB, date ",string[d0]," ==";

fit:{[s]
  t0:.z.p;
  t:1e-9*"j"$exec time-0D09:30:00 from quote where date=d0,sym=s;   / secondes depuis l'ouverture
  N:count t;
  -1"  ",string[s]," : N=",string[N]," evenements ...";
  x0:(log 0.5*N%T;0f;log 10f);                                   / depart volontairement loin
  r:nelderMead[nllX[t;T];x0;1e-9;1000];
  -1"    -> ",string[r`iter]," iterations, ",string[`second$.z.p-t0]," (hh:mm:ss)";
  p:toNat r`x;
  se:sqrt {x[y;y]}[inv hess[nllP[t;T];p]] each til 3;
  lr:2*neg[r`f]-(N*log N%T)-N;                                   / LR vs Poisson homogene
  ks:ksStat comp[t] . p;
  tr:par s;
  `sym`N`iter`mu`muHat`muSE`alpha`alphaHat`alphaSE`beta`betaHat`betaSE`n`nHat`LRvsPoisson`KS`KS5pct!
   (s;N;r`iter;tr`mu;p 0;se 0;tr`alpha;p 1;se 1;tr`beta;p 2;se 2;tr[`alpha]%tr`beta;p[1]%p 2;lr;ks;1.358%sqrt N)}

res:fit each exec sym from par
-1"";
show res;

/ ---------- resume ----------
z:update zMu:(muHat-mu)%muSE, zAlpha:(alphaHat-alpha)%alphaSE, zBeta:(betaHat-beta)%betaSE from res
-1"\n== ecarts standardises (vrai - estime)/SE : attendu |z| < ~2 ==";
show select sym,zMu,zAlpha,zBeta from z;
-1"LR vs Poisson : chi2(2) a 5% = 5.99 -> l'auto-excitation est ",$[all res[`LRvsPoisson]>5.99;"significative pour tous";"NON significative pour certains"];
-1"KS sur residus avec parametres estimes : ok pour ",string[sum res[`KS]<res`KS5pct],"/",string count res;
