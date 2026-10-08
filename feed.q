/ =====================================================================
/ Feed : publie une journee de quotes Hawkes simulees vers le tickerplant
/ Lancement : q feed.q [vitesse]     (vitesse = acceleration, defaut 60
/             -> 60 : une session de 6h30 rejouee en 6 min 30 s)
/ Le tickerplant doit tourner sur le port 5010.
/ =====================================================================

/ ---------- fonctions reprises de hawkes_quotes.q (sans ecriture HDB) ----------
normal:{sqrt[-2*log 1-x?1f]*cos 6.283185307179586*x?1f}

hawkesBranch:{[mu;alpha;beta;T]
  if[alpha>=beta; '"non stationnaire"];
  n:alpha%beta;
  cdf:sums exp[neg n]*(n xexp til 40)%1f,prds 1f+til 39;
  lam0:mu*T;
  n0:floor 0|0.5+lam0+sqrt[lam0]*first normal 1;
  cur:T*n0?1f;
  res:cur;
  while[count cur;
    k:1+cdf bin count[cur]?1f;
    prn:cur where k;
    ch:prn+(neg log 1-count[prn]?1f)%beta;
    cur:ch where ch<T;
    res,:cur];
  asc res}

par:([sym:`AAPL`MSFT`GOOG`AMZN`TSLA]
  mu:   1.0 0.8 0.5 0.6 1.2;
  alpha:40  35  30  35  45f;
  beta: 50  50  50  50  50f;
  p0:   150 300 140 130 250f;
  sig:  0.0001 0.0001 0.0001 0.0001 0.0002)

T:23400f

mkQuote:{[s]
  r:par s;
  tk:hawkesBranch[r`mu;r`alpha;r`beta;T];
  n:count tk;
  m:floor 100*r[`p0]*exp sums r[`sig]*normal n;
  spr:1+n?3;
  ([] time:`timespan$0D09:30:00+"n"$1000000*"j"$1e3*tk;
      sym:n#s;
      bid:0.01*m;
      ask:0.01*m+spr;
      bsize:`long$100*1+n?10;
      asize:`long$100*1+n?10)}

/ ---------- connexion au tickerplant ----------
h:@[hopen;`::5010;{-2"tickerplant injoignable sur le port 5010 : ",x; exit 1}]
spd:$[count .z.x; "F"$first .z.x; 60f]

day:`time xasc raze mkQuote each exec sym from par;      / tri global par heure (stable)
-1"quotes a publier : ",string[count day],"  (acceleration x",string[spd],")";

/ ---------- rejeu accelere : toutes les 100 ms, envoie les quotes deja "echus" ----------
pos:0
t0:.z.P
.z.ts:{
  tsim:0D09:30:00+"n"$spd*"j"$.z.P-t0;                   / heure simulee
  j:day[`time] bin tsim;                                 / dernier quote <= heure simulee
  if[j>=pos;
    neg[h](".u.upd";`quote;value flip day pos+til 1+j-pos);
    pos::j+1];
  if[pos>=count day; -1"journee publiee : ",string count day; system"t 0"]}
\t 100
