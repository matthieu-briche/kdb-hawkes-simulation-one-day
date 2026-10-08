/ =====================================================================
/ Quotes bid/ask simulees sur instants de Hawkes (horodatage a la ms)
/ + verifications statistiques + ecriture HDB "hdbq" + grille 1 ms
/ Schema HDB standard (kdb+tick) :
/   date (partition) | sym s (p) | time n | bid f | ask f | bsize j | asize j
/ Lancement :  q hawkes_quotes.q
/ =====================================================================

\S 42

/ ---------- outils ----------
normal:{sqrt[-2*log 1-x?1f]*cos 6.283185307179586*x?1f}

/ Hawkes exponentiel par representation en grappes (vectorise).
/ immigrants ~ Poisson(mu*T) (approx normale, negligeable pour mu*T grand) ;
/ chaque evenement a Poisson(alpha/beta) enfants, delais ~ Exp(beta).
/ Renvoie les temps tries en secondes dans [0,T).
hawkesBranch:{[mu;alpha;beta;T]
  if[alpha>=beta; '"non stationnaire"];
  n:alpha%beta;
  cdf:sums exp[neg n]*(n xexp til 40)%1f,prds 1f+til 39;   / cdf de Poisson(n)
  lam0:mu*T;
  n0:floor 0|0.5+lam0+sqrt[lam0]*first normal 1;           / immigrants
  cur:T*n0?1f;
  res:cur;
  while[count cur;
    k:1+cdf bin count[cur]?1f;                             / nb d'enfants par parent
    prn:cur where k;
    ch:prn+(neg log 1-count[prn]?1f)%beta;
    cur:ch where ch<T;
    res,:cur];
  asc res}

/ ---------- verifications ----------
ksStat:{s:asc x; n:count s; F:1-exp neg s; i:1+til n;
  max[(i%n)-F] | max F-(i-1)%n}
acf:{[x;k] y:x-avg x; (sum (k _ y)*(neg k) _ y)%sum y*y}

/ ---------- parametres par symbole (temps en secondes) ----------
/ mu: base (evts/s) ; alpha,beta: excitation (memoire ~ 1/beta = 20 ms)
/ rapport de branchement n=alpha/beta ; cadence moyenne = mu/(1-n)
par:([sym:`AAPL`MSFT`GOOG`AMZN`TSLA]
  mu:   1.0 0.8 0.5 0.6 1.2;
  alpha:40  35  30  35  45f;
  beta: 50  50  50  50  50f;
  p0:   150 300 140 130 250f;
  sig:  0.0001 0.0001 0.0001 0.0001 0.0002)

T:23400f                               / session 09:30-16:00

checkHawkes:{[s]
  r:par s;
  mu:r`mu; al:r`alpha; be:r`beta;
  tk:hawkesBranch[mu;al;be;T];
  d:deltas tk;
  e:exp neg[be]*d;
  S:{1+x*y}\[0f;e];
  inc:(mu*d)+(al%be)*(1-e)*0f,-1_S;                       / compensateur: doit etre Exp(1)
  w:10f;
  c:deltas 1+tk bin w*1+til floor T%w;                     / comptages par fenetre de 10 s
  `sym`ticks`expected`mean`var`KS`KS5pct`acf1`fano`fanoTheo!(s;count tk;floor mu*T%1-al%be;avg inc;var inc;ksStat inc;1.358%sqrt count inc;acf[inc;1];(var c)%avg c;1%(1-al%be)xexp 2)}

-1"== Verifications Hawkes (une simulation par symbole) ==";
checks:checkHawkes each exec sym from par
show checks;
-1"attendu: mean~1, var~1, KS<KS5pct, |acf1|<~0.02, fano>>1 (~fanoTheo, un peu moins)";
-1"KS ok pour tous: ",string[all checks[`KS]<checks`KS5pct]," (5 tests a 5% : un echec isole est normal)";

/ ---------- table quote pour un jour et un symbole ----------
/ time = heure dans la journee (timespan) ; la date est portee par la partition
mkQuote:{[s]
  r:par s;
  tk:hawkesBranch[r`mu;r`alpha;r`beta;T];
  n:count tk;
  m:floor 100*r[`p0]*exp sums r[`sig]*normal n;            / mid en cents (entier)
  spr:1+n?3;                                               / spread 1 a 3 cents
  ([] time:`timespan$0D09:30:00+"n"$1000000*"j"$1e3*tk;    / timespan, resolution 1 ms
      sym:n#s;
      bid:0.01*m;
      ask:0.01*m+spr;
      bsize:`long$100*1+n?10;
      asize:`long$100*1+n?10)}

/ ---------- generation multi-jours et ecriture HDB ----------
system"rm -rf hdbq";                                       / repart d'une base propre
hdb:`:hdbq
dts:2024.01.01+til 14
dts:dts where 1<dts mod 7                                  / jours ouvres

genDay:{[d]
  `quote set raze mkQuote each exec sym from par;
  if[not all value exec all 0<=deltas time by sym from quote; '"temps non croissants"];
  .Q.dpft[hdb;d;`sym;`quote];
  -1 string[d]," : ",string[count quote]," quotes";}

genDay each dts;
delete quote from `.;

/ ---------- relecture ----------
\l hdbq
-1"== HDB chargee ==";
show meta quote;                                           / attendu: time n, sym a=p
if[not "n"=exec first t from meta quote where c=`time; '"time doit etre un timespan (n)"];
-1"ex-aequo a la ms (quotes au meme instant), par symbole, 1er jour:";
show select ties:sum 0=1_deltas time by sym from quote where date=first dts;

/ ---------- valeur de cotation a CHAQUE ms (as-of) : 1 minute d'AAPL ----------
d0:first dts
g:([] sym:60000#`AAPL; time:0D09:30:00+"n"$1000000*til 60000)   / timespan, comme la HDB
q1:update sym:`AAPL from select time,bid,ask from quote where date=d0,sym=`AAPL
grid:aj[`sym`time;g;q1]
grid:update time:d0+time from grid                         / affichage: date+heure -> timestamp
-1"== grille 1 ms (premieres lignes apres le 1er quote) ==";
show 5#select from grid where not null bid;
-1"lignes de la grille: ",string count grid;
