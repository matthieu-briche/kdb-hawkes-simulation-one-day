/ =====================================================================
/ Bid/ask quotes simulated on Hawkes event times (1 ms timestamps)
/ + statistical checks + "hdbq" HDB + 1 ms as-of grid
/ Standard kdb+tick HDB layout:
/   date (partition) | sym s (p) | time n | bid f | ask f | bsize j | asize j
/ Usage:  q hawkes_quotes.q
/ =====================================================================

\S 42
\l lib/hawkes.q

/ ---------- checks on one continuous-time simulation per symbol ----------
checkHawkes:{[s]
  r:par s;
  mu:r`mu; al:r`alpha; be:r`beta;
  tk:hawkesBranch[mu;al;be;T];
  inc:comp[tk;mu;al;be];                                   / must be i.i.d. Exp(1)
  w:10f;
  c:deltas 1+tk bin w*1+til floor T%w;                     / counts per 10 s window
  `sym`ticks`expected`mean`var`KS`KS5pct`acf1`fano`fanoTheo!(s;count tk;floor mu*T%1-al%be;avg inc;var inc;ksStat inc;1.358%sqrt count inc;acf[inc;1];(var c)%avg c;1%(1-al%be)xexp 2)}

-1"== Hawkes checks (one simulation per symbol) ==";
checks:checkHawkes each exec sym from par
show checks;
-1"expected: mean~1, var~1, KS<KS5pct, |acf1|<~0.02, fano>>1 (~fanoTheo, slightly below: see docs/fano-factor.md)";
-1"KS passes for ",string[sum checks[`KS]<checks`KS5pct],"/",string[count checks]," symbols (5 tests at 5%: an isolated failure is normal)";
-1"";
-1"markdown:";
-1 mdTable select sym,ticks,mean,var,KS,KS5pct,acf1,fano,fanoTheo from checks;

/ ---------- multi-day generation and HDB write ----------
system"rm -rf hdbq";                                       / start from a clean database
hdb:`:hdbq
hol:2024.01.01 2024.01.15                                  / NYSE holidays: New Year, MLK Day
cand:2024.01.01+til 30
dts:10#cand where (1<cand mod 7)&not cand in hol           / first 10 NYSE trading days

genDay:{[d]
  `quote set raze mkQuote each exec sym from par;
  if[not all value exec all 0<=deltas time by sym from quote; '"timestamps not non-decreasing"];
  .Q.dpft[hdb;d;`sym;`quote];
  -1 string[d]," : ",string[count quote]," quotes";}

genDay each dts;
delete quote from `.;

/ ---------- reload ----------
\l hdbq
-1"== HDB loaded ==";
show meta quote;                                           / expected: time n, sym a=p
if[not "n"=exec first t from meta quote where c=`time; '"time must be a timespan (n)"];
-1"events sharing their timestamp with the previous one (1 ms ties), per symbol, first day:";
show select ties:sum 0=1_deltas time by sym from quote where date=first dts;

/ ---------- quote value at EVERY ms (as-of): one minute of AAPL ----------
d0:first dts
g:([] sym:60000#`AAPL; time:sessOpen+"n"$1000000*til 60000)  / timespan, as in the HDB
q1:update sym:`AAPL from select time,bid,ask from quote where date=d0,sym=`AAPL
grid:aj[`sym`time;g;q1]
grid:update time:d0+time from grid                         / display: date+time -> timestamp
-1"== 1 ms grid (first rows after the first quote) ==";
show 5#select from grid where not null bid;
-1"grid rows: ",string count grid;
