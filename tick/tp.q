/ =====================================================================
/ Minimal tickerplant, self-contained (no KX kdb+tick scripts needed).
/ Usage (from the repository root):  q tick/tp.q [logdir] -p 5010
/   logdir: directory of the journal, default tplog
/ Same interface as kdb+tick for publishers and subscribers:
/   .u.upd[t;x]  x = list of columns (or one row), time column supplied by the feed
/   .u.sub[t;s]  t = table or ` (all), s = symbol list or ` (all);
/                returns (t;empty schema), or a list of them for t = `
/   .u.i, .u.L   message count and journal path, for replay with -11!(.u.i;.u.L)
/ Zero latency: every update is journalled, then published at once.
/ Not covered: end of day (.u.end), batching, chained tickerplants.
/ =====================================================================

\l sym.q

.u.t:tables`.
.u.w:.u.t!(count .u.t)#enlist ();                          / per table: list of (handle;syms)

/ ---------- journal: one (`upd;t;x) message per update ----------
.u.L:`$":",($[count .z.x;.z.x 0;"tplog"]),"/sym",string .z.D
if[()~key .u.L; .[.u.L;();:;()]];                          / create an empty journal (and its directory)
.u.i:{$[-7h=type x; x; first x]} -11!(-2;.u.L)            / messages already in it (valid prefix if truncated)
.u.l:hopen .u.L
-1"tickerplant: tables ",(", " sv string .u.t),", journal ",(1_string .u.L)," (",string[.u.i]," messages)";

/ ---------- subscriptions ----------
.u.del:{[t;h] .u.w[t]:.u.w[t] where not h=first each .u.w t;}
.u.sub:{[t;s]
  if[t~`; :.u.sub[;s] each .u.t];
  if[not t in .u.t; '"unknown table: ",string t];
  .u.del[t;.z.w];
  .u.w[t]:.u.w[t],enlist(.z.w;s);
  (t;value t)}
.z.pc:{[h] .u.del[;h] each .u.t;}

/ ---------- publish ----------
/ x: list of columns in schema order; s: ` (everything) or a symbol list
.u.sel:{[t;x;s] $[s~`; x; x[;where x[cols[t]?`sym] in s]]}
.u.pub:{[t;x]
  {[t;x;w] d:.u.sel[t;x;w 1]; if[count first d; neg[w 0](`upd;t;flip cols[t]!d)]}[t;x] each .u.w t;}

.u.upd:{[t;x]
  if[not t in .u.t; '"unknown table: ",string t];
  if[0>type first x; x:enlist each x];                     / one row -> one-row columns
  .u.l enlist(`upd;t;x); .u.i+:1;
  .u.pub[t;x];}
