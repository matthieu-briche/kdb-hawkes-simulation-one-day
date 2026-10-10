/ =====================================================================
/ Minimal real-time database for tick/tp.q (same subscribe/replay protocol as kdb+tick's r.q).
/ Usage (from the repository root):  q tick/rdb.q [host:port] -p 5011
/   host:port of the tickerplant, default :5010
/ On start: subscribes to every table and symbol, then replays the journal
/ up to the subscription point, so no update is lost or counted twice.
/ Not covered: end of day (no save to the HDB, the day stays in memory).
/ =====================================================================

upd:{[t;x] t insert x;}

tpAddr:`$":",$[count .z.x; .z.x 0; ":5010"]
tph:@[hopen;tpAddr;{-2"tickerplant unreachable at ",(1_string tpAddr),": ",x; exit 1}]

/ subscribe and read (.u.i;.u.L) in one synchronous call: updates published after it
/ arrive on the handle and are processed after the replay below
r:tph"(.u.sub[`;`];.u `i`L)"
(.[;();:;].) each r 0;                                     / empty tables with the published schema
if[not null last r 1; -11!r 1];                            / replay the first .u.i journal messages
-1"rdb: subscribed to ",(1_string tpAddr),", tables ",(", " sv string tables`.),", ",string[sum count each value each tables`.]," rows after replay";

/ ---------- a few live queries ----------
/ q)stats[]                               rows, last time, average spread per symbol
/ q)last1[`TSLA]                          latest quote
stats:{select rows:count i, last time, avgSpreadCents:100*avg ask-bid by sym from quote}
last1:{[s] select from quote where sym=s, i=last i}
