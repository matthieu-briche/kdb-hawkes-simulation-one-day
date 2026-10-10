/ =====================================================================
/ Feed: replays one day of the "hdbq" HDB into the tickerplant in
/ accelerated time, so the real-time stream is exactly a stored day.
/ Usage: q feed.q [speed] [date]
/   speed: acceleration, default 60 -> a 6h30 session replayed in 6 min 30 s
/   date : HDB partition to replay (YYYY.MM.DD), default the first day
/ Requires hdbq (q hawkes_quotes.q) and a tickerplant on port 5010.
/ =====================================================================

h:@[hopen;`::5010;{-2"tickerplant unreachable on port 5010: ",x; exit 1}]
spd:$[count .z.x; "F"$.z.x 0; 60f]
\l lib/hawkes.q

system"l hdbq";
dt:$[1<count .z.x; "D"$.z.x 1; first date]
if[not dt in date; -2"date not in the HDB: ",string dt; exit 1]

/ one day, all symbols, in global time order (xasc is stable); de-enumerate sym for IPC
day:`time xasc update sym:value sym from select time,sym,bid,ask,bsize,asize from quote where date=dt
-1"quotes to publish for ",string[dt],": ",string[count day],"  (speed x",string[spd],")";

/ ---------- accelerated replay: every 100 ms, send the quotes already "due" ----------
pos:0
t0:.z.P
.z.ts:{
  tsim:sessOpen+"n"$spd*"j"$.z.P-t0;                     / simulated clock
  j:day[`time] bin tsim;                                 / last quote <= simulated clock
  if[j>=pos;
    neg[h](".u.upd";`quote;value flip day pos+til 1+j-pos);
    pos::j+1];
  if[pos>=count day; -1"day published: ",string count day; system"t 0"]}
\t 100
