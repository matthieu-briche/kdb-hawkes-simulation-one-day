/ =====================================================================
/ Checks that the RDB holds exactly the HDB day replayed by feed.q.
/ Usage (from the repository root):  q tests/tick_check.q rdbPort date [timeout_s]
/ Polls the RDB until its row count reaches the HDB count, then compares
/ the two tables row by row. Exit code 0 on success, 1 otherwise.
/ Called by tests/tick_smoke.sh.
/ =====================================================================

rdbPort:.z.x 0
dt:"D"$.z.x 1
tmo:$[2<count .z.x; "J"$.z.x 2; 120]

system"l hdbq";
exp0:`time xasc update sym:value sym from select time,sym,bid,ask,bsize,asize from quote where date=dt
-1"tick_check: ",string[count exp0]," quotes expected for ",string dt;

h:@[hopen;`$"::",rdbPort;{-2"rdb unreachable on port ",rdbPort,": ",x; exit 1}]
t0:.z.P
while[(count[exp0]>n:h"count quote")&tmo>`long$(.z.P-t0)%1e9; system"sleep 1"];
-1"tick_check: rdb has ",string[n]," rows after ",string[`long$(.z.P-t0)%1e9]," s";

got:h"select time,sym,bid,ask,bsize,asize from quote"
strip:{flip {`#x} each flip x}                             / compare values, not attributes (g#, p#, s#)
ok:(count[exp0]=count got)&strip[exp0]~strip `time xasc got
-1"tick_check: ",$[ok;"PASS, the RDB matches the HDB day";"FAIL, the RDB differs from the HDB day"];
exit $[ok;0;1]
