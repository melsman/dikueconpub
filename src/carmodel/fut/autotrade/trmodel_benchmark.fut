import "trmodel"

module trm = trmodel f64

-- ==
-- entry: bench_bellmanJ
-- input { 2i64 2i64 2i64 1i64 }
-- output { [[8.670580092671319f64, 0.22769411217320684f64, 8.670580092671319f64,
--           0.22769411217320684f64, 0.12769411068309072f64], [0.0147f64, 0.0798f64, 0.0147f64, 0.0399f64, 0.8010f64]] }
-- input { 2i64 2i64 2i64 2i64 }
-- output { [[13.43965002284233f64, 5.039207394202223f64, 13.43965002284233f64,
--           5.039207394202223f64, 4.939207392712107f64], [0.4705f64, 0.0008f64, 0.4705f64, 0.0004f64, 0.0077f64]]}

entry bench_bellmanJ (n:i64) (c:i64) (Ax:i64) (N:i64) : ?[ns].[2][ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp (replicate c 100.0f64)
  let p = trm.simple_prices mp 0.85
  let tau = 0
  let util : trm.utility [ns][nd] = trm.utility mp p tau
  let tr = trm.age_transition mp
  let (ev, dev) =
    loop (ev, _) = (trm.ev0 mp, (replicate ns (replicate ns 0.0f64))) for _i < N do
      let (ev', dev') = trm.bellmanJ mp util tr ev
      in (ev', dev')
  in [ev, dev[0]]

entry bench_bellmanJ_with_full (n:i64) (c:i64) (Ax:i64) (N:i64) : ?[ns].([ns][ns]f64) =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp (replicate c 100.0f64)
  let p = trm.simple_prices mp 0.85
  let tau = 0
  let util : trm.utility [ns][nd] = trm.utility mp p tau
  let tr = trm.age_transition mp
  let (_, dev) =
    loop (ev, _) = (trm.ev0 mp, (replicate ns (replicate ns 0.0f64))) for _i < N do
      let (ev', dev') = trm.bellmanJ mp util tr ev
      in (ev', dev')
  in dev