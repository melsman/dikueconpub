import "../trmodel"
import "../../lib/github.com/diku-dk/linalg/dpsolve"
import "../../lib/github.com/diku-dk/linalg/lu"

module trm = trmodel f64
module dps = mk_dpsolve_dense f64
module lu = mk_lu f64

local module bench_ols = {
  type t = f64
  type mat [n] = [n][n]t
  def blksz : i64 = 16
  def ols [n] (A:mat[n]) (b:[n]t) : [n]t = lu.ols blksz A b
  def eye n = tabulate_2d n n (\i j -> f64.bool (i == j))
  def sub a b = map2 (map2 (f64.-)) a b
}

-- ==
-- entry: bench_linear_solver
-- input { 1i64 30i64 25i64 }
-- output { [158.0146f64, 114020.0f64] }
entry bench_linear_solver (n:i64) (c:i64) (Ax:i64) : [2]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp (replicate c 100.0f64)
  let p = trm.simple_prices mp 0.85
  let u_0 = tabulate_2d n c
    (\i j -> f64.(i64 5 + i64 2 * (i64 i + i64 j) / (i64 n + i64 c)))
  let mp = trm.set_u_0 u_0 mp
  let tau = 0
  let util : trm.utility [ns][nd] = trm.utility mp p tau
  let tr = trm.age_transition mp
  let ev0 = trm.ev0 mp
  let (ev1, dV) = trm.bellmanJ mp util tr ev0
  let F = bench_ols.sub (bench_ols.eye ns) dV
  let x = map2 (f64.-) ev0 (bench_ols.ols F (map2 (f64.-) ev0 ev1))
  in [reduce f64.max (-f64.inf) x, reduce (+) 0.0 x]

-- ==
-- entry: bench_newton_single
-- input { 1i64 20i64 40i64 40i64 0.0f64 }
-- output { [178.7272f64, 40.0f64] }


---- n is superfluous, but keeping it to make it easier to create a version with multiple consumer types if I want to do that later
entry bench_newton_single (n:i64) (c:i64) (Ax:i64) (iter:i64) (tol:f64) : [2]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp (replicate c 100.0f64)
  let p = trm.simple_prices mp 0.85
  let u_0 = tabulate_2d n c
		(\i j -> f64.(i64 5 + i64 2 * (i64 i + i64 j) / (i64 n + i64 c)))
  let mp = trm.set_u_0 u_0 mp
  let param = dps.default with pi_max = iter
                            with pi_tol = tol
  let tau = 0
  let util : trm.utility [ns][nd] = trm.utility mp p tau
  let tr = trm.age_transition mp
  let f = trm.bellmanJ mp util tr
  let ev0 = trm.ev0 mp
  let {res=ev,jac=_,conv=_,iter=i,tol=_} = dps.nk f ev0 param
  in [reduce f64.max (0.0) ev, (f64.i64 i)]

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