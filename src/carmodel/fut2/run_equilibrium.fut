
import "autotrade/trmodel"
import "autotrade/equilibrium"
import "lib/github.com/diku-dk/linalg/lu"

module mk_run (R:real) = {

  type t = R.t

  module trm = trmodel R
  module eqb = equilibrium R trm
  module lu = mk_lu R

  def newton_step [n][c][Ax][ns][nd]
      (mp: trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) (damp:t) (p: trm.prices[c][Ax])
      : (trm.prices[c][Ax], [c][Ax-1]t, t) =
    let (ed2d, ded4d) : ([c][Ax-1]t, [c][Ax-1][c][Ax-1]t) =
      eqb.ed_ded_price_all mp sa_max p

    let ed_flat  = flatten ed2d
    let ded_flat = map flatten (flatten ded4d)
    let blksz : i64 = 16
    let dp       = lu.ols blksz ded_flat ed_flat
    let dp2d : [c][Ax-1]t = unflatten dp

    -- Global alpha: largest step keeping all prices >= pscrap, capped at damp
    let alpha =
      tabulate_2d c (Ax-1) (\ct a ->
        let a' = a + 1
        in if R.(dp2d[ct][a] > i64 0)
           then R.((p[ct][a'] - mp.pscrap[ct]) / dp2d[ct][a])
           else R.highest)
      |> flatten |> reduce R.min damp

    let p_new : [c][Ax]t =
      tabulate_2d c Ax (\ct a ->
        if a == 0 then p[ct][0]
        else let a' = a - 1 in R.(p[ct][a] - alpha * dp2d[ct][a']))
    let max_dp = reduce R.max (R.i64 0) (map R.abs (map (\x -> R.(alpha * x)) dp))
    in (p_new, ed2d, max_dp)

  def newton [n][c][Ax][ns][nd] (mp: trm.mp[n][c][Ax][ns][nd]) (p0: trm.prices[c][Ax]) (sa_max:i64) (damp:t) (tol:t) (max_iter:i64)
      : (trm.prices[c][Ax], [c][Ax-1]t, t, i64, bool) =
    let ed0 : [c][Ax-1]t = replicate c (replicate (Ax-1) (R.i64 0))
    let (p, ed, max_dp, iter) =
      loop (p, ed, max_dp, iter) = (p0, ed0, R.highest, 0i64)
      while iter < max_iter && R.(max_dp > tol) do
        let (p', ed', max_dp') = newton_step mp sa_max damp p
        in (p', ed', max_dp', iter + 1)
    in (p, ed, max_dp, iter, R.(max_dp <= tol))
}

module R = f64
module r = mk_run R

entry solve (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (damp:R.t) (tol:R.t) (max_iter:i64): ([c][Ax]f64, [c][Ax-1]f64, f64, i64, bool) =
  let [ns][nd] mp : r.trm.mp [n][c][Ax][ns][nd] = r.trm.mk n c Ax
  let mp = r.trm.set_newprices mp (replicate c (R.i64 200))
  let mp = r.trm.set_acc_0 (replicate c acc0) mp
  let mp = r.trm.set_transcost transcost mp
  let mp = r.trm.set_mum mum mp
  let p0 = r.trm.simple_prices mp (R.f32 0.85)
  let (p, ed, max_dp, iter, conv) = r.newton mp p0 sa_max damp tol max_iter
  in (p, ed, max_dp, iter, conv)

--- Returns ed and flattened ded at the initial prices, before any Newton step.
entry test_ed_ded (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t)
    : ([][]R.t, [][][]R.t) =
  let [ns][nd] mp : r.trm.mp [n][c][Ax][ns][nd] = r.trm.mk n c Ax
  let mp = r.trm.set_newprices mp (replicate c (R.i64 200))
  let mp = r.trm.set_transcost transcost mp
  let p0 = r.trm.simple_prices mp (R.f32 0.85)
  let (ed, ded) = r.eqb.ed_ded_price_all mp sa_max p0
  in (ed, map flatten ded)

--- Returns dp (Newton step) and updated prices from a single step.
entry test_step (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t)
    : ([][]R.t, [][]R.t) =
  let [ns][nd] mp : r.trm.mp [n][c][Ax][ns][nd] = r.trm.mk n c Ax
  let mp = r.trm.set_newprices mp (replicate c (R.i64 200))
  let mp = r.trm.set_acc_0 (replicate c (R.f32 (-5))) mp
  let mp = r.trm.set_transcost transcost mp
  let mp = r.trm.set_mum mum mp
  let p0 = r.trm.simple_prices mp (R.f32 0.85)
  let (p1, _, _) = r.newton_step mp sa_max (R.f64 0.5) p0
  in (p0, p1)

--- Get demand and supply at equilibrium prices.
entry test_demand_supply_eq (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (tol:R.t) (max_iter:i64)
    : ([][]R.t, [][]R.t) =
  let [ns][nd] mp : r.trm.mp [n][c][Ax][ns][nd] = r.trm.mk n c Ax
  let mp = r.trm.set_newprices mp (replicate c (R.i64 200))
  let mp = r.trm.set_acc_0 (replicate c (R.f32 (-5))) mp
  let mp = r.trm.set_transcost transcost mp
  let mp = r.trm.set_mum mum mp
  let p0 = r.trm.simple_prices mp (R.f32 0.85)
  let (p_eq, _, _, _, _) = r.newton mp p0 sa_max (R.f64 0.5) tol max_iter
  in r.eqb.demand_supply_all mp sa_max p_eq

--- Demand and supply at initial prices
entry test_demand_supply_p0 (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t)
    : ([][]R.t, [][]R.t) =
  let [ns][nd] mp : r.trm.mp [n][c][Ax][ns][nd] = r.trm.mk n c Ax
  let mp = r.trm.set_newprices mp (replicate c (R.i64 200))
  let mp = r.trm.set_acc_0 (replicate c (R.f32 (-5))) mp
  let mp = r.trm.set_transcost transcost mp
  let mp = r.trm.set_mum mum mp
  let p0 = r.trm.simple_prices mp (R.f32 0.85)
  in r.eqb.demand_supply_all mp sa_max p0

-- ==
-- entry: solve
-- input { 2i64 1i64 25i64 20i64 0f64 [0.1f64, 0.3f64] 1e-6f64 20i64 }
-- output { ... }
