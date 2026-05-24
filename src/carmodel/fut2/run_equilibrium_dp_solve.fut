
import "autotrade/trmodel"
import "autotrade/equilibrium"
import "lib/github.com/diku-dk/linalg/dpsolve"

module R = f64
module trm = trmodel R
module eqb = equilibrium R trm
module dps = mk_dpsolve_dense R

def solve_nk [n][c][Ax][ns][nd] (mp:trm.mp [n][c][Ax][ns][nd]) (p0: trm.prices[c][Ax]) (sa_max:i64) (tol:R.t) (max_iter:i64): ([c][Ax]f64, [c][Ax-1]f64, f64, i64, bool, [n]i64, [n]i64, [n]i64) =
  let endo_prices_from_flat (ps: [c*(Ax-1)]R.t) : trm.prices[c][Ax] =
      let p2d : [c][Ax-1]R.t = unflatten ps
      in tabulate_2d c Ax (\ct a -> if a==0 then p0[ct][0] else p2d[ct][a-1])
  let prices_to_flat (p: trm.prices[c][Ax]) : [c*(Ax-1)]R.t =
      flatten (tabulate_2d c (Ax-1) (\ct a -> p[ct][a+1]))

  let p_minus_ed (pf: [c*(Ax-1)]R.t)
      : ([c*(Ax-1)]R.t, [c*(Ax-1)][c*(Ax-1)]R.t) =
      let p = endo_prices_from_flat pf
      let (ed2d, ded4d, _, _, _) = eqb.ed_ded_price_all mp sa_max p
      let ed_flat  = flatten ed2d
      let ded_flat = map flatten (flatten ded4d)
      let p_minus_ed_flat = map2 (\price ed -> R.(price - ed)) pf ed_flat
      let JG = tabulate_2d (c*(Ax-1)) (c*(Ax-1)) (\i j ->
          let eye = R.bool (i==j)
          in R.(eye - ded_flat[i,j]))
      in (p_minus_ed_flat, JG)
  
  let param = dps.default with pi_max = max_iter
                            with pi_tol = tol

  let {res = p_flat, jac = _, conv, iter, tol = _} =
      dps.nk p_minus_ed (prices_to_flat p0) param

  let ed : [c][Ax-1]R.t = replicate c (replicate (Ax-1) (R.nan))
  let max_dp = R.nan
  let sa_iters_tot = replicate n (-1i64)
  let nk_iters_tot = replicate n (-1i64)
  let rtrips_tot   = replicate n (-1i64)

  let p = endo_prices_from_flat p_flat
  in (p, ed, max_dp, iter, conv, sa_iters_tot, nk_iters_tot, rtrips_tot)

--- Benchmarking format: n-c-Ax-acc0-transcost (acc0 is positive, negated inside)
-- ==
-- entry: bench_solve
-- input { 2i64 1i64 25i64 5f64 0f64 }
-- input { 2i64 5i64 25i64 5f64 0f64 }
-- input { 2i64 9i64 25i64 5f64 0f64 }
-- input { 2i64 13i64 25i64 5f64 0f64 }
-- input { 2i64 17i64 25i64 5f64 0f64 }
-- input { 2i64 21i64 25i64 5f64 0f64 }
entry bench_solve (n:i64) (c:i64) (Ax:i64) (acc0:R.t) (transcost:R.t) : f64 =
  let mum = tabulate n (\i -> if i == 0 then 0.1f64 else 0.3f64)
  let pnew = tabulate c (\i -> if i == 0 then 200.0f64 else 260.0f64)
  let u_0 = tabulate_2d n c (\_ j -> if j == 0 then 6.0f64 else 6.5f64)
  let u_a = tabulate_2d n c (\_ j -> if j == 0 then -0.5f64 else -0.475f64)
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp pnew
  let mp = trm.set_acc_0 (replicate c (R.neg acc0)) mp
  let mp = trm.set_transcost transcost mp
  let mp = trm.set_mum mum mp
  let mp = trm.set_u_0 u_0 mp
  let mp = trm.set_u_a u_a mp
  let p0 = eqb.spp_price_solve mp 100
  let (p, _, _, _, _, _, _, _) = solve_nk mp p0 20 (R.f64 1e-13) 20
  -- excluding pnew (age 0) to match MATLAB average
  let p_used = flatten (map (\row -> row[1:]) p)
  let tot_p = reduce (\x y -> R.(x + y)) (R.i64 0) p_used
  let count = R.i64 (c * (Ax - 1))
  in R.(tot_p/count)

--- Validation: same setup as bench_solve, but extends the output with
--- diagnostics that depend on the equilibrium solution. Returns:
---   p          equilibrium prices [c][Ax]
---   max|ed|    max absolute excess demand at p
---   stat_res   max over (tau, state) of |q_tau @ ctp_tau - q_tau|
---   norm_err   |sum(q_pop) - 1|, where q_pop = sum_tau tw_tau * q_tau
---   min_q      min entry of q_pop (catches negative / disconnected ctp)
---   iter       Newton iterations performed
---   conv       convergence flag
---   sa_tot, nk_tot, rtr_tot   per-household total SA / NK / round-trip counts

entry validate_solve (n:i64) (c:i64) (Ax:i64) (acc0:R.t) (transcost:R.t)
    : ([c][Ax]f64, f64, f64, f64, f64, i64, bool, [n]i64, [n]i64, [n]i64) =
  let mum = tabulate n (\i -> if i == 0 then 0.1f64 else 0.3f64)
  let pnew = tabulate c (\i -> if i == 0 then 200.0f64 else 260.0f64)
  let u_0 = tabulate_2d n c (\_ j -> if j == 0 then 6.0f64 else 6.5f64)
  let u_a = tabulate_2d n c (\_ j -> if j == 0 then -0.5f64 else -0.475f64)
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let sa_max = 20


  let mp = trm.set_newprices mp pnew
  let mp = trm.set_acc_0 (replicate c (R.neg acc0)) mp
  let mp = trm.set_transcost transcost mp
  let mp = trm.set_mum mum mp
  let mp = trm.set_u_0 u_0 mp
  let mp = trm.set_u_a u_a mp
  let p0 = eqb.spp_price_solve mp 100

  let (p, _, _, iter, conv, sa_iters_tot, nk_iters_tot, rtrips_tot) =
    solve_nk mp p0 sa_max (R.f64 1e-13) 20

  let (ed_real, _, _, _, _) = eqb.ed_ded_price_all mp sa_max p
  let max_abs_ed = reduce R.max (R.i64 0) (map R.abs (flatten ed_real))

  -- Per-tau ergodic check at the converged prices.
  let tr = trm.age_transition mp
  let ev0 = trm.ev0 mp
  let param = eqb.dps.default
  let (qs, stat_resids) =
    unzip (#[sequential_outer] map (\tau ->
      let utils = trm.utility mp p tau
      let f = trm.bellmanJ mp utils tr
      let {res=ev, jac=_, conv=_, iter_sa=_, iter_nk=_, rtrips=_, tol=_} =
        eqb.dps.poly f ev0 param (R.i64 0)
      let ctp = eqb.ctp_from_utils mp tr utils ev
      let q = eqb.ergodic ctp
      let q_at_ctp =
        map (\col -> reduce (R.+) (R.i64 0) (map2 (R.*) q col)) (transpose ctp)
      let res = reduce R.max (R.i64 0)
                  (map R.abs (map2 (R.-) q_at_ctp q))
      in (q, res)) (iota n))
  let stat_res = reduce R.max (R.i64 0) stat_resids
  
  -- Population ownership distribution: q_pop = sum_tau tw_tau * q_tau.
  let q_pop =
    reduce (map2 (R.+)) (replicate ns (R.i64 0))
           (map2 (\q tw -> map (\x -> R.(x * tw)) q) qs mp.tw)
  let sum_q = reduce (R.+) (R.i64 0) q_pop
  let norm_err = R.(abs (sum_q - i64 1))
  let min_q = reduce R.min R.highest q_pop

  in (p, max_abs_ed, stat_res, norm_err, min_q,
      iter, conv, sa_iters_tot, nk_iters_tot, rtrips_tot)