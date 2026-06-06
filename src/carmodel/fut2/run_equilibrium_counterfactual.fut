import "autotrade/trmodel"
import "autotrade/equilibrium"
import "lib/github.com/diku-dk/linalg/lu"
import "lib/github.com/diku-dk/linalg/dpsolve"

import "run_equilibrium_module"

module R = f64
module trm = trmodel R
module eqb = equilibrium R trm
module dps = mk_dpsolve_dense R

module ed_ded_default = {
  module R   = R
  module trm = trm
  let ed_ded_price_all = eqb.ed_ded_price_all_man
}

module r = mk_run ed_ded_default

--- Build mp with heterogeneous pnew[j]=pnew_base+pnew_slope*j,
--- u_0[tau,j]=xi_0+eta*j and constant u_a[tau,j]=-gam.
def build_mp_het [n_]
    (c:i64) (Ax:i64)
    (transcost:R.t) (acc0:R.t) (mum:[n_]R.t)
    (xi_0:R.t) (eta:R.t) (gam:R.t)
    (pnew:[c]R.t)
    : ?[ns][nd]. r.trm.mp [n_][c][Ax][ns][nd] =
  let [ns][nd] mp : r.trm.mp [n_][c][Ax][ns][nd] = r.trm.mk n_ c Ax
  let u_0 = tabulate_2d n_ c (\_ j -> R.(xi_0 + eta * i64 j))
  let u_a = tabulate_2d n_ c (\_ _ -> R.(i64 0 - gam))
  let mp = r.trm.set_newprices mp pnew
  let mp = r.trm.set_u_0 u_0 mp
  let mp = r.trm.set_u_a u_a mp
  let mp = r.trm.set_transcost transcost mp
  let mp = r.trm.set_acc_0 (replicate c acc0) mp
  let mp = r.trm.set_mum mum mp
  in mp

--- Per-tau stationary distribution and choice probabilities at the given prices p.
def diagnostics_at_p [n_][c][Ax][ns][nd]
    (mp: r.trm.mp[n_][c][Ax][ns][nd]) (sa_max:i64) (p: r.trm.prices[c][Ax])
    : ([n_][ns]R.t, [n_][ns][nd]R.t) =
  let solve_evs (tau:i64) : ([ns]R.t, [ns][nd]R.t) =
    let utils = r.trm.utility mp p tau
    let tr = r.trm.age_transition mp
    let ev_init = r.trm.ev0 mp
    let f = r.trm.bellmanJ mp utils tr
    let param = dps.default with sa_max = sa_max
    let {res=ev, jac=_, conv=_, iter_sa=_, iter_nk=_, rtrips=_, tol=_} =
      dps.poly f ev_init param (R.i64 0)
    let (_, v) = r.trm.bellman0 mp utils tr ev
    in (ev, v)
  let evs = #[sequential_outer] map solve_evs (iota n_)
  let qcc (tau:i64) =
    let (ev, v) = evs[tau]
    let tr = r.trm.age_transition mp
    let ccp : [ns][nd]R.t = r.trm.ccp_tau mp v ev
    let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
    let (delta, _, _) = r.trm.trade_transition mp ccp
    let ctp = r.trm.ctp_tau tr delta
    let q_tau = eqb.ergodic ctp
    in (q_tau, ccp)
  let pairs = #[sequential_outer] map qcc (iota n_)
  let (q_taus, ccps) = unzip pairs
  in (q_taus, ccps)

--- Solve equilibrium under heterogeneous mp; return p, q_tau, ccp, tw plus
--- convergence info.
-- ==
-- entry: solve_eqb_het
-- input { 2i64 21i64 25i64 20i64 0f64 -5f64 0.5f64 1e-6f64 50i64
--         [0.1f64,0.3f64] 6.0f64 0.025f64 0.5f64 200f64 10f64 }
entry solve_eqb_het
    (n:i64) (c:i64) (Ax:i64) (sa_max:i64)
    (transcost:f64) (acc0:f64) (damp:f64) (tol:f64) (max_iter:i64)
    (mum:[n]f64)
    (xi_0:f64) (eta:f64) (gam:f64)
    (pnew_base:f64) (pnew_slope:f64)
    : ([c][Ax]f64, [n][]f64, [n][][]f64, [n]f64, i64, bool) =
  let pnew = tabulate c (\j -> pnew_base + pnew_slope * f64.i64 j)
  let mp = build_mp_het c Ax transcost acc0 mum xi_0 eta gam pnew
  let p0 = r.eqb.spp_price_solve mp 100
  let (p, _ed, _max_dp, iter, conv, _, _, _) =
    r.newton mp p0 sa_max damp tol max_iter
  let (q_taus, ccps) = diagnostics_at_p mp sa_max p
  in (p, q_taus, ccps, mp.tw, iter, conv)

--- Returns per-tau stationary distribution and choice probabilities with caller-supplied pnew
--- (under the assumption that mp.pnew matches p_in[:, 0]).
entry eval_at_p
    (n:i64) (c:i64) (Ax:i64) (sa_max:i64)
    (transcost:f64) (acc0:f64)
    (mum:[n]f64)
    (xi_0:f64) (eta:f64) (gam:f64)
    (pnew:[c]f64) (p_in:[c][Ax]f64)
    : ([n][]f64, [n][][]f64, [n]f64) =
  let mp = build_mp_het c Ax transcost acc0 mum xi_0 eta gam pnew
  let (q_taus, ccps) = diagnostics_at_p mp sa_max p_in
  in (q_taus, ccps, mp.tw)
