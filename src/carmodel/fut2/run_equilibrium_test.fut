
import "autotrade/trmodel"
import "autotrade/equilibrium"
import "lib/github.com/diku-dk/linalg/lu"

module R = f64
module trm = trmodel R
module eqb = equilibrium R trm
module lu = mk_lu R

def mk_mp [n][c] (Ax:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (newprices:[c]R.t) =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let mp = trm.set_acc_0 (replicate c acc0) mp
  let mp = trm.set_transcost transcost mp
  let mp = trm.set_mum mum mp
  in mp

--- Compute all intermediates of ded_dprice_tau at given prices, for a given tau.
def intermediates_at_price [n][c][Ax] (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64) =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp (replicate c (R.i64 200))
  let mp = trm.set_acc_0 (replicate c acc0) mp
  let mp = trm.set_transcost transcost mp
  let mp = trm.set_mum mum mp
  let utils = trm.utility mp p tau
  let tr = trm.age_transition mp
  let ev0 = trm.ev0 mp
  let f = trm.bellmanJ mp utils tr
  let param = eqb.dps.default with sa_max = sa_max
  let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (R.i64 0)
  let (ev, v) = trm.bellman0 mp utils tr ev
  let ctp = eqb.ctp_from_utils mp tr utils ev
  let ccp_scrap = trm.ccp_scrap_tau mp p tau
  let ccp = trm.ccp_tau mp v ev
  let du = eqb.utility_dprice_man mp ccp_scrap tau
  let dbellman = eqb.dbellman_prices_man ccp du
  let dev = eqb.dev_fixed_point_dprice mp ctp dbellman
  let dctp = eqb.dctp_dprices_from_du_dev mp tr utils ev du dev
  let dccp = eqb.dccp_dprices_from_du_dev mp tr utils ev du dev
  let (erg, erg_dprice) = eqb.ergodic_with_dprice ctp dctp
  let dccp_scrap = eqb.ccp_scrap_dprice_man mp p tau
  let ded = eqb.ded_dprice ccp dccp erg erg_dprice ccp_scrap dccp_scrap
  in (dbellman, dev, dctp, dccp, erg_dprice, ded)

--- dbellman[c*(Ax-1)][ns] at given prices and tau
entry test_at_p_dbellman (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64)
    : ?[ns].[c*(Ax-1)][ns]R.t =
  let (dbellman, _, _, _, _, _) = intermediates_at_price sa_max transcost mum acc0 p tau
  in flatten dbellman

--- dev[c*(Ax-1)][ns] at given prices and tau
entry test_at_p_dev (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64)
    : ?[ns].[c*(Ax-1)][ns]R.t =
  let (_, dev, _, _, _, _) = intermediates_at_price sa_max transcost mum acc0 p tau
  in flatten dev

--- dctp[c*(Ax-1)][ns][ns] at given prices and tau
entry test_at_p_dctp (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64)
    : ?[ns].[c*(Ax-1)][ns][ns]R.t =
  let (_, _, dctp, _, _, _) = intermediates_at_price sa_max transcost mum acc0 p tau
  in flatten dctp

--- dccp[c*(Ax-1)][ns][nd] at given prices and tau
entry test_at_p_dccp (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64)
    : ?[ns][nd].[c*(Ax-1)][ns][nd]R.t =
  let (_, _, _, dccp, _, _) = intermediates_at_price sa_max transcost mum acc0 p tau
  in flatten dccp

--- erg_dprice[c*(Ax-1)][ns] at given prices and tau
entry test_at_p_erg_dprice (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64)
    : ?[ns].[c*(Ax-1)][ns]R.t =
  let (_, _, _, _, erg_dprice, _) = intermediates_at_price sa_max transcost mum acc0 p tau
  in flatten erg_dprice

--- ded[c][Ax-1][c][Ax-1] at given prices and tau (single type, before weighting)
entry test_at_p_ded (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t) (tau:i64)
    : [c][Ax-1][c][Ax-1]R.t =
  let (_, _, _, _, _, ded) = intermediates_at_price sa_max transcost mum acc0 p tau
  in ded

--- Full ed_ded_price_all at given prices (aggregated over all tau)
entry test_at_p_ed_ded (n:i64) (c:i64) (Ax:i64) (sa_max:i64) (transcost:R.t) (mum:[n]R.t) (acc0:R.t) (p:[c][Ax]R.t)
    : ([c][Ax-1]R.t, [c][Ax-1][c][Ax-1]R.t) =
  let mp = mk_mp Ax transcost mum acc0 (replicate c (R.i64 200))
  in eqb.ed_ded_price_all mp sa_max p
