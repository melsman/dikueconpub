import "equilibrium"
import "trmodel"

module trm = trmodel f64
module eqb = equilibrium f64 trm

-- ==
-- entry: test_edf_deriv
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[[-0.023859654165286f64], [0.023198151962656f64]]], 
--             [[[0.023198151962656f64], [-0.023859654165286f64]]]]}

entry test_edf_deriv [c] (n:i64) (newprices:[c]f64) (Ax:i64) : [c][Ax-1][c][Ax-1]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
  
    let deriv_ta (tau:i64) : [c][Ax-1][c][Ax-1]f64 =
      let utils : trm.utility [ns][nd] = trm.utility mp p tau
      let tr = trm.age_transition mp
      let ev0 = trm.ev0 mp
      let f = trm.bellmanJ mp utils tr
      let param = eqb.dps.default
      let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
      let (ev, v) = trm.bellman0 mp utils tr ev
      in eqb.ded_dprice_tau mp tau ev v p

    let zeroes : [c][Ax-1][c][Ax-1]f64 = tabulate_2d c (Ax-1) (\_ _ -> replicate c (replicate (Ax-1) 0.0f64))
    let edfs : [n][c][Ax-1][c][Ax-1]f64 = #[sequential_outer] map deriv_ta (iota n)
    let edfs : [n][c][Ax-1][c][Ax-1]f64 = map2 (\edf tw -> map (map (map (map (\x -> f64.(x * tw))))) edf) edfs mp.tw
    in reduce (map2 (map2 (map2 (map2 (\x y -> f64.(x + y)))))) zeroes edfs

--- No input here, just for testing.
  entry test_edf_deriv_with_mum [n][c] (mum:[n]f64) (newprices:[c]f64) (Ax:i64) : [c][Ax-1][c][Ax-1]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let mp = trm.set_mum mum mp
    let p = trm.simple_prices mp 0.85
  
    let deriv_ta (tau:i64) : [c][Ax-1][c][Ax-1]f64 =
      let utils : trm.utility [ns][nd] = trm.utility mp p tau
      let tr = trm.age_transition mp
      let ev0 = trm.ev0 mp
      let f = trm.bellmanJ mp utils tr
      let param = eqb.dps.default
      let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
      let (ev, v) = trm.bellman0 mp utils tr ev
      in eqb.ded_dprice_tau mp tau ev v p

    let zeroes : [c][Ax-1][c][Ax-1]f64 = tabulate_2d c (Ax-1) (\_ _ -> replicate c (replicate (Ax-1) 0.0f64))
    let edfs : [n][c][Ax-1][c][Ax-1]f64 = #[sequential_outer] map deriv_ta (iota n)
    let edfs : [n][c][Ax-1][c][Ax-1]f64 = map2 (\edf tw -> map (map (map (map (\x -> f64.(x * tw))))) edf) edfs mp.tw
    in reduce (map2 (map2 (map2 (map2 (\x y -> f64.(x + y)))))) zeroes edfs

-- ==
-- entry: test_edf_tau_deriv
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[[-0.023859654165286f64], [0.023198151962656f64]]], 
--             [[[0.023198151962656f64], [-0.023859654165286f64]]]]}

entry test_edf_tau_deriv [c] (n:i64) (newprices:[c]f64) (Ax:i64) : [c][Ax-1][c][Ax-1]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    let (ev, v) = trm.bellman0 mp utils tr ev
    in eqb.ded_dprice_tau mp tau ev v p

-- ==
-- entry: test_edf_tau_man_deriv
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[[-0.023859654165286f64], [0.023198151962656f64]]], 
--             [[[0.023198151962656f64], [-0.023859654165286f64]]]]}

entry test_edf_tau_man_deriv [c] (n:i64) (newprices:[c]f64) (Ax:i64) : [c][Ax-1][c][Ax-1]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    let (ev, v) = trm.bellman0 mp utils tr ev
    in eqb.ded_dprice_tau_man mp tau ev v p

-- ==
-- entry: test_ccp_scrap_dprice_man
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[-0.00000001011f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64]], 
--             [[0.0f64, 0.0f64, -0.00000001011f64, 0.0f64, 0.0004f64]]]}

entry test_ccp_scrap_dprice_man [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[c][Ax-1][ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 0
  in eqb.ccp_scrap_dprice_man mp p tau

-- ==
-- entry: test_ergodic_with_dprice
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64]], 
--             [[-0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64]]]}

entry test_ergodic_with_dprice [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[c][Ax-1][ns]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    let ctp = eqb.ctp_from_utils mp tr utils ev
    let dctp = eqb.ctp_dprice_full mp p tau ev
    in eqb.ergodic_with_dprice ctp dctp|>(.1)

entry test_ctp_from_utils [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[ns][ns]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    in eqb.ctp_from_utils mp tr utils ev


-- Look at sparsity/simplification for the below later.
-- ==
-- entry: test_ctp_deriv_ad
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[[0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64], 
--              [0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64], 
--              [0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64],
--              [0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64],
--              [0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64]]], 
--            [[[-0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64], 
--              [-0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64], 
--              [-0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64],
--              [-0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64],
--              [-0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64]]]]}

entry test_ctp_deriv_ad [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[c][Ax-1][ns][ns]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    in eqb.ctp_dprice_full mp p tau ev

-- ==
-- entry: test_ccp_deriv_ad
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output {  [[[[-0.0001f64, 0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64], 
--              [0.0000f64, 0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64], 
--              [-0.0000f64,0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64],
--              [0.0000f64, 0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64],
--              [0.0000f64, 0.0238f64, -0.0001f64, -0.0232f64, -0.0000f64, -0.0004f64]]], 
--            [[[-0.0000f64, -0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64], 
--              [0.0000f64, -0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64], 
--              [-0.0000f64, -0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64],
--              [0.0000f64, -0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64],
--              [0.0000f64, -0.0233f64, -0.0000f64, 0.0238f64, -0.0001f64, -0.0004f64]]]]}

entry test_ccp_deriv_ad [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns][nd].[c][Ax-1][ns][nd]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    in eqb.ccp_dprice_full mp p tau ev

-- ==
-- entry: test_dev_prices
-- input { 2i64 [100f64, 100f64] 2i64 }
-- output { [[[1.0393f64, 0.9393f64, 0.9393f64, 0.9393f64, 0.9393f64]],
--           [[0.9393f64, 0.9393f64, 1.0393f64, 0.9393f64, 0.9393f64]]]}

entry test_dev_prices [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[c][Ax-1][ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 0
  let utils = trm.utility mp p tau
  let du = eqb.utility_dprice_ad mp p tau
  let ev0 = trm.ev0 mp
  let tr = trm.age_transition mp
  let f = trm.bellmanJ mp utils tr
  let param = eqb.dps.default
  let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
  let (ev, v) = trm.bellman0 mp utils tr ev
  let ctp = eqb.ctp_from_utils mp tr utils ev
  let ccp : [ns][nd]f64 = trm.ccp_tau mp v ev
  let dbellman = eqb.dbellman_prices_man ccp du
  in eqb.dev_fixed_point_dprice mp ctp dbellman

-- ==
-- entry: test_dbellman_prices_man
-- input { 1i64 [100f64] 3i64 }
-- output { [[[0.0919f64, -0.0039f64, -0.0042f64, -0.0042f64],
--           [-0.0088f64, 0.0832f64, -0.0092f64, -0.0092f64]]]}

entry test_dbellman_prices_man [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[c][Ax-1][ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 0
  let utils = trm.utility mp p tau
  let du = eqb.utility_dprice_ad mp p tau
  let ev0 = trm.ev0 mp
  let tr = trm.age_transition mp
  let (ev, v) = trm.bellman0 mp utils tr ev0
  let ccp : [ns][nd]f64 = trm.ccp_tau mp v ev
  let ccp = map (map (\x -> if f64.isnan x then f64.i64 0 else x)) ccp
  in eqb.dbellman_prices_man ccp du

-- ==
-- entry: test_dbellman_prices_ad
-- input { 1i64 [100f64] 3i64 }
-- output { [[[0.0919f64, -0.0039f64, -0.0042f64, -0.0042f64],
--           [-0.0088f64, 0.0832f64, -0.0092f64, -0.0092f64]]]}

entry test_dbellman_prices_ad [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[c][Ax-1][ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 0
  let utils = trm.utility mp p tau
  let du = eqb.utility_dprice_ad mp p tau
  let ev0 = trm.ev0 mp
  let tr = trm.age_transition mp
  in eqb.dbellman_prices_ad mp utils tr ev0 du

-- ==
-- entry: test_utility_dprice_ad
-- input { 1i64 [100f64] 3i64 }
-- output { [
--          [[[0.0f64, 0.1f64, 0.0f64, 0.1f64, 0.1f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64]],
--          [[0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.1f64, 0.1f64, -0.0f64, 0.1f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64]]]
--          ] }
entry test_utility_dprice_ad [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns][nd].[c][Ax-1][ns][nd]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 0
  in eqb.utility_dprice_ad mp p tau

-- ==
-- entry: test_utility_dprice_man
-- input { 1i64 [100f64] 3i64 }
-- output { [
--          [[[0.0f64, 0.1f64, 0.0f64, 0.1f64, 0.1f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64]],
--          [[0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.1f64, 0.1f64, -0.0f64, 0.1f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64]]]
--          ] }
entry test_utility_dprice_man [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns][nd].[c][Ax-1][ns][nd]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 0
  let ccp_scrap_tau = trm.ccp_scrap_tau mp p tau
  in eqb.utility_dprice_man mp ccp_scrap_tau tau