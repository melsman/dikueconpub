import "equilibrium"
import "trmodel"

module trm = trmodel f64
module eqb = equilibrium f64 trm


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
  let tau = 0
  in eqb.utility_dprice_man mp tau