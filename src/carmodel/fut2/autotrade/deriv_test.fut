import "equilibrium"
import "trmodel"

module trm = trmodel f64
module eqb = equilibrium f64 trm

-- ==
-- entry: test_utility_dprice_ad
-- input { 1i64 [100f64] 3i64 }
-- output { [
--          [[[0.0f64, 0.1f64, 0.0f64, 0.1f64, 0.1f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64],
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64]]],
--          [[[0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.1f64, 0.1f64, -0.0f64, 0.1f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64]]]
--          ] }
entry test_utility_dprice_ad [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns][nd].[Ax-1][c][ns][nd]f64 =
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
--           [0.0f64, 0.0f64, -0.1f64, 0.0f64, 0.0f64]]],
--          [[[0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.1f64, 0.1f64, -0.0f64, 0.1f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64],
--           [0.0f64, 0.0f64, 0.0f64, -0.1f64, 0.0f64]]]
--          ] }
entry test_utility_dprice_man [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns][nd].[Ax-1][c][ns][nd]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let tau = 0
  in eqb.utility_dprice_man mp tau