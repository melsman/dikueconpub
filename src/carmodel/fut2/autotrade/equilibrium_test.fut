import "equilibrium"
import "trmodel"

module trm = trmodel f64
module eqb = equilibrium f64 trm

-- ==
-- entry: test_ergodic
-- input {[ [0.4953f64, 0.0009f64, 0.4953f64, 0.0005f64, 0.0081f64],
--          [0.4955f64, 0.0005f64, 0.4955f64, 0.0005f64, 0.0081f64],
--          [0.4953f64, 0.0005f64, 0.4953f64, 0.0009f64, 0.0081f64],
--          [0.4955f64, 0.0005f64, 0.4955f64, 0.0005f64, 0.0081f64],
--          [0.4955f64, 0.0005f64, 0.4955f64, 0.0005f64, 0.0081f64], ] }
-- output { [0.4953f64, 0.0007f64, 0.4953f64, 0.0007f64, 0.0081f64] }

entry test_ergodic [ns] (ccp:[ns][ns]f64) : [ns]f64 =
    eqb.ergodic ccp

-- ==
-- entry: test_utils_single 
-- input {2i64 2i64 2i64 }
-- output { [[5.5000f64, 4.5000f64, 5.5000f64, 4.5000f64, 5.5000f64, 8.5000f64],
--           [f64.nan, -3.9000f64, -2.9000f64, -3.9000f64, -2.9000f64, 0.1000f64],
--           [5.5000f64, 4.5000f64, 5.5000f64, 4.5000f64, 5.5000f64, 8.5000f64],
--           [f64.nan, -3.9000f64, -2.9000f64, -3.9000f64, -2.9000f64, 0.1000f64],
--           [f64.nan, -4.0000f64, -3.0000f64, -4.0000f64, -3.0000f64, 0.0000f64]]}

entry test_utils_single  (n:i64) (c:i64) (Ax:i64) : ?[ns][nd].[ns][nd]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.utils_single mp p 1 sa_max

-- ==
-- entry: test_solve_single 
-- input {2i64 2i64 2i64 }
-- output { [68.0597f64, 59.6592f64, 68.0597f64, 59.6592f64, 59.5592f64] }

entry test_solve_single  (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.solve_single mp p 1 sa_max

-- ==
-- entry: test_edf_q_delta
-- input {2i64 2i64 2i64 }
-- output { [ [0.4953f64, 0.0009f64, 0.4953f64, 0.0005f64, 0.0081f64],
--          [0.4955f64, 0.0005f64, 0.4955f64, 0.0005f64, 0.0081f64],
--          [0.4953f64, 0.0005f64, 0.4953f64, 0.0009f64, 0.0081f64],
--          [0.4955f64, 0.0005f64, 0.4955f64, 0.0005f64, 0.0081f64],
--          [0.4955f64, 0.0005f64, 0.4955f64, 0.0005f64, 0.0081f64], ] }

entry test_edf_q_delta (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns][ns]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let u_0 = tabulate_2d n c
            (\i j -> f64.(i64 5 + i64 2 * (i64 i + i64 j) / (i64 n + i64 c)))
    let mp = trm.mk n c Ax |> trm.set_u_0 u_0
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.edf_q_delta mp p 1 sa_max
