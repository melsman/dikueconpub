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
