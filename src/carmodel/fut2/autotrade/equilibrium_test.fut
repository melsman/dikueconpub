import "equilibrium"
import "trmodel"

module trm = trmodel f64
module eqb = equilibrium f64 trm

entry test_spp_price_solve [n][c] (mum:[n]f64) (newprices:[c]f64) (Ax:i64) (acc0:f64) (Axb:i64) : [c][Ax]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let mp = trm.set_mum mum mp
    let mp = trm.set_acc_0 (replicate c acc0) mp
    in eqb.spp_price_solve mp Axb


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

entry test_ergodic_dist [c][n] (mum:[n]f64) (newprices:[c]f64) (Ax:i64) (acc0:f64) : ?[ns].[ns]f64 =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let mp = trm.set_mum mum mp
    let mp = trm.set_acc_0 (replicate c acc0) mp
    let p = trm.simple_prices mp 0.85
    let tau = 0
    let utils : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp utils tr
    let param = eqb.dps.default
    let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = eqb.dps.poly f ev0 param (f64.f32 0)
    let ctp = eqb.ctp_from_utils mp tr utils ev
    in eqb.ergodic ctp

-- ==
-- entry: test_utils_single 
-- input {2i64 2i64 2i64 }
-- output { [[5.5000f64, 4.5000f64, 5.5000f64, 4.5000f64, 5.5000f64, 8.5000f64],
--           [f64.nan, -3.9000f64, -2.9000f64, -3.9000f64, -2.9000f64, 0.1000f64],
--           [5.5000f64, 4.5000f64, 5.5000f64, 4.5000f64, 5.5000f64, 8.5000f64],
--           [f64.nan, -3.9000f64, -2.9000f64, -3.9000f64, -2.9000f64, 0.1000f64],
--           [f64.nan, -4.0000f64, -3.0000f64, -4.0000f64, -3.0000f64, 0.0000f64]]}

entry test_utils_single  (n:i64) (c:i64) (Ax:i64) : ?[ns][nd].[ns][nd]f64 =
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.utils_single mp p 1

-- ==
-- entry: test_solve_single 
-- input {2i64 2i64 2i64 }
-- output { [104.0525f64, 95.6521f64, 104.0525f64, 95.6521f64, 95.6521f64] }

entry test_solve_single  (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.solve_single mp p 1 sa_max

-- ==
-- entry: test_edf_ccp
-- input {2i64 2i64 2i64 }
-- output { [[0.0005f64, 0.4951f64, 0.0005f64, 0.4951f64, 0.0005f64, 0.0084f64],
--          [0.0f64, 0.4951f64, 0.0005f64, 0.4951f64, 0.0005f64, 0.0084f64],
--          [0.0005f64, 0.4951f64, 0.0005f64, 0.4951f64, 0.0005f64, 0.0084f64],
--          [0.0f64, 0.4951f64, 0.0005f64, 0.4951f64, 0.0005f64, 0.0084f64],
--          [0.0f64, 0.4951f64, 0.0005f64, 0.4951f64, 0.0005f64, 0.0084f64]] }


entry test_edf_ccp (n:i64) (c:i64) (Ax:i64) : ?[ns][nd].[ns][nd]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.edf_ccp_tau mp p 1 sa_max

-- ==
-- entry: test_edf_delta
-- input {2i64 2i64 2i64 }
-- output { [[0.0009f64, 0.4951f64, 0.0005f64, 0.4951f64, 0.0084f64],
--          [0.0005f64, 0.4953f64, 0.0005f64, 0.4953f64, 0.0084f64],
--          [0.0005f64, 0.4951f64, 0.0009f64, 0.4951f64, 0.0084f64],
--          [0.0005f64, 0.4953f64, 0.0005f64, 0.4953f64, 0.0084f64],
--          [0.0005f64, 0.4953f64, 0.0005f64, 0.4953f64, 0.0084f64]] }

entry test_edf_delta (n:i64) (c:i64) (Ax:i64) : ?[ns][nd].[ns][nd]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    let ccp = eqb.edf_ccp_tau mp p 1 sa_max
    let ccp = map (map (\x -> if f64.isnan x then f64.i64 0 else x)) ccp
    let (delta, _, _) = trm.trade_transition mp ccp
    in delta

-- ==
-- entry: test_edf_q_tau
-- input {2i64 2i64 2i64 }
-- output { [0.4951f64, 0.0007f64, 0.4951f64, 0.0007f64, 0.0084f64] }

entry test_edf_q_tau (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.edf_q_tau mp p 1 sa_max

-- ==
-- entry: test_edf_q_tau_2
-- input {2i64 2i64 2i64 }
-- output { [0.4951f64, 0.0007f64, 0.4951f64, 0.0007f64, 0.0084f64] }

entry test_edf_q_tau_2 (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.edf_q_tau_2 mp p 1 sa_max

-- ==
-- entry: test_edf_tau_test
-- input {2i64 2i64 2i64 }
-- output { [-0.4944f64, -0.4944f64] }

entry test_edf_tau_test (n:i64) (c:i64) (Ax:i64) : ?[np].[np]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.edf_tau_test mp p 1 sa_max

-- ==
-- entry: test_edf_test
-- input {2i64 2i64 2i64 }
-- output { [-0.4944f64, -0.4944f64] }

entry test_edf_test (n:i64) (c:i64) (Ax:i64) : ?[np].[np]f64 =
    let sa_max = 3
    let pnews : [c]f64 = replicate c 100.0f64
    let mp = trm.mk n c Ax
    let mp = trm.set_newprices mp pnews
    let p = trm.simple_prices mp (f64.f32 0.85)
    in eqb.edf_test mp p sa_max