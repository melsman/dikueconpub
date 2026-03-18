import "trmodel"
import "../lib/github.com/diku-dk/linalg/lup"
import "../lib/github.com/diku-dk/linalg/lu"
import "../lib/github.com/diku-dk/linalg/dpsolve"

module equilibrium (R:real) (trm:trmodel with t = R.t) = {

    module lu = mk_lu R
    module lup = mk_lup R
    module dps = mk_dpsolve_dense R
    
    type t = R.t
    type mp [n][c][Ax][ns][nd] = trm.mp [n][c][Ax][ns][nd]
    type~ transition[ns] = trm.transition [ns]

    type~ sol [n][ns][nd][np] = 
        {p             : [np]t,         -- prices
         ed            : [np]t,         -- excess demand vector
         F             : transition[ns], -- struct with physical transition matrices
         ev_tau        : [n][ns]t,       -- expected values for all types
         ccp_tau       : [n][ns][nd]t,   -- conditional choice probabilities for all choices
         ccp_scrap_tau : [n][ns]t,       -- scrapping probabilities dependent on not having a car
         ctp_tau       : [n][ns][ns]t,   -- total transition matrix delta_tau*F.notrade
         delta_tau     : [n][ns][ns]t,   -- trade transition probability matrices for all consumer types
         h_tau         : [n][ns]t       -- post-trade equilibrium holdings for all consumer types
        }
    

    --- ergodic: finds the invariant distribution for an NxN Markov transition probability, ccp
    --- inputs
    --- ctp: a total transition matrix of size ns x ns
    --- outputs: The ergodic distribution [ns]t q, that satisfies q@ccp=q
    --- Note - could probably save some time here by seeing if I could implement some sort of solver than doesn't
    --- require sparse matrices?
    def ergodic [ns] (ctp:[ns][ns]t) : [ns]t =
        let ap = tabulate_2d (ns+1) (ns+1) (\i j -> if (i==ns || j==ns) then R.i64 1 else if i==j then R.(i64 1-ctp[j][i]) else R.(i64 0-ctp[j][i]))
        let ed0 = tabulate (ns+1) (\i -> if (i==ns) then R.i64 2 else R.i64 1)
        let res' = lup.ols ap ed0
        in map (\i->res'[i]) (iota ns)

    def ed_tau [n][c][Ax][ns][nd] (q_tau:[ns]t) (deltaK_diag:[ns]t) (deltaT_tau:[ns][ns]t) (_:mp[n][c][Ax][ns][nd]) (ccp_scrap_tau:[ns]t) : ?[np].[np]t = 
        let deltaT_trans = transpose deltaT_tau
        let demand = flatten (
            tabulate_2d c (Ax-1) 
                (\ct a-> let ind = ct*Ax+a
                 in map2 (\x y->R.(x*y)) deltaT_trans[ind] q_tau|>reduce (R.+) (R.i64 0)
                )
            )
        let sell_p = map (\x->R.(i64 1-x)) deltaK_diag
        let no_scrap = map (\x->R.(i64 1-x)) ccp_scrap_tau
        let supply = flatten (
            tabulate_2d c (Ax-1)
            (\ct a->let ind = ct*Ax+a
             in R.(sell_p[ind]*no_scrap[ind]*q_tau[ind])
            )
        )
        in map2 (\d s->R.(d-s)) demand supply

    -------------------------------
    ------ Derivatives - AD -------
    -------------------------------

    def price_basis [c][Ax] (cartype:i64) (age:i64) : trm.prices[c][Ax] =
        tabulate_2d c Ax (\cartype' age' -> R.bool (cartype == cartype' && age == age'))

    def utility_basis [ns][nd] (s:i64) (d:i64) : [ns][nd]t =
        tabulate_2d ns nd (\s' d' -> R.bool (s == s' && d == d'))

    ---- There is highly likely sparsity here that I am not taking advantage of. In particular,
    ---- the Jacobian is likely very sparse, since a chance in price for a car will only affect the utility of consumers who own that car
    ---- or buy that car
    def utility_dprice_ad [n][c][Ax][ns][nd]
        (mp: trm.mp[n][c][Ax][ns][nd])
        (p: trm.prices[c][Ax])
        (tau: i64)
        : [c][Ax-1][ns][nd]t =
        let f = \p' -> trm.utility mp p' tau
        -- Note that the prices for age-zero cars are known, and thus not endogenous, so we only take the derivative with respect
        -- to the other ages
        let du =
            tabulate_2d c (Ax-1) (\j a ->
            jvp f p (price_basis j (a+1)))
        in du

    def ctp_from_mp_p [n][c][Ax][ns][nd]
        (mp: mp[n][c][Ax][ns][nd])
        (p: trm.prices[c][Ax]) (tau:i64) (ev_s:trm.ev[ns]) : [ns][ns]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let (ev, v) = trm.bellman0 mp utils tr ev_s
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, _, _) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
        let ctp : [ns][ns]t = trm.ctp_tau tr delta
        in ctp
        
    def ctp_from_utils [n][c][Ax][ns][nd]
        (mp: trm.mp[n][c][Ax][ns][nd])
        (tr: trm.transition[ns])
        (utils: trm.utility[ns][nd])
        (ev_s: trm.ev[ns]) : [ns][ns]t =
        let (ev, v) = trm.bellman0 mp utils tr ev_s
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, _, _) = trm.trade_transition mp ccp
        in trm.ctp_tau tr delta


    --------------------------------
    ------ Derivaties - non-AD -----
    --------------------------------

    --- It doesn't seem like like utility of prices depend on the current price of the car - i.e. the derivative is linear.
    --- Will need to modify this if I implement ptransport.
    def utility_dprice_man [n][c][Ax][ns][nd]
            (mp: trm.mp[n][c][Ax][ns][nd])
            (tau: i64) : [c][Ax-1][ns][nd]t =
        let mu = mp.mum[tau]
        in tabulate_2d c (Ax-1) (\cartype age ->
            let row = cartype*Ax + age
            --- The first 1 added in col is because the first decision is keeping, and the reason (age+1) is used
            --- is because that in decisions, cars range from age 0 to Ax-1, but in states, they range from 1 + Ax.
            let col = 1 + cartype*Ax + (age+1)
            in tabulate_2d ns nd (\s d ->
                    let row_val = (if ((s == row) && not (d == 0)) then mu else R.i64 0)
                    let col_val = (if d == col then R.(i64 0-mu) else R.i64 0)
                    in R.(row_val + col_val)
                    )
                )

    def dbellman_prices_man [Ax][c][ns][nd]
        (ccp : [ns][nd]t)
        (du  : [c][Ax-1][ns][nd]t) : [c][Ax-1][ns]t =
        tabulate_2d c (Ax-1) (\cartype age ->
            tabulate ns (\s ->
            reduce (\i j -> R.(i+j)) (R.i64 0)
                (map2 (\i j -> R.(i*j)) ccp[s] du[cartype][age][s])
                )
            )

    def dev_fixed_point_dprice [n][c][Ax][ns][nd]
            (mp: mp[n][c][Ax][ns][nd])
            (ctp: [ns][ns]t)
            (dbellman: [c][Ax-1][ns]t) : [c][Ax-1][ns]t =
        let par_mat : [ns][ns]t =
            tabulate_2d ns ns (\i j ->
            let eye = R.bool (i == j)
            in R.(eye - mp.bet * ctp[i,j]))
        let blksz : i64 = 16
        in tabulate_2d c (Ax-1) (\cartype age ->
            lu.ols blksz par_mat dbellman[cartype][age])
        




    -- def ctp_dprice_via_utils_ad [n][c][Ax][ns][nd]
    --    (mp: trm.mp[n][c][Ax][ns][nd])
    --    (p: trm.prices[c][Ax])
    --    (tau: i64)
    --    (ev_s: trm.ev[ns])
    --    : [c][Ax-1][ns][ns]t =
    --    let tr = trm.age_transition mp
    --    let utils = trm.utility mp p tau
    --    let du_dp = utility_dprice_man mp tau
    --    let g = \u -> ctp_from_utils mp tr u ev_s
    --    in tabulate_2d c (Ax-1) (\cartype age ->
    --        jvp g utils du_dp[cartype][age])
        
    ---- test functions
     def utils_single [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) : [ns][nd]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        in utils

    def solve_single [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) (sa_max:i64) : [ns]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        in ev

    def edf_ccp_tau [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) (sa_max:i64) : [ns][nd]t  =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        let (_, v) = trm.bellman0 mp utils tr ev
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        in ccp

    def edf_q_delta [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) (sa_max:i64) : [ns][ns]t  =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        let (_, v) = trm.bellman0 mp utils tr ev
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, _, _) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
        in delta
    
    def edf_q_tau [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) (sa_max:i64) : [ns]t  =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        let (_, v) = trm.bellman0 mp utils tr ev
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, _, _) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
        let ctp : [ns][ns]t = trm.ctp_tau tr delta
        let q_tau : [ns]t = ergodic ctp
        in q_tau

    def edf_q_tau_2 [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) (sa_max:i64) : [ns]t  =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        let ctp = ctp_from_utils mp tr utils ev 
        let q_tau : [ns]t = ergodic ctp
        in q_tau

    def edf_tau_test [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (tau:i64) (sa_max:i64) : ?[np].[np]t  =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let ccp_scrap_tau = trm.ccp_scrap_tau mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        let (_, v) = trm.bellman0 mp utils tr ev
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, deltaK_diag, deltaT) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
        let ctp : [ns][ns]t = trm.ctp_tau tr delta
        let q_tau : [ns]t = ergodic ctp
        in ed_tau q_tau deltaK_diag deltaT mp ccp_scrap_tau

    def edf_test [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[c][Ax]t) (sa_max:i64) : ?[np].[np]t  =
        let comp (tau:i64) : [ns]t =
            let utils : trm.utility [ns][nd] = trm.utility mp p tau
            let tr = trm.age_transition mp
            let ev0 = trm.ev0 mp
            let f = trm.bellmanJ mp utils tr
            let param = dps.default with sa_max=sa_max
            let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
            in ev
        let evs : [n][ns]t =
            #[sequential_outer] map comp (iota n)
        let edf (ev:[ns]t) (tau:i64) =
            let utils : trm.utility [ns][nd] = trm.utility mp p tau
            let tr = trm.age_transition mp
            let ccp_scrap_tau = trm.ccp_scrap_tau mp p tau
            let (_, v) = trm.bellman0 mp utils tr ev
            let ccp : [ns][nd]t = trm.ccp_tau mp v ev
            let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
            let (delta, deltaK_diag, deltaT) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
            let ctp : [ns][ns]t = trm.ctp_tau tr delta
            let q_tau : [ns]t = ergodic ctp
            in ed_tau q_tau deltaK_diag deltaT mp ccp_scrap_tau
        let edfs = map2 edf evs (iota n)
        let tw = trm.get_tw mp
        let edfs_scaled = map2 (\x w-> map (\y->R.(y*w)) x) edfs tw
        in edfs_scaled|>reduce (map2 (\x y->R.(x+y))) (replicate (c*(Ax-1)) (R.i64 0))
}