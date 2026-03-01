import "trmodel"
import "../lib/github.com/diku-dk/linalg/lup"
import "../lib/github.com/diku-dk/linalg/dpsolve"

module equilibrium (R:real) (trm:trmodel with t = R.t) = {

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

    ------- edf functions
    --- this one is mostly for testing
    --- a bit annoying that I have to calculate ccp and delta again outside bellman, but it seems to be necessary since
    --- dpsolver requires two outputs
    --def edf_single [n][c][Ax][ns][nd][np] (mp:mp[n][c][Ax][ns][nd]) (p:[np]t) (t:tau) : ((trm.utility[ns][nd]), [ns][nd]t, [ns][ns]t, [ns]t, [np]t)  =
    --    let utils = trm.utility [ns][nd] = trm.utility mp p tau
    --    let tr = trm.age_transition mp
    --    let ev0 = trm.ev0 mp
    --    let f = trm.bellmanJ mp util tr
    --    let param = dps.default with sa_max=sa_max
    --    let {res=ev,jac=_,conv=b,iter_sa=i,iter_nk=j,rtrips=k,tol} = dps.poly f ev0 param (R.f32 0)
    --    let (_, v) = trm.bellman0 mp utils tr ev
    --    let ccp : [ns][nd]t = map2 (\r x -> map (R.exp <-< (R./ mp.sigma) <-< (R.- x)) r) v ev
    --    let delta : [ns][ns]t = trm.trade_transition mp ccp
    --    let ctp : [ns][ns]t = ctp_tau tr delta
    --    let q_tau : [ns][t] = ergodic ctp
    --    let 

    ---- test functions
     def utils_single [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[Ax][c]t) (tau:i64) (sa_max:i64) : [ns][nd]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        in utils

    def solve_single [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[Ax][c]t) (tau:i64) (sa_max:i64) : [ns]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ev0 = trm.ev0 mp
        let f = trm.bellmanJ mp utils tr
        let param = dps.default with sa_max=sa_max
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.f32 0)
        in ev

    def edf_ccp_tau [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[Ax][c]t) (tau:i64) (sa_max:i64) : [ns][nd]t  =
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

    def edf_q_delta [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[Ax][c]t) (tau:i64) (sa_max:i64) : [ns][ns]t  =
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
    
    def edf_q_tau [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:[Ax][c]t) (tau:i64) (sa_max:i64) : [ns]t  =
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
}