import "trmodel"
import "../lib/github.com/diku-dk/linalg/lup"

module type equilibrium = {

    type t
    type mp [n][c][Ax][ns][nd]
    type~transition[ns]

    type~ sol[n][ns][nd][np]

    val ergodic [ns] : [ns][ns]t ->[ns]t
}

module equilibrium (R:real) (trm:trmodel with t = R.t) : equilibrium with t = R.t = {

    module lup = mk_lup R
    
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
    
    def ergodic [ns] (ccp:[ns][ns]t) : [ns]t =
        let ap = tabulate_2d (ns+1) (ns+1) (\i j -> if (i==ns || j==ns) then R.i64 1 else if i==j then R.(i64 1-ccp[j][i]) else R.(i64 0-ccp[j][i]))
        let ed0 = tabulate (ns+1) (\i -> if (i==ns) then R.i64 2 else R.i64 1)
        let res' = lup.ols ap ed0
        in map (\i->res'[i]) (iota ns)

}