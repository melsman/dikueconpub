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
    
    ------------------------------------------------
    ------- Social planner functions ---------------
    ------------------------------------------------

    def bellman_spp_with_deriv_p [Axb][n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd]) (tau: i64) (car: i64) (V:[Axb]t)
                                                                                            : ([Axb]t, [Axb][Axb]t, [Axb]i64) =
        let end = Axb - 1
        let apr : [Axb]t = tabulate Axb (\a -> trm.acc_prob_j mp (i32.i64 a) (i32.i64 car) (i32.i64 Axb))
        -- Utility from replacing with a new car (age 0)
        let urep : t = R.(trm.u_car mp tau {cartype=car, age=0} - mp.mum[tau] * (mp.pnew[car] - mp.pscrap[car]))
        -- Utility from keeping a car of age a, for a=0..end-1.
        let ukeep : [Axb-1]t = tabulate (Axb-1) (\a -> trm.u_car mp tau {cartype=car, age=a})
        -- Value of replacement.
        let vrep : t = R.(urep + mp.bet * ((i64 1 - apr[0]) * V[1] + apr[0] * V[end]))
        -- Value of keeping age-a car.
        let vkeep : [Axb-1]t = tabulate (Axb-1) (\a -> 
                let a_next = a+1
                in R.(ukeep[a] + mp.bet * ((i64 1 - apr[a]) * V[a_next] + apr[a] * V[end]))
            )

        -- Updated Bellman value.
        let Vnew : [Axb]t =
            tabulate Axb (\a ->
                if a == end then vrep
                else R.max vrep vkeep[a]
        )
        -- Optimal replacement indicator:
        -- P[a] = 1 if replace, 0 if keep.
        let P : [Axb]i64 =
            tabulate Axb (\a ->
                if a == end then 1
                else i64.bool R.(vrep > vkeep[a])
            )

        -- Transition matrix
        let tpm : [Axb][Axb]t =
            tabulate_2d Axb Axb (\i j ->
            let Pi = R.i64 P[i]
            let diag = if (i < end && i+1==j) then R.((i64 1-Pi)*(i64 1-apr[i])) else R.i64 0
            let to_end = if j==end then R.(Pi*apr[0]+(i64 1-Pi)*apr[i]) else R.i64 0
            let to_age1 = if j==1 then R.(Pi*(i64 1 - apr[0])) else R.i64 0
            in R.(diag + to_end + to_age1)
            )

        let dV : [Axb][Axb]t = map (map (R.* mp.bet)) tpm
        in (Vnew, dV, P)

    def bellman_spp_with_deriv [Axb][n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd]) (tau: i64) (car: i64) (V:[Axb]t) : ([Axb]t, [Axb][Axb]t) =
        let (Vnew, dv, _) =  bellman_spp_with_deriv_p mp tau car V
        in (Vnew, dv)
    
    def solve_spp_single [Axb][n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (tau:i64) (car:i64) (V0:[Axb]t) : ([Axb]t) =
        let param = dps.default
        let f = bellman_spp_with_deriv mp tau car
        let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f V0 param (mp.bet)
        in ev
    
    ----- Note that this is to some extent a combination of spp.solve_spp and equilibrium.price_init, in 
    ----- particular since some changes had to be made to prevent irregular parallelism.
    def spp_price_solve [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (Axb:i64) : [c][Ax]t =
        let V0 = replicate Axb (R.i64 0)
        let Vs = #[sequential_outer] tabulate_2d n c (\t ct -> solve_spp_single mp t ct V0)
        let policies = tabulate_2d n c  (\t ct ->
            let (_, _, P) = bellman_spp_with_deriv_p mp t ct Vs[t][ct]
            in P
        )
        let abar_spp = tabulate_2d n c (\t ct->
            let first_replace_age = reduce (i64.min) Axb (map (\j-> if policies[t][ct][j]==1 then j else Axb) (iota (Axb-1)))
            in if first_replace_age == Axb then Ax else first_replace_age
        )
        let p_spp = tabulate_3d n c Axb (\t ct a ->
            R.(mp.pnew[ct]-(Vs[t][ct][0]-Vs[t][ct][a])/mp.mum[t])
        )
        let abar_and_tau_j = tabulate c (\ct ->
            let tau_j = reduce (\tc tn-> if abar_spp[tn][ct]>abar_spp[tc][ct] then tn else tc) 0 (iota n)
            in (abar_spp[tau_j][ct], tau_j)
        )
        let initial_prices = tabulate_2d c Ax (\ct a->
            let (abar, t) = abar_and_tau_j[ct]
            in if a==0 then mp.pnew[ct]
            else if a>=abar then mp.pscrap[ct]
            else p_spp[t][ct][a]
        )
        in initial_prices



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

    def demand_supply_tau [n][c][Ax][ns][nd] (q_tau:[ns]t) (deltaK_diag:[ns]t) (deltaT_tau:[ns][ns]t) (_:mp[n][c][Ax][ns][nd]) (ccp_scrap_tau:[ns]t) : ([c][Ax-1]t, [c][Ax-1]t) =
        let deltaT_trans = transpose deltaT_tau
        let demand = tabulate_2d c (Ax-1) (\ct a->
                let ind = ct*Ax+a
                in map2 (\x y->R.(x*y)) deltaT_trans[ind] q_tau|>reduce (R.+) (R.i64 0)
                )
        let sell_p = map (\x->R.(i64 1-x)) deltaK_diag
        let no_scrap = map (\x->R.(i64 1-x)) ccp_scrap_tau
        let supply = tabulate_2d c (Ax-1) (\ct a->
            let ind = ct*Ax+a
            in R.(sell_p[ind]*no_scrap[ind]*q_tau[ind])
            )
        in (demand, supply)

    def ed_tau [n][c][Ax][ns][nd] (q_tau:[ns]t) (deltaK_diag:[ns]t) (deltaT_tau:[ns][ns]t) (mp:mp[n][c][Ax][ns][nd]) (ccp_scrap_tau:[ns]t) : [c][Ax-1]t =
        let (demand, supply) = demand_supply_tau q_tau deltaK_diag deltaT_tau mp ccp_scrap_tau
        in map2 (map2 (\x y->R.(x-y))) demand supply


    ----------------------------------------------
    ------ Helper functions for derivatives ------
    ----------------------------------------------

    def ctp_from_mp_p [n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd])
        (p: trm.prices[c][Ax]) (tau:i64) (ev_s:trm.ev[ns]) : [ns][ns]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let (ev, v) = trm.bellman0 mp utils tr ev_s
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, _, _) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
        let ctp : [ns][ns]t = trm.ctp_tau tr delta
        in ctp
        
    def ctp_from_utils [n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd]) (tr: trm.transition[ns])
        (utils: trm.utility[ns][nd])  (ev_s: trm.ev[ns]) : [ns][ns]t =
        let (ev, v) = trm.bellman0 mp utils tr ev_s
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        let (delta, _, _) = trm.trade_transition mp ccp
        in trm.ctp_tau tr delta

    def ccp_from_utils [n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd]) (tr: trm.transition[ns])
        (utils: trm.utility[ns][nd]) (ev_s: trm.ev[ns]) : [ns][nd]t =
        let (ev, v) = trm.bellman0 mp utils tr ev_s
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev
        let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
        in ccp

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



    --------------------------------
    ------ Derivaties - non-AD -----
    --------------------------------

    --- It doesn't seem like like utility of prices depend on the current price of the car - i.e. the derivative is linear.
    --- Will need to modify this if I implement ptransport and ptc_sale
    def utility_dprice_man [n][c][Ax][ns][nd]
            (mp: trm.mp[n][c][Ax][ns][nd])
            (ccp_scrap: [ns]t)
            (tau: i64) : [c][Ax-1][ns][nd]t =
        let mu = mp.mum[tau]
        in tabulate_2d c (Ax-1) (\cartype age ->
            let row = cartype*Ax + age
            --- The first 1 added in col is because the first decision is keeping, and the reason (age+1) is used
            --- is because that in decisions, cars range from age 0 to Ax-1, but in states, they range from 1 + Ax.
            let col = 1 + cartype*Ax + (age+1)
            --- sell-price enters utility both directly (mum*psell) and via ev_scrap (whose derivative is -mum*ccp_scrap).
            let sell_mu = R.(mu * (i64 1 - ccp_scrap[row]))
            in tabulate_2d ns nd (\s d ->
                    let row_val = (if ((s == row) && not (d == 0)) then sell_mu else R.i64 0)
                    let col_val = (if d == col then R.(i64 0-mu) else R.i64 0)
                    in R.(row_val + col_val)
                    )
                )

    def ccp_scrap_dprice_man [n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd]) (p:trm.prices[c][Ax]) (tau:i64) : [c][Ax-1][ns]t =
        let mu = mp.mum[tau]
        let ccp_scrap_tau = trm.ccp_scrap_tau mp p tau
        in tabulate_2d c (Ax-1) (\cartype age ->
            let row = cartype*Ax+age
            in tabulate ns (\state ->
                if state!=row then R.i64 0
                else
                    let dev_scrap =  R.((i64 0-ccp_scrap_tau[state]*mu))
                    in R.((i64 1-ccp_scrap_tau[state])*dev_scrap/mp.sigma_s)
            )
        )


    def dbellman_prices_man [Ax][c][ns][nd] (ccp : [ns][nd]t) (du  : [c][Ax-1][ns][nd]t) : [c][Ax-1][ns]t =
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
    
    def dccp_dprices_from_du_dev [n][c][Ax][ns][nd]
            (mp: trm.mp[n][c][Ax][ns][nd])
            (tr: trm.transition[ns])
            (utils: trm.utility[ns][nd])
            (ev: trm.ev[ns])
            (du: [c][Ax-1][ns][nd]t)
            (dev: [c][Ax-1][ns]t) : [c][Ax-1][ns][nd]t =
        let g = \(u, e) -> ccp_from_utils mp tr u e
        in tabulate_2d c (Ax-1) (\cartype age ->
            jvp g (utils, ev) (du[cartype][age], dev[cartype][age]))

    def ccp_dprice_full [n][c][Ax][ns][nd] (mp: trm.mp[n][c][Ax][ns][nd]) (p: trm.prices[c][Ax])
        (tau: i64) (ev: trm.ev[ns])  : [c][Ax-1][ns][nd]t =
        let tr = trm.age_transition mp
        let utils : trm.utility[ns][nd] = trm.utility mp p tau
        let ctp : [ns][ns]t = ctp_from_utils mp tr utils ev
        let ccp_scrap = trm.ccp_scrap_tau mp p tau

        let du : [c][Ax-1][ns][nd]t =
            utility_dprice_man mp ccp_scrap tau

        let (ev1, v) = trm.bellman0 mp utils tr ev
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev1

        let dbellman : [c][Ax-1][ns]t =
            dbellman_prices_man ccp du

        let dev : [c][Ax-1][ns]t =
            dev_fixed_point_dprice mp ctp dbellman

        in dccp_dprices_from_du_dev mp tr utils ev du dev
        

    def dctp_dprices_from_du_dev [n][c][Ax][ns][nd] (mp: mp[n][c][Ax][ns][nd]) (tr: trm.transition[ns]) (utils: trm.utility[ns][nd])
            (ev: trm.ev[ns]) (du: [c][Ax-1][ns][nd]t) (dev: [c][Ax-1][ns]t) : [c][Ax-1][ns][ns]t =
        let g = \(u, e) -> ctp_from_utils mp tr u e
        in tabulate_2d c (Ax-1) (\cartype age ->
            jvp g (utils, ev) (du[cartype][age], dev[cartype][age]))

    def ctp_dprice_full [n][c][Ax][ns][nd] (mp: trm.mp[n][c][Ax][ns][nd]) (p: trm.prices[c][Ax])
        (tau: i64) (ev: trm.ev[ns])  : [c][Ax-1][ns][ns]t =
        let tr = trm.age_transition mp
        let utils : trm.utility[ns][nd] = trm.utility mp p tau
        let ctp : [ns][ns]t = ctp_from_utils mp tr utils ev
        let ccp_scrap = trm.ccp_scrap_tau mp p tau

        let du : [c][Ax-1][ns][nd]t =
            utility_dprice_man mp ccp_scrap tau

        let (ev1, v) = trm.bellman0 mp utils tr ev
        let ccp : [ns][nd]t = trm.ccp_tau mp v ev1

        let dbellman : [c][Ax-1][ns]t =
            dbellman_prices_man ccp du

        let dev : [c][Ax-1][ns]t =
            dev_fixed_point_dprice mp ctp dbellman

        in dctp_dprices_from_du_dev mp tr utils ev du dev

    def ergodic_with_dprice [c][Ax][ns] (ctp:[ns][ns]t) (dctp: [c][Ax-1][ns][ns]t) : ([ns]t, [c][Ax-1][ns]t) =
        let ap = tabulate_2d (ns+1) (ns+1) (\i j -> if (i==ns || j==ns) then R.i64 1 else if i==j then R.(i64 1-ctp[j][i]) else R.(i64 0-ctp[j][i]))
        let ed0 = tabulate (ns+1) (\i -> if (i==ns) then R.i64 2 else R.i64 1)
        let blksz = 16
        let erg' = lu.ols blksz ap ed0
        let erg = map (\i->erg'[i]) (iota ns)
        let erg_dprice = tabulate_2d c (Ax-1) (\cartype age ->
            let da : [ns+1][ns+1]t = tabulate_2d (ns+1) (ns+1) (\i j-> if (i==ns || j==ns) then R.i64 0 else R.(i64 0-dctp[cartype][age][j][i]))
            let minus_da_inva_ed = map (\row -> reduce (R.+) (R.i64 0) (map2 (\x y -> R.(i64 0-x*y)) row erg')) da
            let res =  lu.ols blksz ap minus_da_inva_ed
            in map (\i->res[i]) (iota ns)
        )
        in (erg, erg_dprice)


    def ded_dprice [c][Ax][ns][nd] (ccp: [ns][nd]t) (dccp: [c][Ax-1][ns][nd]t)  (q: [ns]t)
        (dq: [c][Ax-1][ns]t) (ccp_scrap: [ns]t) (dccp_scrap: [c][Ax-1][ns]t) : [c][Ax-1][c][Ax-1]t =

        tabulate_2d c (Ax-1) (\m_ct m_age ->
            let owned_car = m_ct*Ax + m_age
            let buycol = 1 + m_ct*Ax + (m_age+1)
            let keep_survive = R.(i64 1 - ccp[owned_car,0])
            let no_scrap = R.(i64 1 - ccp_scrap[owned_car])

            in tabulate_2d c (Ax-1) (\p_ct p_age ->

                let dD =
                    reduce (R.+) (R.i64 0)
                        (tabulate ns (\s ->
                            R.( q[s] * dccp[p_ct][p_age][s][buycol]
                            + ccp[s][buycol] * dq[p_ct][p_age][s] )
                            )
                        )

                let dED_supply_keep =
                    R.(dccp[p_ct][p_age][owned_car][0] * no_scrap * q[owned_car])

                let dED_supply_q =
                    R.(i64 0 - keep_survive * no_scrap * dq[p_ct][p_age][owned_car])

                let dED_supply_scrap =
                    R.(keep_survive * dccp_scrap[p_ct][p_age][owned_car] * q[owned_car])

                in R.(dD + dED_supply_keep + dED_supply_q + dED_supply_scrap)
            )
        )

    def ded_dprice_tau [n][c][Ax][ns][nd] (mp: trm.mp[n][c][Ax][ns][nd]) (tau: i64) (ev:[ns]t) 
        (v:[ns][nd]t) (p:trm.prices[c][Ax]) : [c][Ax-1][c][Ax-1]t =
        let utils : trm.utility [ns][nd] = trm.utility mp p tau
        let tr = trm.age_transition mp
        let ctp = ctp_from_utils mp tr utils ev
        let ccp_scrap = trm.ccp_scrap_tau mp p tau
        let ccp = trm.ccp_tau mp v ev

        let du = utility_dprice_man mp ccp_scrap tau
        let dbellman = dbellman_prices_man ccp du
        let dev = dev_fixed_point_dprice mp ctp dbellman
        let dctp = dctp_dprices_from_du_dev mp tr utils ev du dev
        let dccp = dccp_dprices_from_du_dev mp tr utils ev du dev
        let (erg, erg_dprice) = ergodic_with_dprice ctp dctp
        let dccp_scrap = ccp_scrap_dprice_man mp p tau

        in ded_dprice ccp dccp erg erg_dprice ccp_scrap dccp_scrap

    -------- Flattening functions for testing --------
    def flatten_ded [c][Ax] (ded: [c][Ax-1][c][Ax-1]t) : [c*(Ax-1)][c*(Ax-1)]t =
        map flatten (flatten ded)

    -------- Function from getting ed and ed_dprice from mp, p and sa_max --------
    def ed_ded_price_all [n][c][Ax][ns][nd] (mp:trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) (p:trm.prices[c][Ax]) : ([c][Ax-1]t, [c][Ax-1][c][Ax-1]t) =
        let sol_evs (tau:i64) : ([ns]t, [ns][nd]t) =
            let utils : trm.utility [ns][nd] = trm.utility mp p tau
            let tr = trm.age_transition mp
            let ev0 = trm.ev0 mp
            let f = trm.bellmanJ mp utils tr
            let param = dps.default
            let param = param with sa_max = sa_max
            let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.i64 0)
            let (_, v) = trm.bellman0 mp utils tr ev
            in (ev, v)

        let evs = #[sequential_outer] map sol_evs (iota n)
        let deds = #[sequential_outer] map2 (\evs tau -> 
            let (ev, v) = evs
            in ded_dprice_tau mp tau ev v p) evs (iota n)


        let deds : [n][c][Ax-1][c][Ax-1]t = map2 (\ded tw -> map (map (map (map (\x -> R.(x * tw))))) ded) deds mp.tw
        let zeroes : [c][Ax-1][c][Ax-1]t = tabulate_2d c (Ax-1) (\_ _ -> replicate c (replicate (Ax-1) (R.i64 0)))
        let ded = reduce (map2 (map2 (map2 (map2 (\x y -> R.(x + y)))))) zeroes deds

        let edf (evs:([ns]t, [ns][nd]t)) (tau:i64) =
            let (ev, v) = evs
            let tr = trm.age_transition mp
            let ccp_scrap_tau = trm.ccp_scrap_tau mp p tau
            let ccp : [ns][nd]t = trm.ccp_tau mp v ev
            let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
            let (delta, deltaK_diag, deltaT) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
            let ctp : [ns][ns]t = trm.ctp_tau tr delta
            let q_tau : [ns]t = ergodic ctp
            in ed_tau q_tau deltaK_diag deltaT mp ccp_scrap_tau

        let edfs = map2 edf evs (iota n)
        let edfs_scaled = map2 (\x w -> map (map (\y -> R.(y * w))) x) edfs mp.tw |> reduce (map2 (map2 (\x y -> R.(x + y)))) (replicate c (replicate (Ax-1) (R.i64 0)))

        in (edfs_scaled, ded)

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
        in ed_tau q_tau deltaK_diag deltaT mp ccp_scrap_tau|>flatten

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
            in ed_tau q_tau deltaK_diag deltaT mp ccp_scrap_tau|>flatten
        let edfs = map2 edf evs (iota n)
        let tw = trm.get_tw mp
        let edfs_scaled = map2 (\x w-> map (\y->R.(y*w)) x) edfs tw
        in edfs_scaled|>reduce (map2 (\x y->R.(x+y))) (replicate (c*(Ax-1)) (R.i64 0))

    -- Returns tw-weighted demand and supply separately, for diagnosing whether
    -- demand[ct][a] and supply[ct][a] refer to the same car-age market.
    def demand_supply_all [n][c][Ax][ns][nd] (mp:trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) (p:trm.prices[c][Ax]) : ([c][Ax-1]t, [c][Ax-1]t) =
        let sol_evs (tau:i64) : ([ns]t, [ns][nd]t) =
            let utils : trm.utility [ns][nd] = trm.utility mp p tau
            let tr = trm.age_transition mp
            let ev0 = trm.ev0 mp
            let f = trm.bellmanJ mp utils tr
            let param = dps.default with sa_max = sa_max
            let {res=ev,jac=_,conv=_,iter_sa=_,iter_nk=_,rtrips=_,tol=_} = dps.poly f ev0 param (R.i64 0)
            let (_, v) = trm.bellman0 mp utils tr ev
            in (ev, v)
        let evs = #[sequential_outer] map sol_evs (iota n)
        let dsf (evs:([ns]t, [ns][nd]t)) (tau:i64) =
            let (ev, v) = evs
            let tr = trm.age_transition mp
            let ccp_scrap_tau = trm.ccp_scrap_tau mp p tau
            let ccp : [ns][nd]t = trm.ccp_tau mp v ev
            let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
            let (delta, deltaK_diag, deltaT) : ([ns][ns]t, [ns]t, [ns][ns]t) = trm.trade_transition mp ccp
            let ctp : [ns][ns]t = trm.ctp_tau tr delta
            let q_tau : [ns]t = ergodic ctp
            in demand_supply_tau q_tau deltaK_diag deltaT mp ccp_scrap_tau
        let ds_per_type = map2 dsf evs (iota n)
        let zero : [c][Ax-1]t = replicate c (replicate (Ax-1) (R.i64 0))
        let (demands, supplies) = unzip ds_per_type
        let demand = map2 (\d tw -> map (map (\x -> R.(x * tw))) d) demands mp.tw
                     |> reduce (map2 (map2 (\x y -> R.(x + y)))) zero
        let supply = map2 (\s tw -> map (map (\x -> R.(x * tw))) s) supplies mp.tw
                     |> reduce (map2 (map2 (\x y -> R.(x + y)))) zero
        in (demand, supply)
}