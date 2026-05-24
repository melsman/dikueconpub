
import "autotrade/trmodel"
import "autotrade/equilibrium"
import "lib/github.com/diku-dk/linalg/lu"

module type ed_ded_function = {
    module R   : real
    module trm : trmodel with t = R.t

    val ed_ded_price_all [n][c][Ax][ns][nd] :
         trm.mp[n][c][Ax][ns][nd]
      -> (sa_max:i64)
      -> trm.prices[c][Ax]
      -> ([c][Ax-1]R.t, [c][Ax-1][c][Ax-1]R.t, [n]i64, [n]i64, [n]i64)
  }


module mk_run (E: ed_ded_function) = {

  type t = E.R.t

  module R = E.R
  module trm = E.trm
  module eqb = equilibrium R trm
  module lu = mk_lu R

  def ed_ded_price_all [n][c][Ax][ns][nd]
      (mp: trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) (p: trm.prices[c][Ax])
      : ([c][Ax-1]t, [c][Ax-1][c][Ax-1]t, [n]i64, [n]i64, [n]i64) =
    --- Using #[unsafe] to deal with a sink bug on the Futhark servers, which use an old Futhark version. 
    --- Remove this once futhark servers get Futhark version 0.25.3, as the code compiles perfectly well on my local machine with that newer version.
    #[unsafe] E.ed_ded_price_all mp sa_max p

  def newton_step [n][c][Ax][ns][nd]
      (mp: trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) (damp:t) (p: trm.prices[c][Ax]) (sa_iters_tot:[n]i64) (nk_iters_tot:[n]i64) 
      (rtrips_tot: [n]i64): (trm.prices[c][Ax], [c][Ax-1]t, t, [n]i64, [n]i64, [n]i64) =
    let (ed2d, ded4d, sa_iters, nk_iters, rtrips) : ([c][Ax-1]t, [c][Ax-1][c][Ax-1]t, [n]i64, [n]i64, [n]i64) =
      ed_ded_price_all mp sa_max p

    let sa_iters_tot = map2 (\x y -> x + y) sa_iters sa_iters_tot
    let nk_iters_tot = map2 (\x y -> x + y) nk_iters nk_iters_tot
    let rtrips_tot = map2 (\x y -> x + y) rtrips rtrips_tot

    let ed_flat  = flatten ed2d
    let ded_flat = map flatten (flatten ded4d)
    let blksz : i64 = 16
    let dp       = lu.ols blksz ded_flat ed_flat
    let dp2d : [c][Ax-1]t = unflatten dp

    -- Ensure that prices do not go below scrap value.
    let p_new : [c][Ax]t =
      tabulate_2d c Ax (\ct a ->
        if a == 0 then p[ct][0]
        else let a' = a - 1
             in R.(max (p[ct][a] - damp * dp2d[ct][a']) mp.pscrap[ct]))
    let max_ed = reduce R.max (R.i64 0) (map R.abs ed_flat)
    in (p_new, ed2d, max_ed, sa_iters_tot, nk_iters_tot, rtrips_tot)

  def newton [n][c][Ax][ns][nd] (mp: trm.mp[n][c][Ax][ns][nd]) (p0: trm.prices[c][Ax]) (sa_max:i64) (damp:t) (tol:t) (max_iter:i64)
      : (trm.prices[c][Ax], [c][Ax-1]t, t, i64, bool, [n]i64, [n]i64, [n]i64) =
    let ed0 : [c][Ax-1]t = replicate c (replicate (Ax-1) (R.i64 0))
    let start_iters = replicate n (0i64)
    let (p, ed, max_ed, iter, sa_iters_tot, nk_iters_tot, rtrips_tot) =
      loop (p, ed, max_ed, iter, sa_iters_tot, nk_iters_tot, rtrips_tot) = (p0, ed0, R.highest, 0i64, start_iters, start_iters, start_iters)
      while iter < max_iter && R.(max_ed > tol) do
        let (p', ed', max_ed', sa_iters_tot', nk_iters_tot', rtrips_tot') = newton_step mp sa_max damp p sa_iters_tot nk_iters_tot rtrips_tot
        in (p', ed', max_ed', iter + 1, sa_iters_tot', nk_iters_tot', rtrips_tot')
    in (p, ed, max_ed, iter, R.(max_ed <= tol), sa_iters_tot, nk_iters_tot, rtrips_tot)
}