"""g_path_trials.py: the trials ledger of the G prediction.

WHAT THIS MODULE COMPUTES.  The engine predicts M_Pl, and so G, by
inverting the electron's closure (root.py, _F_e_of):

    M_Pl = m_e / [C_EFF e^{-S} Delta_e^2 (1 + c2 alpha/2pi)
                  (1 - w4 (alpha/2pi)^2 Delta_e)],
    S = 2 pi / (b0 alpha_G2),  alpha_G2 = P h / 2pi = 1/(24 pi).

No continuous parameter is tuned on this path, and several factors
are chosen from finite menus. A path counts as a closure when the
closure is explained: why it holds, and why its alternatives fail on
grounds other than the target. This module places every choice on the
G path beside its menu. It solves for M_Pl again with each alternative,
the other choices held at their canonical values, and states the
explanation each choice has.

THE CHOICES (the explanation each has in the repository):
  theta = h(fund) = 2/9   EXPLAINED, blind: the positivity and
                          nondegeneracy gates leave it alone
                          (tests/probes/theta_menu.py)
  b0 = 32/3               EXPLAINED by matter content (d11 Dirac in
                          the 7, couplings.py); Weyl counting gives 38/3,
                          the labelled bridge light 19/3, light G2 scalars 21/2
                          (boundary_content.py B2)
  boundary M_b            IDENTIFIED (root.py): the scale where
                          alpha_G2 = 1/(24 pi) holds is taken to be the
                          full Planck mass. With G = 1/(c M_b^2), c is
                          the regulator t0 M_b^2 under another name, and
                          G is the only observable that depends on it.
                          The lane excludes readings that place the
                          regulator at a chosen scale. Of the factors
                          of general relativity, the gates as wired
                          pass k = 1 alone, with no measured G, and G
                          tests it (G10). That exclusion rests on the
                          wired threshold. Its independent ground is
                          the coupling variable, an identification
                          (tests/probes/coupling_normalization.py)
  loop order              EXPLAINED for G by the boundary principle
                          (the map is read whole at M_b). The datum
                          1/(24 pi) belongs to no perturbative scheme,
                          so it enters as the closure action
                          S0 = 2 pi/(b0 alpha). That is the only reading
                          with no scheme constant. Taken as a coupling
                          in an analytic scheme, it presupposes that
                          scheme and its result varies with the scheme
                          (M_Pl x0.22 MS-bar, x0.23 DR-bar, x0.30 V;
                          G6). The point mu* where the action matches
                          MS-bar is alpha_s's identification (IH):
                          k = (b1/2b0^2) ln 12 is the two loop Lambda
                          factor there (G8)
  alpha_G2 = P h/2pi      wiring.py's channel coupling definition. The
                          D^2_eff = 12 reading is the alternative
  C_EFF = 1/2             the closure's common mode taken whole:
                          C_EFF P Delta^2 P = 1 with P Delta^2 P =
                          1 + P v^2 P = 2, equal to the fundamental
                          channel's weight 1/d10. The menu
                          {4/9, 1/4, 2/9} needs Q = 3/4, 4/3, 3/2 and
                          fails the forgetting or the Koide bound
                          (tests/probes/frame_coefficients.py F1);
                          G_F tests it
  c2 = -1 in (1 + c2 alpha/2pi)
                          the bridge vertex that 1/alpha couples to,
                          -2 pi R(h_bridge). The response rule admits no
                          other weight of the web (frame_coefficients.py
                          F2); G_F tests it
  w4 = 27/2               the one channel share C_EFF of the unital
                          Albert traversal 27 (frame_coefficients.py F3,
                          echo_naturality.py); m_mu tests it on the
                          29-entry menu (vertex_pinning.py)

STATEMENTS (the numbers below check each one).
  G1  Every choice except the depth-2 coefficient and the depth-4
      weight moves M_Pl by an O(1) factor, so its menu supports the
      0.1% statement: given the menus, M_Pl from the four integers
      comes out within -0.11% at depth 1.
  G2  The step from -1140 ppm to +21.5 ppm is the depth-2 coefficient
      alone. Its menu is spaced by alpha/2pi = 1161 ppm in M_Pl
      (581 ppm for half integers). For an unexplained coefficient, a
      match within the achieved 22 ppm has a 3.8% (7.7%) chance.
  G3  The depth-4 weight moves M_Pl by 0.7 ppm. Across its whole
      menu it spans 0 to 26 ppm.
  G4  The same coefficient enters the electron in the ratio v_EW/m_e,
      where M_Pl cancels: v_EW = M_Pl e^{-S_quark} with
      S_quark = S0 - C2(26) + N (pi/512)/2pi - w16 (alpha/2pi)^2. G_F
      fixes v_EW to 0.26 ppm. The c2 menu is the rival readings of
      frame_coefficients.py F2, c2 = -2 pi R(h) = -h over the conformal
      weights of the web: h_bridge = 1, h_7 = 2/5, h_26 = 3/5, h11 = 1/2
      and h10 = 2/9. With N = 30 held, over those five readings x the
      29-entry w16 menu, 145 combinations, exactly one lies within
      2 sigma: (c2, N, w16) = (-1, 30, 16), at +0.019 ppm. The next is
      2.7 ppm away (11 sigma), and 465 ppm among the readings at
      w16 = 16. The vertex unit pi/512 and the depth-2 unit alpha(0)
      are incommensurate, so no integer trade between N and c2 exists.
  G5  The chance that a random G_F target finds a combination within
      2 sigma is 0.22% with N = 30 and w16 = 16 taken as derived
      (N = Tr+(Q^2) = h-dual(E8), w16 = charge trace x C2(26)), and
      4.0% with N alone. On the derived N and w16, c2 = -1 is fixed
      by G_F, an observable independent of G. The coefficient is
      therefore explained by the electroweak closure, and the -55 ppm
      on G is a prediction.
  G6  C_EFF enters v_EW/m_e as 1/C_EFF, so G_F fixes it too: every
      menu alternative moves v_EW by 11% or more. The exponent
      S0 = 2 pi/(b0 alpha_G2) cancels from v_EW/m_e, and mu* = v_EW
      e^(S_quark - S0 + C2(26)) carries no S0 either, so alpha_s does
      not test it. G is the only precision test of the absolute
      exponent: its -55 ppm is S0 at 6.1e-7. Of the choices S0 rests
      on, b0 = 32/3 (matter content) and alpha_G2 = P h/2pi (wiring)
      are structural. The one loop truncation treats the exponent as
      an exact RG invariant. Its alternative moves M_Pl by a factor
      near 5. That alternative is not a single reading. It treats the
      boundary datum as a coupling in an analytic scheme, and each
      scheme gives its own M_Pl (x0.218 MS-bar, x0.232 DR-bar, x0.303
      V scheme), because the two loop relation between coupling and
      scale carries the scheme's constant. The boundary principle
      reads the datum in no scheme ("at the boundary the map is read
      whole", boundary.py). Every analytic reading is therefore
      excluded at once. The closure action, alpha = 2 pi/(b0 S), is
      the only reading that needs no scheme constant. For G the loop
      order is explained without the measured G, and G tests it.
  G7  The same one loop exactness fixes alpha_s. One loop running from
      alpha_G2(M_Pl) = 1/(24 pi) over ln(M_Pl/mu*) = S0 - C2(26) gives
      1/alpha(mu*) = 24 pi - (b0/2pi)(S0 - 6) = 32/pi exactly, the
      couplings.py boundary pi/32 = 2pi/(b0 C2(26)). Two loop G2
      running over the same interval gives alpha_G2(mu*) = 0.1112 and
      alpha_s(M_Z) = 0.1374, 21.5 sigma from PDG. One loop gives
      0.11839 (+0.4 sigma). So 1/(24 pi) is not the MS-bar coupling at
      M_Pl. The G2 coupling is the closure action, alpha = 2pi/(b0 S).
      Actions add (S(mu*) = C2(26)), and G tests the exponent at 5e-7.
  G8  The action reading meets MS-bar at mu* only by identification.
      With the standard two loop MS-bar Lambda taken as the action's
      Lambda, alpha_MS(mu*) = 0.0855 and alpha_s(M_Z) = 0.1007
      (-19.2 sigma). Scheme constants of +-0.5 give -8.4 and
      -28.3 sigma. The engine's alpha_s holds for alpha_MS(mu*) = pi/32
      exactly, a scheme constant k = -(b1/2b0^2) ln(b0/128) = +0.8445,
      Lambda_action = 2.327 Lambda_MSbar. Since b0/128 = 1/12 is
      b0 alpha/4pi at alpha = pi/32, k is the two loop Lambda factor
      (b0 alpha/4pi)^(-b1/2b0^2) at mu*, not a finite conversion. The
      action and MS-bar couplings coincide at mu*, and the G6
      alternative is the same factor evaluated at M_b. The loop order
      and k are one identification. The alpha_s agreement is
      therefore conditional on that identification at the matching
      scale, which the engine does not derive. A derivation of the
      action to MS-bar conversion that gives any other constant kills
      the alpha_s reading.
  G9  The band is narrow: alpha_s within 2 sigma allows
      alpha_MS(mu*) between -1.57% and +1.01% of pi/32. Read as the
      static potential (V) scheme coupling, a1 = (31/9)C_A -
      (20/9)T_F n_f = 64/9 for G2 with three Dirac 7s, pi/32 gives
      alpha_s = 0.1111 (-7.7 sigma). Read in MS without the 4pi e^-gamma
      rescaling, 0.1465 (+31.6 sigma). The identification is
      therefore specific: the closure action's coupling at the
      matching scale is the MS-bar coupling, to about one percent.
  G10 The boundary. The inversion computes M_b, the scale at which
      alpha_G2 = 1/(24 pi) holds, and reads M_b as the full Planck
      mass sqrt(hbar c/G). Read as the reduced mass sqrt(hbar c/8 pi G)
      it moves M_Pl by sqrt(8 pi) = 5.01. As the regulator scale of the
      cutoff identification t0 = d11 l_Pl^2, by sqrt(d11) = 1.73. With
      1/G = Sigma/(2 pi t0) the full Planck reading is
      t0 = Sigma/(2 pi) M_b^-2 = 4.112 l_Pl^2 on the engine's one list,
      4.297 on the full ledger (boundary_content.py B4). The
      regulator is the boundary written in proper time form, and the
      cutoff menu of sakharov_normalisation.py FACT 3 is its menu.
      Write G = 1/(c M_b^2). The anchored regulator is the root of
      Sigma(t0) = 2 pi c t0, so choosing c is choosing t0 M_b^2. c is
      the regulator under another name, and no source derives the
      regulator (sakharov_normalisation.py FACT 6, 7). Every other
      dimensionful output of the engine is computed from M_b, so c
      affects G alone (and rho_crit, which is G times H0^2).
      Gate S, with no measured G, asks t0 >= M_b^-2, so c <= 4.11 /
      6.98 / 9.84 at contacts +1 / 0 / -1.
      Two readings place the regulator: the cutoff scale, c = d11
      (t0 = d11 l_Pl^2, FACT 7's diagnostic), and the regulator at the
      boundary, c = Sigma/(2 pi) (t0 = M_b^-2, on S's edge, M_Pl x2.03).
      Each obtains G from the ledger with the regulator placed at a
      chosen scale, which the lane excludes. The regulator's
      normalization is scheme dependent (FACT 6, 7), so identifying
      the regulator with a chosen scale has no scheme free content.
      The remaining readings relate M_b to G itself,
      M_b^2 = 1/(k G), with k a factor with which G enters general
      relativity: 1 (force law), 4 pi (Gauss), 8 pi (field equations),
      16 pi (action), 32 pi (graviton coupling).
      S's threshold is soft, so it is evaluated under each convention:
      the boundary at t0 = 0.267, 1 or 3.74 M_b^-2 (a row at M_b at
      half weight, the nominal point, its mirror), with a sharp or a
      smooth regulator (FACT 7's 2 ln 2). k = 1 passes under every
      convention; 8 pi, 16 pi and 32 pi fail under every convention
      (at the nominal threshold 8 pi would need 758 to 974 bridge
      units across the contacts, against the 65-376 window); 4 pi
      passes under two of six. At the wired threshold, the nominal
      point at every contact (gravity.py gate S), 4 pi fails at each
      contact and 2 pi at the covariant one, and gate M, the anchored
      regulator within 0.6 to 15 M_b^-2, fails 4 pi at two contacts.
      Off the menu, c = 2 or 3 passes under four of six and at the
      wired threshold. On the menu the gates as wired leave k = 1
      alone, with no measured G, and G tests it. That exclusion of
      2 pi and 4 pi rests on the wired threshold. Its independent
      ground is the coupling variable, an identification. 1/alpha(0)
      matches the force law coefficient, a factor 4 pi from either
      rationalized one, and a single definition gives the common mode
      in the same variable (tests/probes/coupling_normalization.py).
      The menu itself, the factors of general relativity, is a
      judgment. The wired gates hold c between 0.656 and 4.11.
      Why c stays an identification. c has one observable in nature as
      well as in the engine: general relativity has one coupling
      besides Lambda, so every gravitational observable measures G.
      The ledger fixes c only to order one (S, M), and the lane forbids
      sharpening it through the regulator. No present premise can fix
      c, and no observable except G can test it. So c is an
      identification by necessity, and G is its test.
      The parameter free values that give the 1 have rivals: unity in
      the engine's coupling variable against the amount's 1/(2 pi),
      and the q = 0 second moment Pv2P = 1 against the per mode 1/2.
      Each rivalry sets a sum against one of its shares. The common
      mode carries no label that could pick out a share, so it takes
      the sum. On the menu of the engine's own shares, c = 1 alone
      passes S under every convention (tests/probes/boundary_whole.py).
      On the whole lattice of those shares, 216 readings each with c
      at least one, it is the only fusion closed set of closures that
      holds every closure (tests/probes/boundary_lattice.py). That is
      the common mode at depth one. On the depth four web it is the
      count of the lepton circulant (tests/probes/common_mode_count.py):
      G M_b^2 = 3 Q_l - 1, c = 1/(3 Q_l - 1). This is the
      identification that G tests (root.py, THE COMMON MODE).

Usage: python3 tests/audit/g_path_trials.py
"""

import math
import os
import sys
from fractions import Fraction

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_here), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402
import vertex_pinning as vp  # noqa: E402

CODATA = r.M_PL_CODATA_GEV
SIGMA_G = 22e-6


def m_pl(c_eff=None, S=None, c2=-1.0, w4=None):
    a = 1.0 / r.inv_alpha_phys
    e = a / (2 * math.pi)
    De = r._Delta_e
    c_eff = r.C_EFF if c_eff is None else c_eff
    S = r.S_lepton if S is None else S
    w4 = float(Fraction(r.d11**3, r.d10)) if w4 is None else w4
    return r.M_E_ANCHOR_MEV / (c_eff * math.exp(-S) * De**2 * (1 + c2 * e)
                               * (1 - w4 * e * e * De)) / 1e3


def ppm(M):
    return 1e6 * (M / CODATA - 1)


# The rival readings of c2 (frame_coefficients.py F2): c2 = -2 pi R(h) = -h
# over the conformal weights of the web. The response rule admits h_bridge alone.
C2_READINGS = {"h_bridge": -float(r.h_bridge), "h_7": -float(r.h_7), "h_26": -float(r.h_26),
               "h11": -float(r.h11), "h10": -float(r.h10)}
C2_MENU = tuple(C2_READINGS.values())
N_TRACE = 30                          # N = Tr+(Q^2) = h-dual(E8), held as derived


def gf_lock(win_sigma=2.0, span=1500.0, samples=20000, seed=1):
    """v_EW/m_e against G_F over (c2, w16) at N = 30, and the chance rate of a match."""
    import random
    from interference.reference import PDG_EW
    a = 1.0 / r.inv_alpha_phys
    e = a / (2 * math.pi)
    De, w4 = r._Delta_e, r.d11**3 / r.d10
    v_meas = (math.sqrt(2) * PDG_EW["G_F"])**-0.5
    sv = 0.5 * PDG_EW["G_F_err"] / PDG_EW["G_F"]
    S0, C2 = r.S_lepton, float(r.C2_26)

    def res(c2, N, w16):
        Fe = r.C_EFF * math.exp(-S0) * De**2 * (1 + c2 * e) * (1 - w4 * e * e * De)
        Sq = S0 - C2 + N * (math.pi / 512) / (2 * math.pi) - w16 * e * e
        return 1e6 * (r.M_E_ANCHOR_MEV / 1e3 * math.exp(-Sq) / Fe / v_meas - 1)
    win = win_sigma * 1e6 * sv
    w16_menu = [float(w) for w in vp.menu()]
    N = N_TRACE
    full = {(c2, N, w): res(c2, N, w) for c2 in C2_MENU for w in w16_menu}
    hits = [(k, round(v, 3)) for k, v in full.items() if abs(v) < win]
    runner = sorted(abs(v) for v in full.values())[1]
    scopes = {"N, w16 derived": [res(c2, N, 16.0) for c2 in C2_MENU],
              "N derived": list(full.values())}
    rng = random.Random(seed)
    targets = [rng.uniform(-span, span) for _ in range(samples)]
    chance = {}
    for k, vals in scopes.items():
        vs = sorted(vals)
        import bisect
        n = 0
        for t in targets:
            i = bisect.bisect_left(vs, t - win)
            n += i < len(vs) and vs[i] < t + win
        chance[k] = n / samples
    assert hits == [((-1.0, 30, 16.0), 0.019)] and runner > 2.5
    assert len(full) == 145 and chance["N, w16 derived"] < 0.005 < chance["N derived"] < 0.05
    return dict(n_combos=len(full), hits=hits, runner_up_ppm=runner, chance=chance)


def loop_order():
    """alpha_s(M_Z) with one loop and two loop G2 running from M_Pl to mu*."""
    import contextlib
    import io
    from interference import couplings as cpl
    from interference import masses
    with contextlib.redirect_stdout(io.StringIO()):
        m = masses.derive(r.derive())
    b0, b1 = 32 / 3, 232 / 3
    a_pl, M = 1 / (24 * math.pi), r.M_Pl_GeV
    mu_star = M * math.exp(-(r.S_lepton - float(r.C2_26)))
    L = math.log(M / mu_star)

    def down(a, two):
        n, h = 20000, -L / 20000
        f = lambda x: -(b0 / (2 * math.pi)) * x * x - two * (b1 / (8 * math.pi**2)) * x**3
        for _ in range(n):
            k1 = f(a); k2 = f(a + h * k1 / 2); k3 = f(a + h * k2 / 2); k4 = f(a + h * k3)
            a += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return a
    v0, mt = m["v_EW_pred_GeV"], m["m_t"] / 1e3

    def alpha_MZ(a_star):
        MV = math.sqrt(4 * math.pi * a_star) * v0 * math.sqrt(float(r.E_v2) / 3)
        lam3 = 1 - r.d11 * r.n7 * math.log(MV / mu_star)
        a = 1 / (1 / a_star - lam3 / (12 * math.pi))
        return cpl._run_SM_2loop(cpl._run_SM_2loop(a, mu_star, mt, 6), mt, 91.1876, 5)
    one, two = down(a_pl, 0), down(a_pl, 1)
    out = dict(inv_alpha_star_one=1 / one, alpha_star_two=two,
               alpha_MZ_one=alpha_MZ(one), alpha_MZ_two=alpha_MZ(two), alpha_MZ_of=alpha_MZ)
    assert abs(1 / one - 32 / math.pi) < 1e-6
    assert abs(out["alpha_MZ_one"] - 0.11839) < 5e-5 and out["alpha_MZ_two"] > 0.135
    return out


def scheme_conversion():
    """alpha_s(M_Z) when the action scheme Lambda is read through the MS-bar two loop relation."""
    from scipy.optimize import brentq
    lo = loop_order()
    b0, b1 = 32 / 3, 232 / 3
    k_needed = -(b1 / (2 * b0**2)) * math.log(b0 / 128)
    out = {"k_needed": k_needed, "Lambda_ratio": math.exp(k_needed)}
    alpha_MZ = lo["alpha_MZ_of"]
    for k in (0.0, -0.5, 0.5):
        f = lambda a: 2 * math.pi / (b0 * a) + (b1 / (2 * b0**2)) * math.log(b0 * a / (4 * math.pi)) + k - 6
        a_ms = brentq(f, 1e-3, 2.0)
        out[k] = (a_ms, alpha_MZ(a_ms))
    assert abs(out[0.0][1] - 0.1007) < 3e-4 and abs(k_needed - 0.8445) < 1e-3
    return out


def scheme_band():
    """The alpha_MS(mu*) band alpha_s allows, and pi/32 read in standard schemes."""
    from scipy.optimize import brentq
    lo_ = loop_order()
    f = lo_["alpha_MZ_of"]
    a_star = math.pi / 32
    lo = brentq(lambda a: f(a) - (0.1180 - 0.0018), 0.05, 0.12)
    hi = brentq(lambda a: f(a) - (0.1180 + 0.0018), 0.05, 0.12)
    b0, gamma = 32 / 3, 0.5772156649
    a1V = (31 / 9) * 4 - (20 / 9) * 1 * 3
    reads = {"V": 1 / (1 / a_star + a1V / (4 * math.pi)),
             "MS": 1 / (1 / a_star - (b0 / (4 * math.pi)) * (math.log(4 * math.pi) - gamma))}
    out = dict(band=(lo / a_star - 1, hi / a_star - 1), a1V=a1V,
               alpha_s={k: f(a) for k, a in reads.items()})
    assert abs(a1V - 64 / 9) < 1e-12
    assert out["alpha_s"]["V"] < 0.1180 - 0.005 and out["alpha_s"]["MS"] > 0.1180 + 0.02
    return out


def datum_readings():
    """G6: the boundary datum read as the closure action or as a coupling of an analytic scheme.

    Read as alpha_X(M_b) with 1/alpha_MS = 1/alpha_X + delta_X and the
    two loop Lambda, each scheme gives its own S. The action reading
    S0 = 2 pi/(b0 alpha) carries no scheme constant.
    """
    b0, b1 = 32 / 3, 232 / 3
    a = 1 / (24 * math.pi)
    S0 = 2 * math.pi / (b0 * a)

    def analytic(delta):
        a_ms = 1 / (1 / a + delta)
        return 2 * math.pi / (b0 * a_ms) + (b1 / (2 * b0**2)) * math.log(b0 * a_ms / (4 * math.pi))
    S = {"action": S0, "MS-bar": analytic(0.0), "DR-bar": analytic(4 / (12 * math.pi)),
         "V": analytic((64 / 9) / (4 * math.pi))}
    base = m_pl()
    shift = {k: m_pl(S=v) / base for k, v in S.items()}
    assert abs(shift["action"] - 1.0) < 1e-9
    assert abs(shift["MS-bar"] - 0.2177) < 1e-3 and abs(shift["DR-bar"] - 0.2316) < 1e-3
    assert abs(shift["V"] - 0.3030) < 1e-3
    return dict(S=S, mpl_shift=shift)


def boundary_gate():
    """G10: the boundary menu tested by gate S on the one list, with no measured G.

    G = 1/(c M_b^2). The anchored regulator solves Sigma(t0) = 2 pi c t0
    with t0 in M_b^-2, and gate S (boundary.inside) asks t0 >= M_b^-2,
    the same test for every reading. The regulator placed at the
    boundary is c = Sigma(M_b^-2)/(2 pi), contact by contact.
    """
    from scipy.optimize import brentq
    from interference import boundary as B
    one, xi = B.rows(r.d10, r.d11, r.n7, r.n26), float(r.xi_bridge)
    sigma_at_Mb = {ct: B.sigma(one, 1.0, ct, xi) for ct in B.CONTACTS}
    c_max = {ct: s / (2 * math.pi) for ct, s in sigma_at_Mb.items()}
    cs = {"full": dict.fromkeys(B.CONTACTS, 1.0),
          "reduced": dict.fromkeys(B.CONTACTS, 8 * math.pi),
          "cutoff": dict.fromkeys(B.CONTACTS, float(r.d11)),
          "regulator at M_b": c_max}
    t0 = {k: {ct: brentq(lambda t: B.sigma(one, t, ct, xi) - 2 * math.pi * c[ct] * t, 1e-9, 100.0)
              for ct in B.CONTACTS} for k, c in cs.items()}
    survives = {k: all(v >= 1 - 1e-9 for v in t0[k].values()) for k in cs}
    # "regulator at M_b" passes by construction (c = Sigma(1)/2 pi puts t0 at 1): it marks S's edge
    assert all(B.inside(one, ct, xi) for ct in B.CONTACTS)
    assert all(abs(t0["full"][ct] - B.anchored(one, ct, xi)[0] / (2 * math.pi)) < 1e-9
               for ct in B.CONTACTS)
    assert survives == {"full": True, "reduced": False, "cutoff": True, "regulator at M_b": True}
    # the menu of regulator free readings, M_b^2 = 1/(k G), read across S's conventions
    from scipy.special import expn
    x_half = brentq(lambda x: expn(2, x) - 0.5, 1e-3, 5.0)        # a row at M_b at half weight
    conventions = [(th, F) for th in (x_half, 1.0, 1.0 / x_half) for F in (1.0, 2 * math.log(2))]
    gr = {"1": 1.0, "4pi": 4 * math.pi, "8pi": 8 * math.pi, "16pi": 16 * math.pi, "32pi": 32 * math.pi,
          "2 (off menu)": 2.0, "3 (off menu)": 3.0}
    # the anchored regulator solves F Sigma(t0) = 2 pi k t0; Sigma falls with t0, so it lies at or
    # beyond the threshold th exactly when F Sigma(th) >= 2 pi k th (the argument of B.inside)
    gr_menu = {name: sum(all(B.sigma(one, th, ct, xi) * F >= 2 * math.pi * k * th for ct in B.CONTACTS)
                         for th, F in conventions) for name, k in gr.items()}
    assert gr_menu["1"] == len(conventions)
    assert gr_menu["8pi"] == gr_menu["16pi"] == gr_menu["32pi"] == 0
    assert 0 < gr_menu["4pi"] < len(conventions)
    # the gates as wired: S at the nominal threshold and M's band, contact by contact
    gr["2pi"] = 2 * math.pi
    t0_gr = {name: [brentq(lambda t: B.sigma(one, t, ct, xi) - 2 * math.pi * k * t, 1e-9, 1e3)
                    for ct in B.CONTACTS] for name, k in gr.items()}
    wired = {name: (all(t >= 1 - 1e-9 for t in ts), all(0.6 < t < 15.0 for t in ts))
             for name, ts in t0_gr.items()}
    assert wired["1"] == (True, True) and wired["2pi"] == (False, True)
    assert wired["4pi"] == wired["8pi"] == (False, False)
    assert [name for name in ("1", "2pi", "4pi", "8pi", "16pi", "32pi") if all(wired[name])] == ["1"]
    c_window = (max(sigma_at_Mb.values()) / (2 * math.pi * 15.0), min(c_max.values()))
    return dict(t0=t0, c_max=c_max, sigma_at_Mb=sigma_at_Mb, survives=survives,
                mpl_shift={k: math.sqrt(c[1.0]) for k, c in cs.items()},
                gr_menu=gr_menu, n_conventions=len(conventions), x_half=x_half,
                wired=wired, c_window=c_window)


def run(report=print):
    alpha_G2 = 1 / (24 * math.pi)
    S1 = 2 * math.pi / (Fraction(32, 3) * alpha_G2)
    assert abs(S1 - r.S_lepton) < 1e-12
    b0, b1 = 32 / 3, 232 / 3
    two_loop = (b1 / (2 * b0**2)) * math.log(b0 * alpha_G2 / (4 * math.pi))
    menus = {
        "C_EFF": {"4/9": 4 / 9, "1/2": 0.5, "1/4": 0.25, "2/9": 2 / 9},
        "b0": {"32/3 (Dirac)": 32 / 3, "38/3 (Weyl)": 38 / 3,
               "19/3 (labelled light bridge)": 19 / 3, "21/2 (light G2 scalars)": 21 / 2},
        "boundary": {"full": 1.0, "reduced": math.sqrt(8 * math.pi), "cutoff": math.sqrt(r.d11)},
        "loop": {"one": 0.0, "two": two_loop},
        "alpha_G2": {"P h/2pi = 1/(24pi)": 1.0, "D^2_eff = 12": 3.0},
        "c2": {str(c): c for c in (0.0, -0.5, -1.0, -1.5, -2.0, 1.0)},
        "w4": {str(v): float(v) for v in vp.menu()},
    }
    out = {}
    out["C_EFF"] = {k: ppm(m_pl(c_eff=v)) for k, v in menus["C_EFF"].items()}
    out["b0"] = {k: ppm(m_pl(S=2 * math.pi / (v * alpha_G2))) for k, v in menus["b0"].items()}
    out["loop"] = {k: ppm(m_pl(S=r.S_lepton + v)) for k, v in menus["loop"].items()}
    out["alpha_G2"] = {k: ppm(m_pl(S=r.S_lepton / v)) for k, v in menus["alpha_G2"].items()}
    out["boundary"] = {k: ppm(v * m_pl()) for k, v in menus["boundary"].items()}
    out["c2"] = {k: ppm(m_pl(c2=v)) for k, v in menus["c2"].items()}
    out["w4"] = {k: ppm(m_pl(w4=v)) for k, v in menus["w4"].items()}
    depth1 = ppm(m_pl(c2=0.0, w4=0.0))
    e = 1.0 / (2 * math.pi * r.inv_alpha_phys)
    spacing = 1e6 * e
    achieved = abs(out["c2"]["-1.0"])
    luck_int, luck_half = 2 * achieved / spacing, 2 * achieved / (spacing / 2)
    w4_span = max(out["w4"].values()) - min(out["w4"].values())
    assert abs(depth1 + 1139.9) < 0.5 and abs(out["c2"]["-1.0"] - 22.2) < 0.1
    assert abs(out["w4"]["27/2"] - out["w4"]["0"] - 0.73) < 0.02
    assert min(abs(v) for k, v in out["C_EFF"].items() if k != "1/2") > 1e5
    assert abs(out["boundary"]["full"] - out["c2"]["-1.0"]) < 1e-9
    assert min(abs(v) for k, v in out["boundary"].items() if k != "full") > 7e5
    g_ppm = 1e6 * (r.COMMON_MODE * (1 + 1e-6 * out['c2']['-1.0'])**-2 - 1)    # G read through the common mode
    out_g = dict(G_ppm=g_ppm, G_sigma=1e-6 * g_ppm / SIGMA_G)
    assert abs(g_ppm + 54.6) < 0.1
    report(f"canon: M_Pl {out['c2']['-1.0']:+.1f} ppm, G {g_ppm:+.1f} ppm "
           f"({1e-6 * g_ppm / SIGMA_G:+.2f} sigma); M_Pl through depth 1 alone {depth1:+.1f} ppm")
    for name, res in out.items():
        if name == "w4":
            continue
        report(f"  {name:9s}: " + ", ".join(f"{k} {v:+.3g}" for k, v in res.items()))
    report(f"  w4 menu  : 29 entries span {w4_span:.1f} ppm in M_Pl; 27/2 sits at "
           f"{out['w4']['27/2'] - out['w4']['0']:+.2f} ppm from 'no term'")
    report(f"depth-2 coefficient: menu spacing {spacing:.0f} ppm (integer), {spacing / 2:.0f} ppm "
           f"(half integer). Chance of a match within {achieved:.1f} ppm: "
           f"{100 * luck_int:.1f}% / {100 * luck_half:.1f}%")
    gf = gf_lock()
    report(f"G_F lock: {gf['n_combos']} combinations of (c2, w16) at N = 30. Within 2 sigma: "
           f"{gf['hits']}. Next nearest {gf['runner_up_ppm']:.0f} ppm")
    report("  chance of a random 2-sigma match: " + ", ".join(
        f"{k} {100 * v:.2f}%" for k, v in gf['chance'].items()))
    lo = loop_order()
    report(f"loop order: 1/alpha_G2(mu*) one loop {lo['inv_alpha_star_one']:.6f} (32/pi); "
           f"alpha_s(M_Z) one loop {lo['alpha_MZ_one']:.5f}, two loop {lo['alpha_MZ_two']:.5f} "
           f"({(lo['alpha_MZ_two'] - 0.1180) / 0.0009:+.1f} sigma)")
    sc = scheme_conversion()
    report(f"scheme: standard MS-bar Lambda gives alpha_MS(mu*) = {sc[0.0][0]:.4f}, alpha_s(M_Z) = "
           f"{sc[0.0][1]:.4f} ({(sc[0.0][1] - 0.1180) / 0.0009:+.1f} sigma). The engine's reading needs "
           f"k = {sc['k_needed']:+.4f}, Lambda_action/Lambda_MSbar = {sc['Lambda_ratio']:.3f}")
    sb = scheme_band()
    report(f"scheme band: alpha_MS(mu*) within {100 * sb['band'][0]:+.2f}% / {100 * sb['band'][1]:+.2f}% "
           f"of pi/32. Pi/32 as V-scheme gives alpha_s {sb['alpha_s']['V']:.4f}, as MS {sb['alpha_s']['MS']:.4f}")
    dr = datum_readings()
    report("boundary datum read as the closure action or as an analytic scheme coupling (M_Pl / canon): "
           + ", ".join(f"{k} {v:.3f}" for k, v in dr["mpl_shift"].items())
           + ". The boundary principle reads it in no scheme, so the action is the only reading left")
    bg = boundary_gate()
    fmt = lambda d: " / ".join(f"{v:.3f}" for v in d.values())
    report(f"boundary under gate S: c_max {fmt(bg['c_max'])} at contacts +1/0/-1. Anchored t0 "
           f"(M_b^-2) full {fmt(bg['t0']['full'])}, reduced {fmt(bg['t0']['reduced'])} (outside), "
           f"cutoff {fmt(bg['t0']['cutoff'])}. Pass S: "
           + ", ".join(k for k, v in bg["survives"].items() if v) + ". The lane rule removes readings that place the regulator; "
           f"GR menu under {bg['n_conventions']} conventions of S: "
           + ", ".join(f"k = {n} {v}/{bg['n_conventions']}" for n, v in bg["gr_menu"].items())
           + ". At the wired threshold S and M pass "
           + ", ".join(f"k = {n}" for n, v in bg["wired"].items() if all(v))
           + f"; c window {bg['c_window'][0]:.3f} to {bg['c_window'][1]:.2f}; G tests k = 1")
    report("RESULT: the depth-2 coefficient is fixed by G_F on the derived N = 30 and w16 = 16, "
           "so the ppm statement on G is a prediction. The 0.1% statement rests on O(1)-spaced "
           "menus: C_EFF is fixed by G_F (G6). The wired gates leave k = 1 alone on the GR menu "
           "without G, and G tests it (G10). The loop order takes the closure action whole. The point "
           "where it matches MS-bar is alpha_s's identification (G6, G8).")
    return dict(menus=out, depth1=depth1, G=out_g, spacing_ppm=spacing, luck_int=luck_int,
                luck_half=luck_half, w4_span=w4_span, gf=gf, loop=lo, scheme=sc, band=sb,
                boundary_gate=bg, datum_readings=dr)


if __name__ == "__main__":
    run()
