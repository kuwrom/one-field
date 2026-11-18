"""g_lane_exhaustion.py: every lane of the G path, crossed, against every observable that reads it.

WHAT THIS MODULE COMPUTES.  g_path_trials.py varies one menu at a time
with the others held at their canonical values. Its G_F lock crosses
(c2, N, w16). A closure claim needs more than that: the canonical lane
must be the only lane in the full product of the menus that closes. This module takes the full product of every menu on the G path.
It counts the lanes that close on each observable, and on all of them
together.

THE MENUS (each as g_path_trials.py states it):
  C_EFF      4/9, 1/2, 1/4, 2/9
  b0         32/3 (Dirac), 38/3 (Weyl), 19/3 (labelled light bridge),
             21/2 (light G2 scalars)
  alpha_G2   P h/2pi = 1/(24 pi), the D^2_eff = 12 reading (x3)
  reading    the boundary datum as the closure action, or as an MS-bar,
             DR-bar or V-scheme coupling with the two loop Lambda (G6;
             b1 = 232/3 held for every b0)
  c          G = 1/(c M_b^2): 1, 4 pi, 8 pi, 16 pi, 32 pi (the factors
             of general relativity), d11 (cutoff scale), Sigma/2 pi
             (regulator at the boundary), 2
  c2         the depth 2 coefficient, the rival readings of
             frame_coefficients.py F2: c2 = -2 pi R(h) = -h over the
             conformal weights h_bridge, h_7, h_26, h11, h10
  w4, w16    the depth 4 weights, the 29-entry vertex menu each
N = 30 (h-dual(E8)) and theta = 2/9 are held: each is explained
(g_path_trials.py G5, theta_menu.py).

THE OBSERVABLES AND WHAT EACH READS:
  G     every menu except w16, through the common mode:
        G = common_mode/(c M_b^2). Window 3 sigma. The canonical lane
        sits at -2.43 sigma.
  G_F   v_EW/m_e. It depends on C_EFF, c2, w4 and w16, and not on the
        exponent, the reading or c. Window 2 sigma. Every lane keeps the
        engine's form S_quark = S - C2(26) + N (pi/512)/2pi
        - w16 (alpha/2pi)^2, so S cancels. Every alternative is read with
        that form.
  m_mu  the gap ratio and the depth 4 vent, which depends on w4 alone.
        Window 2 sigma (vertex_pinning.py).

STATEMENTS (the numbers below check each one).
  E1  G alone does not select a lane. Its fine coefficients (c2, w4)
      can pull a coarse lane that comes out near the target inside
      the window.
  E2  G_F and m_mu pin C_EFF, c2, w4 and w16 without G.
  E3  With those pinned, exactly one lane of the full product closes on
      all three observables: the canonical one. The nearest other lane
      that closes on G_F and m_mu misses G by the margin printed.

Usage: python3 tests/audit/g_lane_exhaustion.py
"""

import itertools
import math
import os
import sys
from fractions import Fraction

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_here), "probes"))
sys.path.insert(0, os.path.dirname(os.path.dirname(_here)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402
from interference import boundary as B  # noqa: E402
import vertex_pinning as vp  # noqa: E402
sys.path.insert(0, _here)
from g_path_trials import C2_READINGS  # noqa: E402

CODATA_MPL = r.M_PL_CODATA_GEV
SIGMA_G = r.G_REL_SIGMA
B0, B1 = 32 / 3, 232 / 3


def _menus():
    one = B.rows(r.d10, r.d11, r.n7, r.n26)
    sigma_cov = B.sigma(one, 1.0, 1.0, float(r.xi_bridge))
    return {
        "C_EFF": {"1/2": 0.5, "4/9": 4 / 9, "1/4": 0.25, "2/9": 2 / 9},
        "b0": {"32/3": 32 / 3, "38/3": 38 / 3, "19/3": 19 / 3, "21/2": 21 / 2},
        "alpha_G2": {"P h/2pi": 1 / (24 * math.pi), "D2_eff=12": 3 / (24 * math.pi)},
        "reading": {"action": None, "MS-bar": 0.0, "DR-bar": 4 / (12 * math.pi), "V": (64 / 9) / (4 * math.pi)},
        "c": {"1": 1.0, "4pi": 4 * math.pi, "8pi": 8 * math.pi, "16pi": 16 * math.pi, "32pi": 32 * math.pi,
              "d11": float(r.d11), "Sigma/2pi": sigma_cov / (2 * math.pi), "2": 2.0},
        "c2": dict(C2_READINGS),
        "w4": {str(v): float(v) for v in vp.menu()},
        "w16": {str(v): float(v) for v in vp.menu()},
    }


def _exponent(b0, a, delta):
    """S for the boundary datum a read as the action (delta None) or as an analytic scheme coupling."""
    if delta is None:
        return 2 * math.pi / (b0 * a)
    a_ms = 1 / (1 / a + delta)
    return 2 * math.pi / (b0 * a_ms) + (B1 / (2 * b0**2)) * math.log(b0 * a_ms / (4 * math.pi))


def run(report=print):
    from interference.reference import PDG_EW
    M = _menus()
    a = 1.0 / r.inv_alpha_phys
    e = a / (2 * math.pi)
    De = r._Delta_e
    me_GeV = r.M_E_ANCHOR_MEV / 1e3
    C2, N = float(r.C2_26), 30
    v_meas = (math.sqrt(2) * PDG_EW["G_F"])**-0.5
    sv = 0.5 * PDG_EW["G_F_err"] / PDG_EW["G_F"]

    # m_mu depends on w4 alone (vertex_pinning node 3)
    gap = {k: 1.0 + r.BA_ratio * math.cos(float(r.h10) + 2 * math.pi * k / 3) for k in range(3)}
    ref_mu, sig_mu = vp.M_MU_CODATA

    def mu_pull(w4):
        vent = lambda g: 1.0 - w4 * e * e * g
        return (r.M_E_ANCHOR_MEV * (gap[2] / gap[1])**2 * vent(gap[2]) / vent(gap[1]) - ref_mu) / sig_mu
    mu_ok = {k: abs(mu_pull(w)) < 2.0 for k, w in M["w4"].items()}

    # G_F depends on (C_EFF, c2, w4, w16). The exponent and c cancel from v_EW/m_e
    def gf_pull(ce, c2, w4, w16):
        Fe = ce * De**2 * (1 + c2 * e) * (1 - w4 * e * e * De)
        v = me_GeV * math.exp(C2 - N * (math.pi / 512) / (2 * math.pi) + w16 * e * e) / Fe
        return (v / v_meas - 1) / sv

    # G depends on every menu except w16: M_Pl(Newton) = sqrt(c) M_b
    expo = {(kb, ka, kr): _exponent(b0, al, d)
            for (kb, b0), (ka, al), (kr, d) in itertools.product(M["b0"].items(), M["alpha_G2"].items(), M["reading"].items())}

    def count(w4):
        """The common mode count P v^2 P, evaluated on the lane's depth four lepton rates."""
        rates = [g * g * (1 - w4 * e * e * g) for g in gap.values()]
        q0 = (sum(math.sqrt(x) for x in rates) / 3) ** 2
        return (sum(rates) / 3 - q0) / q0
    assert abs(count(r.d11**3 / r.d10) / r.COMMON_MODE - 1) < 1e-10

    def g_pull(ce, S, c, c2, w4):
        Mb = me_GeV / (ce * math.exp(-S) * De**2 * (1 + c2 * e) * (1 - w4 * e * e * De))
        return (count(w4) * (CODATA_MPL / (math.sqrt(c) * Mb))**2 - 1) / SIGMA_G     # exact

    n_lanes = 0
    close_G, close_GF, close_mu, close_all = 0, 0, 0, []
    gf_mu_lanes_G = []                                   # G pulls of lanes closing on G_F and m_mu
    gf_cache = {}
    for (kc, ce), (kc2, c2), (kw4, w4) in itertools.product(M["C_EFF"].items(), M["c2"].items(), M["w4"].items()):
        gf_ok_w16 = [k16 for k16, w16 in M["w16"].items()
                     if abs(gf_cache.setdefault((kc, kc2, kw4, k16), gf_pull(ce, c2, w4, w16))) < 2.0]
        for (kS, S), (kcc, c) in itertools.product(expo.items(), M["c"].items()):
            gp = g_pull(ce, S, c, c2, w4)
            g_ok = abs(gp) < 3.0
            n16 = len(M["w16"])
            n_lanes += n16
            close_G += n16 * g_ok
            close_GF += len(gf_ok_w16)
            close_mu += n16 * mu_ok[kw4]
            if gf_ok_w16 and mu_ok[kw4]:
                gf_mu_lanes_G.append((abs(gp), (kc, *kS, kcc, kc2, kw4, tuple(gf_ok_w16))))
                if g_ok:
                    close_all.append((kc, *kS, kcc, kc2, kw4, tuple(gf_ok_w16), round(gp, 2)))
    gf_mu_lanes_G.sort()
    canon = ("1/2", "32/3", "P h/2pi", "action", "1", "h_bridge", "27/2", ("16",), round(g_pull(0.5, r.S_lepton, 1.0, -1.0, 13.5), 2))
    assert close_all == [canon], close_all
    runner = next(g for g, lane in gf_mu_lanes_G if lane[:7] != canon[:7])
    assert runner > 900
    assert close_G > 100                                   # E1: G alone admits many lanes
    pinned = {lane[0] for _, lane in gf_mu_lanes_G}, {lane[5] for _, lane in gf_mu_lanes_G}, \
        {lane[6] for _, lane in gf_mu_lanes_G}, {lane[7] for _, lane in gf_mu_lanes_G}
    assert pinned == ({"1/2"}, {"h_bridge"}, {"27/2"}, {("16",)})

    report(f"G-path lanes crossed: {n_lanes} "
           f"({' x '.join(str(len(v)) for v in M.values())})")
    report(f"E1 close on G alone (3 sigma): {close_G}. On G_F alone (2 sigma): {close_GF}; "
           f"on m_mu alone (2 sigma): {close_mu}")
    report(f"E2 G_F and m_mu together pin C_EFF {sorted(pinned[0])}, c2 {sorted(pinned[1])}, "
           f"w4 {sorted(pinned[2])}, w16 {sorted(pinned[3])}, with no G")
    report(f"E3 close on all three: {len(close_all)}: {close_all}. The nearest other lane closing on "
           f"G_F and m_mu misses G by {runner:.3g} sigma ({gf_mu_lanes_G[1][1][1:5]})")
    return dict(n_lanes=n_lanes, close_G=close_G, close_GF=close_GF, close_mu=close_mu,
                close_all=close_all, runner_sigma=runner, pinned=pinned)


if __name__ == "__main__":
    run()
