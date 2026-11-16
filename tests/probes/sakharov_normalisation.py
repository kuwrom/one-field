"""The Sakharov normalisation, made exact.

The canonical gravity ledger prints the diagnostic
``G_ind/G_N = (2 pi d11)/Sigma`` under the cutoff t0 = d11 l_Pl^2
(FACT 7). With G anchored the regulator is the boundary identification
of root.py.  This probe computes
``Sigma_target = 2 pi d11 = 6 pi``, the cutoff scheme ratios, two
factorizations of the target, the q=0 dispersion and the vacuum rebase.
Seven facts.

FACT 1 (the algebra anchor).  The standard proper time heat kernel
gives 1/G_ind = Sigma/(2 pi t0) with t0 the proper time cutoff and
Sigma the a1-weighted dof sum.  One minimally coupled scalar with a
momentum cutoff (t0 = 1/Lambda^2) reproduces the textbook Sakharov
value 1/G = Lambda^2/(12 pi).  The imported mathematics is standard
(DeWitt Seeley). Nothing engine specific enters here.

FACT 2 (the exact target identity).  G_ind = G_N at Sigma = 6 pi is,
through FACT 1, exactly

    Sigma_target = 2 pi d11,   equivalently   t0 = d11 l_Pl^2.

The target is an expression in d11.

FACT 3 (the cutoff schemes).  The substrate natural cutoffs give
G_ind/G_N = (target implied by t0) / (counted Sigma), in the exact
ratios 1 : 2 : 4 : 12 for t0 = xi0^2, 2 xi0^2, l_Pl^2, d11 l_Pl^2
(with xi0 = l_Pl/2).  With the standard per state weights the counted
Sigma is 27, so the d11 scheme reads 2 pi/9 = 0.698 and the unfactored
schemes read 0.058, 0.116, 0.233.

FACT 4 (one cutoff identity written in two equivalent forms).
  (i)  channel coupling reading, the law behind alpha_G2
       (coupling = closure probability x local response):
           1/G = P * Sigma/(2 pi xi0^2)
       with P = 1/(d10^2 d11) the wiring derived pointed closure
       probability and xi0 = C_EFF l_Pl = l_Pl/d10 the trace lifted
       Planck length.  The d10^2 cancels EXACTLY between P and
       xi0^2, leaving Sigma_target = 2 pi d11 for any d10.
  (ii) common mode share reading: the metric is the q = 0 one of the
       d11 Z3 channels of the gradient bilinear
       (tests/probes/composite_metric.py),
       1/G = (1/d11) Sigma/(2 pi l_Pl^2).
These are the same identity: P = (xi0/l_Pl)^2/d11, so
xi0^2/P = d11 l_Pl^2. Both readings therefore give the same proper time
cutoff t0 = d11 l_Pl^2 and the same target 2 pi d11. No discriminator is
needed.

FACT 5 (rebase and dispersion of a uniform medium, gas sign g0 > 0).
The length xi_0 = C_EFF l_Pl = l_Pl/d10 (root.py) is the length under which
FACT 4's d10 cancellation fires.  The healing length of the
Z3-symmetric uniform state and that of the clock condensate,
mu_clock = (3 g0 - g1) rho0, differ by the rebase factor
sqrt(mu_sym/mu_clock) = 1.3875, computed from the substrate's BdG
constants.  The q=0 phase/phonon dispersion is verified here against
the substrate's BdG: omega^2 = mu p^2 + p^4/4 (the k^4 coefficient
is exactly 1/4 at hbar = m = 1, crossover 2/xi, so xi^2 = 1/mu here;
record_area.py reads xi_0 = l_Pl/d10 in the range convention
xi^2 = 1/(4 mu), half this length). This gapless response
is the sound mode of a uniform medium.  The geometry channel is the
gradient bilinear of the carrier in the empty field
(tests/probes/composite_metric.py).

FACT 6 (no coefficient is read from the ledger; three derivation routes
are closed).  With the standard per state weights the ledger as entered gives
6 pi/Sigma = 0.698 under a covariant regulator and t0 = d11 l_Pl^2 (a diagnostic, FACT 7).  Closing at 6 pi would need
either a per channel weight 0.122 in place of the minimal 1/6, or
t0 = 4.297 l_Pl^2.  None of the places such a number could come from
supplies it:
  (a) Depth.  Every post searchlight term in the web reads
      eq_loop = (alpha/2pi)^2. The depth-4 terms are 5e-7 to 3e-6
      relative and a depth-5 echo would move Sigma by ~4e-8.  The
      gap is 0.43.  No stiffness node reads eq_loop (the bridge carries
      the loop and receives no echo of it), asserted below. The common
      mode's own node reads it once, at 3.6e-5 (root.py).
  (b) The Albert junction.  In the Tits construction the bridge is
      Im(O) x J3(O)_0, traceless in both factors, so a channel has
      zero linear projection on the unit and its quadratic vent into
      the unique (1,1) invariant is its own trace norm.  The depth-4
      transport is an Albert automorphism, block diagonal on 1 + 26,
      so the per channel multiplier is exactly 1 (theorem; the
      explicit J3(O) check is kept with this probe).
  (c) The substrate.  The q=0 channel is shift symmetric, so xi = 0
      exactly. The quantum pressure corrections (1/4) grad^2 ln rho
      and p^4/4 multiply p^2 and are metric corrections, not curvature
      couplings.  The p^4 dispersion softens the induced-1/G
      divergence from 1/t0 to sqrt(pi) Lambda_h/sqrt(t0) and does not
      remove it (asserted below). It belongs to a condensate, not to the
      vacuum. A condensate filling space needs xi > 2.7e6 l_Pl to
      keep c_T = c (carrier_sakharov.py). The record held inside a lump
      regulates only for a Planckian carrier
      (closed_event_scale.py).  t0 stays an identification (Volovik:
      the numerical factor depends on the cut off procedure).
  (d) The spectrum.  Only the bridge count 182 and the 14 G2 vectors
      are forced. The ledger as entered, with the G2 scalars light,
      gives 0.698. The engine's one list, with the G2 scalars as the R6
      composite, gives 0.730.  The amount has no G2 vertex, so no radiative xi
      is available to it (interference/boundary.py, THE AMOUNT'S
      COUPLING).
Result: the bridge's sign repair and the lane are derived. The ledger
serves as a closure lane, and no coefficient is read from it. In
gravity.py the running and the stiffness are both read on the one list,
with the bridge as amount, and the six gates pass
(tests/probes/horizon_amount.py). The ledger as entered passes them too.

FACT 7 (what the ratio is, and what the framework's G prediction is).
G cancels out of G_ind/G_N: the ratio is 2 pi (t0/l_Pl^2)/Sigma, a
dimensionless consistency check between two identifications, the
anchor inversion m_e -> M_Pl (the Planck tie, root.py) and the
proper time cutoff t0 = d11 l_Pl^2.  The framework's prediction of
G is the anchor inversion itself: M_Pl at +22 ppm, G at -55 ppm, with
no gravitational dynamics in it.  The Sakharov ratio adds the
O(1) statement that the induced stiffness is Planckian and attractive.
Its number is regulator dependent.  Asserted below: under the same
rows, the vector contact coefficient +1 (covariant proper time: Kabat's
contact term, which Donnelly-Wall 1412.1895 identify as the entanglement
entropy of the gauge field's edge modes), 0 (hard momentum cutoff on
the physical polarizations) and -1 give 0.698,
0.419 and 0.299, all attractive; a smooth regulator shape
f(u) = (1 - e^-u)^2 rescales 1/G by int f/u^2 = 2 ln 2 and gives 0.504.
What survives every scheme is the sign and the order of magnitude.
Every row couples to one universal tensor metric: the gradient
bilinear of the carrier, written without a background. Each row's
weight a1 is its share of the stiffness of that metric's tensor mode
(tests/probes/composite_metric.py).

CLOSURE STATUS.  The heat kernel arithmetic and coefficient are
computed here, with the q=0 dispersion and rebase of a uniform medium
on the gas sign (FACT 5).  The proper time cutoff t0 = d11 l_Pl^2 is an
identification. This diagnostic compares it with the anchored G. With G
anchored the regulator is the boundary identification of root.py,
t0 = Sigma/(2 pi) l_Pl^2 (g_path_trials.py G10).  Given a cutoff, the
DeWitt Seeley theorem fixes the arithmetic on the composite metric.
The coefficient it returns is FACT 6's.

Usage: python3 sakharov_normalisation.py
"""

import math
import os
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
for _path in (os.path.dirname(_here),
              os.path.join(os.path.dirname(os.path.dirname(_here)),
                           "interference"),
              os.path.join(_here, "interference")):
    if _path not in sys.path:
        sys.path.insert(0, _path)
import _bootstrap  # noqa: E402,F401

from root import (d10, d11, n7, n26,   # noqa: E402
                  xi_bridge, C_EFF)
from nls_soliton import (bogoliubov_omega_sq,        # noqa: E402
                         _bdg_analysis, G1_OVER_G0)

G_RATIO_FROZEN = 0.698131701


def engine_sigma():
    """The venting ledger, rebuilt from the four integers with the
    standard per state weights (scalar +1/6, Weyl +1/12, massless
    vector −1/3, massive vector −1/6: Visser gr-qc/0204062 Table 1)."""
    a_s, a_f, a_g = 1.0 / 6, 1.0 / 12, -1.0 / 3
    n_gen = d11
    rows = [
        ("gauge", 2 * d11 * (d11 + 1), a_g),
        ("fermions", d10 * n_gen * (d11 * d10**2 + d10 + 1), a_f),
        ("coset", d10 * (2 * n7 - (d11**2 - 1)), a_g),
        ("G2-scal", n7, a_s),
        ("bridge", n7 * n26, a_s),            # the amount, minimally coupled (boundary.xi_of)
    ]
    sigma = sum(n * a for _, n, a in rows)
    return sigma, dict((k, n * a) for k, n, a in rows)


def run(report=print):
    ## FACT 1: the anchor #
    # 1/(16 pi G) = a1 Sigma/(32 pi^2 t0)  ->  1/G = Sigma/(2 pi t0);
    # one scalar, a1 = 1/6, t0 = 1/Lambda^2  ->  1/G = Lambda^2/(12 pi)
    one_scalar = (1.0 / 6) * 16.0 * math.pi / (32.0 * math.pi**2)
    assert abs(one_scalar - 1.0 / (12.0 * math.pi)) < 1e-15
    report("FACT 1: 1/G = Sigma/(2 pi t0). One scalar reproduces the")
    report("  textbook Sakharov 1/G = Lambda^2/(12 pi).")

    ## FACT 2: the exact identity ##
    sigma, rows = engine_sigma()
    assert abs(6.0 * math.pi - 2.0 * math.pi * d11) < 1e-18
    g_ratio = 6.0 * math.pi / sigma
    assert abs(g_ratio - G_RATIO_FROZEN) < 1e-8
    report(f"FACT 2: Sigma = {sigma:.6f}, target = 2 pi d11 = 6 pi, "
           f"G_ind/G_N = {g_ratio:.9f}")

    ## FACT 3: the naive schemes ##
    # Sigma_target(t0) = 2 pi t0 M_Pl^2. Ratios quoted against the
    # engine Sigma, with xi0 = l_Pl/2 (the prose identification)
    schemes = {"t0 = xi0^2": 0.25, "t0 = 2 xi0^2": 0.5,
               "t0 = l_Pl^2": 1.0, "t0 = d11 l_Pl^2": 3.0}
    scheme_ratios = {k: 2.0 * math.pi * v / sigma
                     for k, v in schemes.items()}
    assert abs(scheme_ratios["t0 = xi0^2"] - g_ratio / 12) < 1e-9
    assert abs(scheme_ratios["t0 = l_Pl^2"] - g_ratio / 3) < 1e-9
    assert abs(scheme_ratios["t0 = d11 l_Pl^2"] - g_ratio) < 1e-9
    report("FACT 3: naive schemes: "
           + ", ".join(f"{k} -> {v:.4f}" for k, v in
                       scheme_ratios.items()))

    # # FACT 4: the two equivalent factorizations #
    # (i) P = 1/(d10^2 d11), xi0 = C_EFF l_Pl = l_Pl/d10:
    #     Sigma_target = 2 pi (xi0/l)^2 / P = 2 pi d11, d10 cancels
    assert C_EFF == 1.0 / d10
    for d10_test in (1, 2, 3, 5):
        p_test = 1.0 / (d10_test**2 * d11)
        xi_over_l = 1.0 / d10_test
        assert abs(p_test - xi_over_l**2 / d11) < 1e-15
        tgt = 2.0 * math.pi * xi_over_l**2 / p_test
        assert abs(tgt - 2.0 * math.pi * d11) < 1e-12
    # (ii) 1/G = (1/d11) Sigma/(2 pi l^2): Sigma_target = 2 pi d11
    assert abs(2.0 * math.pi * 1.0 * d11 - 6.0 * math.pi) < 1e-18
    report("FACT 4: factorizations (i) P-closure x the area (l_Pl/d10)^2 "
           "(d10 cancels) and (ii) 1/d11 share x Planck area: both "
           "are t0=d11 l_Pl^2 and give 2 pi d11 exactly.")

    # Leading Einstein-Hilbert coefficient:
    # Gamma^(2) = integral sqrt(g) [Sigma/(32 pi^2 t0)] R.
    # In Planck units t0/l_Pl^2=d11. Exact closure gives 1/(16 pi).
    eh_coefficient_lpl2 = sigma / (32.0 * math.pi**2 * d11)
    eh_target_lpl2 = 1.0 / (16.0 * math.pi)
    assert abs(eh_coefficient_lpl2 / eh_target_lpl2
               - 1.0 / g_ratio) < 1e-12
    report(f"  EH coefficient*l_Pl^2 = {eh_coefficient_lpl2:.10f}; "
           f"target 1/(16 pi) = {eh_target_lpl2:.10f}.")

    # # FACT 5: substrate premise audit #
    g0, rho0 = 1.0, 1.0
    g1 = G1_OVER_G0 * g0
    bd = _bdg_analysis(g0, g1, rho0)
    mu_sym = bd["a_q"][0]
    assert abs(mu_sym - (3 * g0 + 2 * g1) * rho0) < 1e-12
    # q=0 dispersion: omega^2 - mu p^2 = p^4/4, exactly, all p
    for p in (0.3, 1.0, 2.7, 8.0):
        w2 = float(bogoliubov_omega_sq(p, g0, g1, rho0, q=0))
        assert abs((w2 - mu_sym * p**2) - p**4 / 4.0) < 1e-9
    xi_sym = 1.0 / math.sqrt(mu_sym)
    mu_clock = (3 * g0 - g1) * rho0
    rebase = math.sqrt(mu_sym / mu_clock)
    assert abs(rebase - 1.3875) < 1e-3
    report(f"FACT 5: q=0 dispersion omega^2 = mu p^2 + p^4/4 exact "
           f"(mu = {mu_sym:.5f}, xi = {xi_sym:.5f}). Clock vacuum "
           f"rebase sqrt(mu_sym/mu_clock) = {rebase:.4f} (uniform medium, "
           f"gas sign). Given the cutoff identification, the standard "
           f"theorems fix the Einstein-Hilbert arithmetic on the composite "
           f"metric (FACT 6).")

    # # FACT 6: the closed derivation routes, executable parts #
    # (a) no stiffness node reads the depth-3 carrier. The depth-4 web
    #     terms are ppm scale against a 0.43 gap
    from root import WEB, Ledger, SEARCHLIGHT_DEPTH
    grav_reads_carrier = any(
        t.parent == "eq_loop"
        for n in ("a1_bridge", "Sigma_vent(UV)", "G_ratio")
        if isinstance(WEB.get(n), Ledger) for t in WEB[n].terms)
    assert not grav_reads_carrier
    st = WEB.state
    post = [abs(t.value(st) / (led.base(st) if callable(led.base)
                                else led.base))
            for led in WEB.values() if isinstance(led, Ledger)
            for t in led.terms
            if t.depth > SEARCHLIGHT_DEPTH and led.mode == "add"
            and (led.base(st) if callable(led.base) else led.base)]
    max_depth4_rel = max(post) if post else 0.0
    gap = abs(1.0 / g_ratio - 1.0)
    assert max_depth4_rel < 1e-4 < gap
    # (c) the p^4 dispersion does not regularize the R coefficient:
    #     C(t0) = int dlam E1(t0 (lam + lam^2/L^2)) ~ sqrt(pi) L/sqrt(t0)
    from scipy.special import exp1
    from scipy.integrate import quad
    #     with subleading term (L^2/2)(ln(t0 L^2) + gamma - 2)
    L = 2.0 / xi_sym                       # Lambda_h = 2/xi0
    consts = []
    for t0 in (1e-4, 1e-6):
        C = quad(lambda lam: exp1(t0 * (lam + lam**2 / L**2)),
                 0.0, np.inf, limit=800)[0]
        pred = (math.sqrt(math.pi) * L / math.sqrt(t0)
                + 0.5 * L**2 * (math.log(t0 * L**2) + np.euler_gamma - 2))
        assert abs(C / pred - 1.0) < 1e-3
        consts.append(C * math.sqrt(t0) / L)
    assert consts[-1] > consts[0]          # approaches sqrt(pi) from below
    report(f"FACT 6: no stiffness node reads eq_loop. Max depth-4 term "
           f"{max_depth4_rel:.1e} relative vs gap {gap:.2f}; p^4 "
           f"dispersion leaves C(t0) sqrt(t0)/Lambda_h -> "
           f"{consts[-1]:.3f} (sqrt(pi) = {math.sqrt(math.pi):.3f}), "
           f"divergent.  The coefficient is an identification.")

    # # FACT 7: regulator dependence of the coefficient, executable #
    # per massless vector state the weight is 1/6 - c/2 with c the
    # contact coefficient (c = 1 covariant, Kabat, the edge mode entropy
    # of Donnelly-Wall 1412.1895; c = 0 hard momentum cutoff on the
    # physical polarizations). The bridge row is unchanged
    n0 = n7
    n_weyl = d11 * (d11 * d10**2 + d10 + 1)
    n_vec = d11 * (d11 + 1) + (2 * n7 - (d11**2 - 1))
    bridge_row = n7 * n26 * (1.0 / 6)
    regulator = {}
    for c in (1.0, 0.0, -1.0):
        sig_c = (n0 / 6.0 + n_weyl / 6.0 + 2 * n_vec * (1.0 / 6 - c / 2.0)
                 + bridge_row)
        regulator[c] = 6.0 * math.pi / sig_c
    assert abs(regulator[1.0] - g_ratio) < 1e-12
    assert all(v > 0 for v in regulator.values())      # attractive throughout
    shape_K = quad(lambda u: (1.0 - math.exp(-u))**2 / u**2, 0.0, np.inf)[0]
    assert abs(shape_K - 2.0 * math.log(2.0)) < 1e-6
    report(f"FACT 7: contact coefficient +1/0/-1 -> G_ind/G_N = "
           f"{regulator[1.0]:.3f}/{regulator[0.0]:.3f}/{regulator[-1.0]:.3f}, "
           f"all attractive. Smooth regulator (1-e^-u)^2 rescales by "
           f"2 ln 2 -> {g_ratio / shape_K:.3f}.  The framework's G "
           f"prediction is the anchor inversion (M_Pl +22 ppm). The "
           f"ratio is a consistency check between identifications.")

    return {
        "sigma": sigma,
        "sigma_rows": rows,
        "g_ratio": g_ratio,
        "gravity_reads_carrier": grav_reads_carrier,
        "max_depth4_relative": max_depth4_rel,
        "p4_asymptotic_constant": consts[-1],
        "regulator_variants": regulator,
        "regulator_shape_factor": shape_K,
        "target_is_2pi_d11": True,
        "scheme_ratios": scheme_ratios,
        "d10_cancellation_exact": True,
        "factorizations_equivalent": True,
        "eh_coefficient_lpl2": eh_coefficient_lpl2,
        "eh_target_lpl2": eh_target_lpl2,
        "mu_symmetric": mu_sym,
        "mu_clock": mu_clock,
        "k4_coefficient": 0.25,
        "xi_symmetric": xi_sym,
        "clock_rebase_factor": rebase,
    }


if __name__ == "__main__":
    run()
