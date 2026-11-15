"""The three coefficients of the ruler: each one an object the web already holds.

THE QUESTION. The anchor inversion reads the electron as
    m_e = M_b C_EFF e^(-S0) Delta_e^2 (1 + c2 alpha/2pi) (1 - w4 L Delta_e),
with L = (alpha/2pi)^2 the closed carrier. Three coefficients stand in
it, C_EFF = 1/2, c2 = -1 and w4 = 27/2. The audits list them as
held by measurement (tests/audit/g_path_trials.py: C_EFF by G_F, c2
by G_F, w4 by m_mu). Is each one an object the web already holds, read
the same way the web reads every other object? And does each
alternative fail by a computation?

PREMISES. P1 the forgetting: P v P = 0 and P v^2 P = 1, the common mode
entering at second order (interference/searchlight/protected_forgetting.py).
P2 the boundary reading: the common mode reads each sum whole
(tests/probes/boundary_lattice.py). P3 the single response definition,
R(h) = h/(2 pi), with h_bridge = h_7 + h_26 = 1 (interference/wiring.py).
P4 the vertex composition rule: "echoes compose existing vertices, they
never invent new ones" (root.py), and the recursive Albert transport
fixes its unit (tests/probes/wiring_theorem_probe.py PART D).
P5 G_F and m_mu/m_e as measured, the only data read, as tests.

STATEMENTS (the numbers below check each one).

  F1 (C_EFF is the closure's common mode read whole). The rates are
     Delta_k^2 = (1 + v_k)^2, and their common mode at second order is
     P Delta^2 P = 1 + 2 P v P + P v^2 P = 2. The frame reads that
     sum whole: C_EFF P Delta^2 P = 1, so C_EFF = 1/(1 + P v^2 P) = 1/2.
     This is the channel weight W(fund) = 1/d10 that root.py names.
     It is C_EFF = 1/(3Q) with Q = 2/3. FAILS for each alternative:
     the first order reading (P Delta P)^2 = 1 gives C_EFF = 1 and
     v_EW/m_e off by a factor two. On the menu {4/9, 1/4, 2/9} the
     Koide value would be 3/4, 4/3, 3/2, and Q above one admits no
     positive circulant. Q = 3/4 is excluded by the lepton ratios,
     not by a G-free gate. C_EFF, E[v^2], w4 and the boundary
     reading's base are one identification, B/A: c = 2/B^2 and
     C_EFF = 1/(1 + B^2/2).
  F2 (c2 is the bridge vertex 1/alpha already reads). The depth one
     edge of 1/alpha is -R(h_bridge) on its base. The lepton's
     depth two echo is the same vertex read through the echoed
     coupling, -alpha R(h_bridge) = -alpha/(2 pi), so
     c2 = -2 pi R(h_bridge) = -1. FAILS for the other conformal weights
     of the web, h_7 = 2/5, h_26 = 3/5, h11 = 1/2 and h10 = 2/9. The
     response is defined only for a dimension one current, and the rule
     enforces that definition (wiring.normalized_current_response). Forced in by hand, G_F
     excludes each by more than 900 standard deviations.
  F3 (w4 is the one channel share of the Albert traversal). The
     lepton's own share is C_EFF of F1. The recursive Albert transport
     includes its unit, dim J3(O) = n26 + 1 = 27, so w4 = C_EFF x 27 =
     27/2. FAILS for the traversal without its unit, 26 (w4 = 13), at
     16 standard deviations of m_mu/m_e, and for the whole traversal,
     27, at 440.
  F4 (the tests). With the three coefficients as read, v_EW/m_e matches
     G_F at +0.08 standard deviations and m_mu/m_e comes out at +0.12.

RESULT. The three coefficients of the ruler are the closure's common
mode read whole, the bridge vertex that 1/alpha reads, and the
one channel share of the Albert traversal. G_F and m_mu test them, and
each alternative fails. c2's alternatives fail by the response's
definition and by G_F. C_EFF's fail by G_F (the first order reading)
and the Koide bound. w4's fail by the lepton ratio. Each coefficient
remains an identification. The frame's reading of the closure's common
mode, C_EFF P Delta^2 P = 1, is a normalization that the data test,
not one that a gate computes. C_EFF and w4 are B/A = sqrt(d10) read
again through P v^2 P = 1, not a second identification. The channel
map of B/A remains an identification (interference/__init__.py, the
selection table).

Scope: the electron's ruler at the depth the engine reads it. The
coefficients are read on the canonical integers, where 1/(1 + d10/2)
and 1/d10 coincide.

Usage: python3 tests/probes/frame_coefficients.py
"""

import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import wiring_theorem_probe as wtp  # noqa: E402
from interference import root as r  # noqa: E402
from interference import wiring  # noqa: E402
from interference.searchlight import protected_forgetting as pf  # noqa: E402

GF = (1.1663788e-5, 0.0000006e-5)            # GeV^-2, PDG 2024
MU_E = (206.7682827, 0.0000046)              # CODATA 2022


def ruler(c_eff=None, c2=-1.0, w4=None):
    """v_EW/m_e and m_mu/m_e for the given coefficients, every other entry as the engine reads it."""
    c_eff = r.C_EFF if c_eff is None else c_eff
    w4 = r.d11**3 / r.d10 if w4 is None else w4
    e, L = r.alpha_phys / (2 * math.pi), r.WEB.state["eq_loop"]
    gap = lambda k: 1 + math.sqrt(r.d10) * math.cos(float(r.h10) + 2 * math.pi * k / 3)
    rate = lambda k: c_eff * math.exp(-r.S_lepton) * gap(k)**2 * (1 + c2 * e) * (1 - w4 * L * gap(k))
    v_over_me = math.exp(-r.WEB.state["S_quark"]) / rate(1)
    return v_over_me, rate(2) / rate(1)


def pulls(c_eff=None, c2=-1.0, w4=None):
    """(G_F pull, m_mu/m_e pull) of the ruler; v = (sqrt2 G_F)^(-1/2), m_e as CODATA."""
    v_over_me, mu_e = ruler(c_eff, c2, w4)
    v_meas = (math.sqrt(2) * GF[0]) ** -0.5 * 1e3 / r.M_E_ANCHOR_MEV
    gf = math.log(v_over_me / v_meas) / (GF[1] / (2 * GF[0]))
    return gf, (mu_e - MU_E[0]) / MU_E[1]


def run(report=print):
    out = {}
    v = pf.harmonic()
    PvP, Pv2P, _ = pf.moments(v)
    delta = 1 + v
    P_delta2 = float((delta**2).mean())
    assert abs(PvP) < 1e-12 and abs(Pv2P - 1) < 1e-12 and abs(P_delta2 - (1 + 2 * PvP + Pv2P)) < 1e-12

    # F1: C_EFF.
    c_eff = 1 / P_delta2
    Q = float(r.Q0)
    assert abs(c_eff - r.C_EFF) < 1e-12 and abs(c_eff - 1 / r.d10) < 1e-12 and abs(c_eff - 1 / (3 * Q)) < 1e-12
    first_order = 1 / float(delta.mean())**2
    menu = {"4/9": 4 / 9, "1/4": 1 / 4, "2/9": 2 / 9}
    menu_Q = {k: 1 / (3 * c) for k, c in menu.items()}
    menu_B2 = {k: 6 * q - 2 for k, q in menu_Q.items()}
    assert abs(first_order - 1) < 1e-12 and abs(ruler(first_order)[0] / ruler()[0] - 0.5) < 1e-12
    assert menu_Q["1/4"] > 1 and menu_Q["2/9"] > 1 and abs(menu_B2["4/9"] - 2.5) < 1e-12
    gf_first = pulls(c_eff=first_order)[0]
    out.update(c_eff=c_eff, P_delta2=P_delta2, first_order=first_order, menu_Q=menu_Q, menu_B2=menu_B2,
               gf_first=gf_first)
    report(f"F1 P Delta^2 P = {P_delta2:.12f}, C_EFF = {c_eff:.12f} = 1/d10 = 1/(3Q). First order reading "
           f"C_EFF = {first_order:.3f} ({gf_first:+.2e} sigma of G_F). Menu Q: "
           + ", ".join(f"{k} -> {q:.3f}" for k, q in menu_Q.items()) + f". Q = 3/4 needs B^2 = {menu_B2['4/9']:.2f}")

    # F2: c2.
    R = wiring.normalized_current_response(float(r.h_bridge))
    edge = r.WEB["inv_alpha"].terms[0].value(r.WEB.state) / r._inv0
    echo = r.WEB["lepton_EM"].terms[0].value(r.WEB.state)
    c2 = -2 * math.pi * R
    assert abs(edge + R) < 1e-15 and abs(echo + r.alpha_phys * R) < 1e-15 and abs(c2 + 1) < 1e-15
    weights = {"h_bridge": float(r.h_bridge), "h_7": float(r.h_7), "h_26": float(r.h_26),
               "h11": float(r.h11), "h10": float(r.h10)}
    admitted = []
    for k, h in weights.items():
        try:
            wiring.normalized_current_response(h)
            admitted.append(k)
        except ValueError:
            pass
    c2_pulls = {k: pulls(c2=-h)[0] for k, h in weights.items()}    # forced in by hand, R = h/(2 pi)
    assert admitted == ["h_bridge"]
    assert abs(c2_pulls["h_bridge"]) < 1 and min(abs(p) for k, p in c2_pulls.items() if k != "h_bridge") > 900
    out.update(R=R, c2=c2, admitted=admitted, c2_pulls=c2_pulls)
    report(f"F2 1/alpha's depth one edge -{R:.6f} and the lepton echo -alpha x {R:.6f}: c2 = {c2:.0f}. The response "
           f"rule admits {', '.join(admitted)} alone. Forced in by hand, G_F pulls: "
           + ", ".join(f"{k} {p:+.1f}" for k, p in c2_pulls.items()))

    # F3: w4.
    albert = wtp.albert_unital_response(n26=r.n26, schur_weight=Fraction(c_eff).limit_denominator())
    w4 = float(albert["response"])
    w4_menu = {"27/2 (share of the unital traversal)": w4, "13 (no unit)": c_eff * r.n26,
               "27 (the whole traversal)": float(r.n26 + 1)}
    w4_pulls = {k: pulls(w4=w)[1] for k, w in w4_menu.items()}
    assert w4 == r.d11**3 / r.d10 and albert["junction"] == 1
    assert abs(w4_pulls["27/2 (share of the unital traversal)"]) < 1
    assert abs(w4_pulls["13 (no unit)"]) > 10 and abs(w4_pulls["27 (the whole traversal)"]) > 100
    out.update(w4=w4, w4_pulls=w4_pulls)
    report("F3 w4 = C_EFF x dim J3(O) = " + f"{w4}" + ". m_mu/m_e pulls: "
           + ", ".join(f"{k} {p:+.1f}" for k, p in w4_pulls.items()))

    # F4: the tests.
    gf, mu = pulls()
    assert abs(gf) < 1 and abs(mu) < 1
    out.update(gf=gf, mu=mu)
    report(f"F4 as read: G_F {gf:+.2f} sigma, m_mu/m_e {mu:+.2f} sigma")
    report("RESULT: C_EFF, c2 and w4 are the closure's common mode read whole, the bridge vertex 1/alpha reads, "
           "and the one channel share of the Albert traversal. G_F and m_mu test them.")
    return out


if __name__ == "__main__":
    run()
