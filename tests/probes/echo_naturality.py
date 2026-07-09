"""
echo_naturality.py: the response rule. One depth-3 carrier, three
depth-4 channel readings, and the kernel regenerates canon without
three independent loop declarations.

SCOPE. The depth-4 responses (2, 16, 27/2) are not three
separately declared edges: the wiring generates them from one loop
object, and this probe exhibits the generation.

WHAT THIS PROBE IS.  A refactoring-equivalence check, not an
independent derivation.  The response weights below are algebraically
equal to canon's depth-4 coefficients, so regeneration is guaranteed
once the factorisation is admitted.  What is established here:
(i) UNIFORMITY: all three depth-4 readings are expressible as ONE loop
amplitude L = (alpha/2pi)^2 read through each node's own channel
response, with no per-node declaration left anywhere; (ii)
EQUIVALENCE: the refactored web's fixed point coincides with the
LIVE canon web (root.py WEB) and with the frozen
precision block. The wiring theorem proves the 2-cycle is the only
depth-3 carrier and derives the unital Albert lift; vertex_pinning.py retains the
independent numerical consequences.

THE RESPONSE RULE (the composition law, stated uniformly):

    node reads L x w_e(n) x w_q(n),  where
      w_e(n) = the node's own e-side channel weight,
      w_q(n) = the trace of the q-sector response operator in the
               node's own channel,
    and the coupling node, whose channel is the loop's own carrier,
    with no matter vertex, reads the unweighted basepointed-orientation
    count of the 2-cycle instead (= 2, the wiring theorem's counting
    function).

    Every weight is wiring-layer data or an already-FORCED vertex:
      coupling  : 2      = cycle orientation count      [wiring thm]
      action    : 8/3    = charge trace                 [the depth-1
                  x 6    = C2 of the branching rep 26    alpha_G2 ->
                                                         alpha_EM
                                                         conversion
                                                         vertex; the
                                                         embedding's
                                                         channel
                                                         Casimir]
      amplitude : 1/2    = W(fund) Schur projector      [the B/A
                  x 27   = dim J3(O) = 1 + n26 = d11^3   theorem,
                                                         projector
                                                         lemma; the
                                                         q-face
                                                         Jordan
                                                         algebra of
                                                         the
                                                         condensation]

THE COMPUTATION (this probe).  A fresh Web is built from canon's own
kernel classes with the three independent loop declarations DELETED and one
new node added: eq_loop, base 0, whose single edge is the loop
amplitude read self-consistently from the state.  Every former
former response becomes a depth-4 READ of eq_loop weighted by the
response rule.  The kernel is solved jointly, and the fixed point is
required to equal

    (a) the LIVE canon web (root.py WEB): inv_alpha, S_quark, M_Pl;
    (b) the frozen precision block:
        1/alpha(0) = 137.035999050
        m_mu       = 105.6583758 MeV
        m_tau      = 1776.9092814 MeV
        v_EW       = 246.219645 GeV

PASSING MEANS: the depth-3 carrier is one loop object under one depth-4
read rule; no per-node carrier declaration is needed anywhere.

STANDING. This file proves refactoring equivalence and uniformity. The
unique-loop/orientation theorem and unital operator closure are in
wiring_theorem_probe.py;
vertex_pinning.py records the data contact.

Usage: python3 echo_naturality.py
"""

import math
import os
import sys
from fractions import Fraction

_here = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.dirname(_here),
          os.path.join(os.path.dirname(os.path.dirname(_here)),
                       "interference"),
          os.path.join(_here, "interference")):
    if p not in sys.path:
        sys.path.insert(0, p)
import _bootstrap  # noqa: E402,F401

# the canonical precision block (frozen values)
TARGETS = dict(inv_alpha=137.035999050, m_mu=105.6583758,
               m_tau=1776.9092814, v_EW=246.219645)


try:
    from probes.wiring_theorem_probe import orientation_multiplicity
except ImportError:
    from wiring_theorem_probe import orientation_multiplicity


def run(report=print):
    from root import (Ledger, Web, d10, d11, n26, h10,
                      h_bridge, charge_trace, C2_26, S_lepton,
                      vertex_charge_trace, alpha_EM, BA_ratio,
                      M_E_ANCHOR_MEV)

    report("ECHO NATURALITY: one loop object, three responses, canon "
           "regenerated")
    report("=" * 70)

    # # the response rule's weights, with provenance asserts ##
    m_orient = orientation_multiplicity(2, True)
    assert m_orient == 2                                  # wiring thm
    assert charge_trace == Fraction(8, 3)                 # forced vertex
    assert C2_26 == d10 * d11 == 6                        # branching rep
    W_fund = Fraction(1, d10)                             # projector lemma
    dim_J3 = 1 + n26
    assert dim_J3 == 27 == d11**3                         # q-face algebra
    R = dict(coupling=Fraction(m_orient),
             action=charge_trace * C2_26,                 # 16
             amplitude=W_fund * dim_J3)                   # 27/2
    report(f"  response weights (all wiring-layer or forced-vertex "
           f"data):")
    report(f"    coupling  = {R['coupling']}  (2-cycle orientation "
           "count, wiring theorem)")
    report(f"    action    = {R['action']}  = charge_trace x C2(26) "
           "(the alpha_G2->alpha_EM conversion vertex x")
    report(f"                the embedding channel's Casimir)")
    report(f"    amplitude = {R['amplitude']} = W(fund) x dim J3(O) "
           "(the B/A Schur projector x the q-face")
    report(f"                Jordan algebra, 27 = 1 + n26 = d11^3)")

    # # one depth-3 carrier, read at recursive depth 4 ##
    web = Web()
    inv0 = 2**9 / math.pi
    L = lambda s: (1.0 / (2.0 * math.pi * s["inv_alpha"]))**2

    # THE one loop object (base 0; its only edge is the loop amplitude
    # read self-consistently from the state)
    web["eq_loop"] = Ledger(
        "eq_loop", 0.0, "add", web_depth=3).echo(
        ["e<->q 2-cycle"], L, 3, "FORCED",
        "the wiring theorem's unique 2-cycle, L = (alpha/2pi)^2")

    # coupling node: depth-1 bridge self-echo + a READ of the loop
    web["inv_alpha"] = (
        Ledger("1/alpha", inv0, "add", web_depth=4)
        .echo(["bridge"], -inv0 * float(h_bridge) / (2 * math.pi), 1,
              "FORCED", "normalized bridge self-response (canon depth-1)")
        .echo(["read:eq_loop"],
              lambda s: -inv0 * float(R["coupling"]) * s["eq_loop"], 4,
              "FORCED", "response rule: orientation count x L",
              parent="eq_loop"))

    # action node: canon depth-1 vents + a READ of the loop
    delta_S = vertex_charge_trace * alpha_EM / (2 * math.pi)
    web["S_quark"] = (
        Ledger("S_quark", S_lepton, "add", web_depth=4)
        .echo(["F4-knot(26)"], -float(C2_26), 1, "FORCED",
              "Casimir vent (canon depth-1)")
        .echo(["EM(alg)"], delta_S, 1, "FORCED",
              "E8 current trace Tr_+(Q^2)=30 (derived in wiring.py)")
        .echo(["read:eq_loop"],
              lambda s: -float(R["action"]) * s["eq_loop"], 4,
              "FORCED", "response rule: charge_trace x C2(26) x L",
              parent="eq_loop"))

    # amplitude nodes: circulant base, canon depth-2 EM echo, and a
    # READ of the loop weighted by the knot's own gap
    theta = float(h10)
    gap = {k: 1.0 + BA_ratio * math.cos(theta + 2 * math.pi * k / 3)
           for k in range(3)}
    for k_idx, name in ((1, "e"), (2, "mu"), (0, "tau")):
        web[f"F_{name}"] = (
            Ledger(f"F_{name}", 1.0, "mul", web_depth=4)
            .echo(["bridge", "EM"],
                  lambda s: -1.0 / (2.0 * math.pi * s["inv_alpha"]), 2,
                  "FORCED", "EM vertex echo (canon depth-2)")
            .echo(["read:eq_loop"],
                  (lambda s, g=gap[k_idx]:
                   -float(R["amplitude"]) * s["eq_loop"] * g), 4,
                  "FORCED", "response rule: W(fund) x dim J3 x L x "
                            "Delta_k", parent="eq_loop"))
    web.solve()

    # # the fixed point must BE canon ##
    x = web.state["inv_alpha"]
    # anchor inversion with the regenerated vent factors
    F = {n: 0.5 * math.exp(-S_lepton) * gap[k]**2 * web.state[f"F_{n}"]
         for k, n in ((1, "e"), (2, "mu"), (0, "tau"))}
    M_Pl_MeV = M_E_ANCHOR_MEV / F["e"]
    got = dict(inv_alpha=x,
               m_mu=M_Pl_MeV * F["mu"],
               m_tau=M_Pl_MeV * F["tau"],
               v_EW=M_Pl_MeV / 1e3 * math.exp(-web.state["S_quark"]))
    tol = dict(inv_alpha=5e-10, m_mu=5e-8, m_tau=5e-7, v_EW=5e-7)
    report("")
    report("  fixed point of the constructive web (one depth-3 carrier; "
           "three depth-4 reads):")
    for key, target in TARGETS.items():
        ok = abs(got[key] - target) < tol[key]
        report(f"    {key:<9s} = {got[key]:<15.11g} canon "
               f"{target:<15.11g} {'REGENERATED' if ok else 'MISMATCH'}")
        assert ok, (key, got[key], target)

    # # equivalence to the LIVE canon web, node by node ##
    # (frozen literals above; the solved engine here, so the check
    #  cannot drift from canon even if the literals are re-frozen)
    from root import WEB as CANON
    assert abs(x - CANON.state["inv_alpha"]) < 1e-9
    assert abs(web.state["S_quark"] - CANON.state["S_quark"]) < 1e-12
    assert abs(M_Pl_MeV / CANON.state["M_Pl_MeV"] - 1.0) < 1e-9
    report("")
    report("  live-canon equivalence: inv_alpha, S_quark, M_Pl match "
           "root.py's solved web.")

    # # exclusivity: the pinning locks ##
    # (frozen separately in test_probes.py::test_vertex_pinning_locks;
    #  run here for the report values only, no duplicate freeze)
    try:
        from probes import vertex_pinning
    except ImportError:
        import vertex_pinning
    pin = vertex_pinning.run(report=lambda *a, **k: None)
    (s1, s2, s3) = pin["survivors"]
    n1, n2, n3 = pin["nearest_dead_sigma"]
    report("")
    report(f"  exclusivity (vertex_pinning): unique survivors "
           f"(2, 16, 27/2); nearest alternatives dead at")
    report(f"  ({n1:.0f}, {n2:.0f}, {n3:.0f}) sigma against "
           "(1/alpha Berkeley, v_EW/G_F, m_mu CODATA).")

    report("-" * 70)
    report("  VERDICT: one loop object, read through the supplied node "
           "responses, reproduces the")
    report("  canon fixed point; no per-node carrier declaration remains.")
    report("  This is equivalence. The unital/inverse-path normalization "
           "is derived in wiring_theorem_probe.py.")
    report("  Kill: this web failing any precision-block value, or "
           "a derivation forcing any")
    report("  response weight off its pinned value.")
    return dict(regenerated=got, responses={k: str(v)
                                            for k, v in R.items()})


if __name__ == "__main__":
    run()
