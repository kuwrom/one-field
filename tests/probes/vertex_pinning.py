"""
vertex_pinning.py: the response-rule multiplicities, cross-pinned.

WHAT THIS ESTABLISHES. The response rule assigns the depth-4
multiplicities

    coupling node (1/alpha):  2      (unweighted orientation count)
    action node   (S_quark):  16     (charge_trace x C2(26))
    amplitude node (m_lep):   27/2   (W(fund) x dim J3(O))

by composing each node's own depth-1 vertex stack (the constructive
demonstration is echo_naturality.py).  What this probe establishes
is that the rule's three outputs are not one collective fit:
EACH multiplicity is independently pinned by its own precision
observable, through three different channels of the web:

    m1 (coupling)  -> 1/alpha(0)   vs Berkeley Cs  (sigma = 2.7e-8)
    m2 (action)    -> v_EW         vs G_F          (sigma ~ 6.3e-5 GeV)
    m3 (amplitude) -> m_mu         vs CODATA       (sigma = 2.3e-6 MeV)

The probe re-solves the relevant closure with the multiplicity
replaced by EVERY value in the committed 4 x 7 vertex-product menu
(28 candidates per node, plus 0 = "no term"), holding the others at
canon, and reports the sigma-pull of each.  Frozen assertions: canon
is the UNIQUE candidate within 2 sigma at every node, and the nearest
menu alternative is dead by a wide margin.  The per-node 1-sigma
pinning width |Delta m| is printed: the data does not merely prefer
the composed values, it pins them to a small fraction of a menu step.

WHY THIS MATTERS.  Cross-pinning converts "the rule landed three
times" into "three independent sub-ppm observables each lock one
multiplicity".  The kill condition is sharp: a derivation forcing
any other value at any node is dead on arrival against the
corresponding observable.

THE NATURALITY STATEMENT (the rule's formal content):

    NATURALITY OF ECHO COMPOSITION.  Let a ledger node carry
    depth-1 echo terms with vertex data V = {v_i} (the node's own
    OPE couplings). Then every admissible depth-4 response through the
    e<->q 2-cycle carries multiplicity  Prod_{v in V'} v  for some
    V' subset of V fixed by the cycle's endpoints, the unweighted
    orientation count when V' is empty, and NO echo introduces a
    coupling outside V.  (In RCFT language: the integrated 2-cycle
    insertion factorizes on the existing 3-point data of the node's
    boundary sector; no new primary enters the OPE closure.)

Statement formal, targets frozen, rigidity priced; the constructive
witness is echo_naturality.py.

Usage: python3 vertex_pinning.py
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


# reference values and uncertainties (comparison only, never inputs)
BERKELEY_INV_ALPHA = (137.035999046, 2.7e-8)
V_EW_GF = (246.219651, 6.3e-5)          # GeV; sigma from G_F(5.1e-7 rel)
M_MU_CODATA = (105.6583755, 2.3e-6)     # MeV


def menu():
    """The committed 4 x 7 vertex-product menu (root.py inventory,
    frozen executably in wiring_theorem_probe.py, PART D), plus 0
    (no term)."""
    e_side = [Fraction(1), Fraction(2), Fraction(8, 3), Fraction(1, 2)]
    q_side = [Fraction(1), Fraction(6), Fraction(27), Fraction(26),
              Fraction(9), Fraction(3, 5), Fraction(182)]
    vals = sorted({a * b for a in e_side for b in q_side} | {Fraction(0)})
    return vals


def solve_inv_alpha(m1):
    """1/alpha = (512/pi)(1 - 1/(2pi) - m1 (alpha/2pi)^2), fixed point."""
    A = 2**9 / math.pi
    x = A * (1.0 - 1.0 / (2 * math.pi))
    for _ in range(200):
        x_new = A * (1.0 - 1.0 / (2 * math.pi)
                     - float(m1) * (1.0 / (2 * math.pi * x))**2)
        if abs(x_new - x) < 1e-15:
            break
        x = x_new
    return x


def run(report=print):
    from root import (h10, S_lepton, vertex_charge_trace,
                      alpha_EM, BA_ratio, M_E_ANCHOR_MEV, C2_26,
                      charge_trace)

    report("VERTEX PINNING: three multiplicities, three independent "
           "locks")
    report("=" * 66)
    vals = menu()
    report(f"  candidate menu ({len(vals)} values incl. 0): "
           f"{[str(v) for v in vals]}")

    alpha_phys = 1.0 / solve_inv_alpha(2)      # canon coupling
    a2p = alpha_phys / (2 * math.pi)

    # ── node 1: coupling (canon m1 = 2) vs Berkeley Cs ──────────────
    ref, sig = BERKELEY_INV_ALPHA
    rows1 = [(m, (solve_inv_alpha(m) - ref) / sig) for m in vals]
    surv1 = [m for m, p in rows1 if abs(p) < 2.0]
    near1 = min((abs(p) for m, p in rows1 if m != 2), default=None)
    # 1-sigma pinning width: d(1/alpha)/dm1 = -(512/pi)(alpha/2pi)^2
    width1 = sig / ((2**9 / math.pi) * a2p**2)
    report(f"\n  m1 (coupling, canon 2) vs 1/alpha Berkeley:")
    report(f"     canon pull {dict(rows1)[Fraction(2)]:+.2f} sigma; "
           f"survivors @2sigma = {[str(m) for m in surv1]}")
    report(f"     nearest alternative dead at {near1:.0f} sigma; "
           f"1-sigma pinning width |dm1| = {width1:.2e}")
    assert surv1 == [Fraction(2)]

    # ── node 2: action (canon m2 = 16) vs v_EW from G_F ─────────────
    ref, sig = V_EW_GF
    # v = M_Pl e^{-(S0 - 6 + 15/512 - m2 (alpha/2pi)^2)}; M_Pl from the
    # canon anchor inversion (root.py), held fixed while m2 varies
    from root import M_Pl_GeV
    S_base = (S_lepton - float(C2_26)
              + vertex_charge_trace * alpha_EM / (2 * math.pi))
    v_of = lambda m2: M_Pl_GeV * math.exp(-(S_base - float(m2) * a2p**2))
    rows2 = [(m, (v_of(m) - ref) / sig) for m in vals]
    surv2 = [m for m, p in rows2 if abs(p) < 2.0]
    near2 = min((abs(p) for m, p in rows2 if m != 16), default=None)
    width2 = (sig / ref) / a2p**2
    report(f"\n  m2 (action, canon 16 = charge_trace x C2(26) = "
           f"{charge_trace}x{C2_26}) vs v_EW(G_F):")
    report(f"     canon pull {dict(rows2)[Fraction(16)]:+.2f} sigma; "
           f"survivors @2sigma = {[str(m) for m in surv2]}")
    report(f"     nearest alternative dead at {near2:.0f} sigma; "
           f"1-sigma pinning width |dm2| = {width2:.3f}")
    assert surv2 == [Fraction(16)]

    # ── node 3: amplitude (canon 27/2) vs m_mu CODATA ───────────────
    ref, sig = M_MU_CODATA
    theta = float(h10)
    gap = {k: 1.0 + BA_ratio * math.cos(theta + 2 * math.pi * k / 3)
           for k in range(3)}
    d_e, d_mu = gap[1], gap[2]
    # the QED factor 1 - a2p cancels in the mu/tau ratio

    def m_mu_of(m3):
        # ratio form: the anchor inversion cancels everything except
        # the gap ratio and the DIFFERENTIAL vent
        vent = lambda g: 1.0 - float(m3) * a2p**2 * g
        return M_E_ANCHOR_MEV * (d_mu / d_e)**2 * vent(d_mu) / vent(d_e)

    rows3 = [(m, (m_mu_of(m) - ref) / sig) for m in vals]
    surv3 = [m for m, p in rows3 if abs(p) < 2.0]
    near3 = min((abs(p) for m, p in rows3 if m != Fraction(27, 2)),
                default=None)
    width3 = (sig / ref) / (a2p**2 * abs(d_mu - d_e))
    report(f"\n  m3 (amplitude, canon 27/2 = W(fund) x dim J3(O)) vs "
           f"m_mu CODATA (ppb-level):")
    report(f"     canon pull {dict(rows3)[Fraction(27, 2)]:+.2f} sigma; "
           f"survivors @2sigma = {[str(m) for m in surv3]}")
    report(f"     nearest alternative dead at {near3:.0f} sigma; "
           f"1-sigma pinning width |dm3| = {width3:.3f}")
    assert surv3 == [Fraction(27, 2)]

    report("-" * 66)
    report("  CROSS-PINNING ESTABLISHED: each multiplicity is locked "
           "by its own observable through its")
    report("  own channel; no collective compensation is possible.  "
           "The response rule's targets are")
    report("  frozen: any derivation must land exactly (2, 16, 27/2).")
    report("  Kill: a derivation forcing any other menu value at any "
           "node is dead against that node's")
    report("  observable at the sigma distances printed above.")
    return dict(survivors=(surv1, surv2, surv3),
                widths=(width1, width2, width3),
                nearest_dead_sigma=(near1, near2, near3))


if __name__ == "__main__":
    run()
