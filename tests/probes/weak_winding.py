"""The weak winding: what it does to the baryon and lepton labels.

THE QUESTION. The label theorem of root.EchoTerm conserves baryon and
lepton number on every edge of the web. The weak stage is chiral, and a
change of the weak winding is not an edge. What does one unit of it do
to the labels, to the proton, and to the neutrino mass term?

PREMISES. P1 the content, from the engine's dictionary (root.py):
d11 = 3 generations, each with a quark doublet in N_c = d11 colours
carrying baryon label 1/N_c and charges Q_u = d10/d11 and -Q_d =
-1/d11, and a lepton doublet carrying lepton label 1 and charges 0 and
-1. P2 the weak stage is chiral: "The chiral weak branch forbids those
masses until symmetry breaking and is therefore the branch that closes"
(tests/probes/zero_mode.py), so the doublets are the left handed halves
and the right handed halves do not see the winding. P3 the count ('t
Hooft, PRL 37, 8 (1976)): one unit of winding changes a label by the
sum of that label over the doublets, and the rate carries the factor
exp(-4 pi/alpha_W). Anomaly matching is forced
(tests/probes/bookkeeping.py). P4 the engine's own M_W and v_EW, so
alpha_W = M_W^2/(pi v_EW^2).

STATEMENTS (the numbers below check each one).

  W1 (the count). One unit of winding changes baryon number by
     d11 N_c (1/N_c) = 3 and lepton number by d11 = 3. Their difference
     does not change.
  W2 (the event is neutral). The electric charge changes by the sum of
     the doublets' mean charges, d11 [N_c (Q_u - Q_d) - 1]/2 = 0. With
     the dictionary's forms N_c (Q_u - Q_d) = d10 - 1, so the event is
     neutral for d10 = 2 alone on the menu d10 = 1 to 4.
  W3 (the proton). Baryon number moves in threes. From one proton the
     nearest reachable value is -2, two antibaryons, heavier than the
     proton. A single proton cannot decay through the winding, so
     root.py's statement stands: no decay at any sensitivity. The
     lightest system the winding can reach holds two baryons.
  W4 (the rate). alpha_W = 0.0339 from the engine's M_W and v_EW, and
     exp(-4 pi/alpha_W) = 10^-161. No measurement separates a web with
     the winding from one without it.
  W5 (the neutrino mass term). A doublet mass or a singlet mass needs
     the lepton label to change by 2 with the baryon label fixed. The
     winding moves both by 3, and no whole number of units gives 2 and
     0. Read with B - L in place of L, the label gate of
     neutrino_mass_term.py passes the same three structures, and the
     survivor is the same.
  W6 (the record). Edges conserve B and L. The weak winding moves them
     together. So no edge and no weak winding changes B - L. Summed
     over the left handed content, B - L is -1 with two right handed
     neutrinos and 0 with three. That sum is what a winding of the
     metric would change B - L by, and the weak winding does not see
     it.

RESULT. The label theorem is a statement about edges, and it stands.
The engine's chiral content adds one event that is not an edge: a unit
of weak winding moves baryon and lepton number together by d11 = 3. It
is neutral for the engine's d10 = 2, it cannot reach a single proton,
it gives no Majorana term, and its rate factor is 10^-161. No edge and
no weak winding changes B - L.

Usage: python3 tests/probes/weak_winding.py
"""

import contextlib
import io
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import couplings, higgs_ew, masses  # noqa: E402
from interference import root as r  # noqa: E402
import neutrino_mass_term as nmt  # noqa: E402


def doublets(d10=r.d10, d11=r.d11):
    """P1, P2: the left handed doublets as (copies, baryon label, lepton label, mean charge)."""
    f = {k: form(Fraction(d10), Fraction(d11)) for k, form in r.CHARGE_FORMS.items()}
    quark = (d11 * f["N_c"], 1 / f["N_c"], Fraction(0), (f["Q_u"] - f["Q_d"]) / 2)
    lepton = (Fraction(d11), Fraction(0), Fraction(1), Fraction(0 - 1, 2))
    return quark, lepton


def per_unit(d10=r.d10, d11=r.d11):
    """P3: the change of each label under one unit of winding."""
    rows = doublets(d10, d11)
    return {"B": sum(n * b for n, b, _, _ in rows), "L": sum(n * l for n, _, l, _ in rows),
            "Q": sum(n * q for n, _, _, q in rows)}


def alpha_weak():
    """P4: alpha_W from the engine's own M_W and v_EW."""
    with contextlib.redirect_stdout(io.StringIO()):
        R = r.derive()
        m = masses.derive(R)
        h = higgs_ew.derive_higgs_ew(m, couplings.derive(R, m))
    return h["M_W_derived"] ** 2 / (math.pi * m["v_EW_pred_GeV"] ** 2)


def run(report=print):
    out = {}

    # W1: the count.
    step = per_unit()
    assert step["B"] == step["L"] == r.d11 == 3 and step["B"] - step["L"] == 0
    out.update(step={k: int(v) for k, v in step.items()})
    report(f"W1 one unit of winding: B {step['B']:+}, L {step['L']:+}, B - L {step['B'] - step['L']:+}")

    # W2: the event is neutral, and only for the engine's d10.
    neutral = [d for d in range(1, 5) if per_unit(d10=d)["Q"] == 0]
    assert step["Q"] == 0 and neutral == [r.d10] == [2]
    assert r.N_c * (r.Q_u - r.Q_d) == r.d10 - 1
    out.update(neutral_d10=neutral)
    report(f"W2 charge change {step['Q']}. Neutral on the menu d10 = 1 to 4 for d10 = {neutral}")

    # W3: the proton.
    reach = sorted((1 - int(step["B"]) * k for k in range(-3, 4) if k), key=abs)
    assert abs(reach[0]) == 2 and all(abs(b) > 1 for b in reach)
    lightest = min(b for b in range(1, 6) if any(abs(b - int(step["B"]) * k) < b for k in (1, -1)))
    assert lightest == 2
    out.update(proton_nearest=reach[0], lightest_system=lightest)
    report(f"W3 from one proton the nearest baryon number is {reach[0]}: heavier than the proton, no decay. "
           f"The lightest system the winding can reach holds {lightest} baryons")

    # W4: the rate.
    a_w = alpha_weak()
    exponent = 4 * math.pi / a_w / math.log(10)
    assert abs(a_w - 0.0339) < 2e-4 and 160 < exponent < 162
    out.update(alpha_W=a_w, log10_rate=-exponent)
    report(f"W4 alpha_W = {a_w:.4f}, exp(-4 pi/alpha_W) = 10^{-exponent:.0f}")

    # W5: the neutrino mass term.
    majorana = [k for k in range(-5, 6) if int(step["L"]) * k == 2 and int(step["B"]) * k == 0]
    b_minus_l = {t: (lambda l, f=f: -f(l)) for t, f in nmt.CHANGE.items()}
    same = [on for on in nmt.structures()
            if [l for l in nmt.LABELS if all(b_minus_l[t](l) == 0 for t in on)]]
    passed = [on for on in nmt.structures() if nmt.conserving_labels(on)]
    assert not majorana and same == passed == [(), ("b",), ("a",)]
    out.update(majorana_units=majorana, label_gate_same=same == passed)
    report("W5 no number of units changes L by 2 with B fixed. The label gate reads the same with B - L")

    # W6: the record, and the sum of B - L over the left handed content.
    def content_sum(singlets):
        n_c, gens = r.N_c, r.d11
        quarks = gens * n_c * (2 - 1 - 1) / n_c          # doublet, and the conjugates of the two right handed quarks
        leptons = gens * (-2 + 1)                        # doublet, and the conjugate of the right handed charged lepton
        return quarks + leptons + singlets               # the conjugate of each right handed neutrino
    assert content_sum(nmt.N_SINGLETS) == -1 and content_sum(3) == 0
    out.update(b_minus_l_sum=int(content_sum(nmt.N_SINGLETS)))
    report(f"W6 no edge and no weak winding changes B - L. Its sum over the content is {int(content_sum(nmt.N_SINGLETS))} "
           "with two right handed neutrinos, 0 with three")

    report("RESULT: edges conserve B and L. One unit of weak winding moves both by 3, is neutral, cannot reach "
           "a single proton, and carries the factor 10^-161. No edge and no weak winding changes B - L.")
    return out


if __name__ == "__main__":
    run()
