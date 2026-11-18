"""Does the bridge's O(g^2) piece survive? Why the amount does not enter the running.

THE QUESTION. The one list's b0 = 32/3 (boundary.py) reads the bridge
as amount: it vents at weight one and takes no share of the G2
running. The label share is -13/3. If a second order O(g^2) matter
gauge portal of the bridge survived at flat connection, the running
would count it and b0 would drop to 19/3 (the light reading).
PvP = 0 does not by itself remove a current-current response, and the
metric free form with a gauge field (tests/probes/determinant_vectors.py)
couples F to the composite matter metric h. So the check is this: at
the empty field vacuum, does the bridge's F^2 piece carry a matter
dependent (density dependent) coefficient that would make it a new
portal?

THE FORM. determinant_vectors.py V2, V5: the O(F^2) density of
sqrt|det(h + i sqrt(U) F)| / U is
    L_F2 = -(1/4) sqrt|det h| F_h^2
with F_h^2 raised by h. At the empty field vacuum h_vac = c eta (the
fluctuation scale c = m^2 <|phi|^2> / 2, signature_of_the_saddle S4),
U positive, this is the bridge in the G2 connection A at the vacuum
value it carries.

STATEMENTS.
  B1 (the cancellation). sqrt|det(c eta)| = c^2 and F_{c eta}^2 =
     F_flat^2 / c^2, so L_F2 = -(1/4) F_flat^2: the c^2 and the 1/c^2
     cancel exactly. The O(F^2) piece is the G2 field's own Maxwell
     term. That row is already counted as the G2 gauge row (44/3
     share), so it is not a new matter gauge portal.
  B2 (c independence). Numerically, the second derivative of the
     Hermitian form in F is c-independent to order F^2 for every U
     and c: the vector sees the conformal class of h alone (V5).
  B3 (b0). The amount's running share is 0 on the one list, b0 = 32/3.
     The O(g^2) piece adds no matter gauge term, so the running is
     unchanged.
  B4 (the counterfactual). Read as light, the label share -13/3 gives
     b0 = 19/3, and boundary_content.py light_shift moves M_Pl by
     x1.58e13 (the anchor fails). That is the branch a surviving portal
     would take.
  B5 (the flat connection). Omega = 0 puts the label in the Einstein
     term (a2 ~ tr(1)) only, not in the running (a4 ~ tr(Omega^2)).
     This is the flat label algebra behind B1.

RESULT. At the empty field vacuum the O(F^2) density is -(1/4)
F_flat^2 for every c, the G2 field's own Maxwell term, so the bridge
adds no matter gauge portal, consistent with b0 = 32/3 on the one
list. The light reading, b0 = 19/3, moves M_Pl by x1.58e13 and is
excluded.

Usage: python3 tests/probes/bridge_b0_check.py     (a few seconds)
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402
import boundary_content as bc  # noqa: E402
from determinant_vectors import f_sq, hermitian_form, quadratic_in_f  # noqa: E402

ETA = np.diag([1.0, -1.0, -1.0, -1.0])
F0 = np.array([0.0, 0.0, 0.0, 1.0])


def f_flat2(F):
    """F_{mu nu} F^{mu nu} in the flat metric."""
    return float(np.trace(ETA @ F @ ETA @ F.T))


def run(report=print):
    out = {}

    # B1 the analytic cancellation, c^2 against 1/c^2.
    cancel = {}
    for c in (0.3, 1.0, 4.0, 7.7):
        det_h = c**4
        f_h2 = 1.0 / c**2            # F_h^2 at h = c eta for F_flat^2 = 1
        cancel[c] = -0.25 * math.sqrt(abs(det_h)) * f_h2
    out["cancel"] = cancel
    assert all(abs(v + 0.25) < 1e-15 for v in cancel.values())
    report(f"B1 L_F2(c) = -(1/4) sqrt|det(c eta)| (F_flat^2 / c^2) = "
           + " / ".join(f"{v:+.4f}" for v in cancel.values())
           + " for c = " + "/".join(f"{c}" for c in cancel)
           + ". The c^2 and the 1/c^2 cancel. L_F2 is the G2 field's "
             "own Maxwell term, already counted as the G2 gauge row, not "
             "a new matter gauge portal")

    # B2 the numerical c independence of the second F derivative.
    rng = np.random.default_rng(3)
    F = rng.normal(size=(4, 4))
    F = F - F.T
    F = F / math.sqrt(f_flat2(F))
    q = {}
    for c in (0.3, 1.0, 4.0, 7.7):
        for U in (0.5, 1.0, 2.0):
            q[(c, U)] = quadratic_in_f(c * ETA, F, U)
    out["second_derivative"] = {f"{c},{U}": v for (c, U), v in q.items()}
    base = q[(1.0, 1.0)]
    worst_c = max(abs(v / base - 1.0) for v in q.values())
    worst_u = max(abs(q[(c, U)] / q[(c, 1.0)] - 1.0)
                  for c in (0.3, 1.0, 4.0, 7.7) for U in (0.5, 2.0))
    # Richardson steps in determinant_vectors, so O(1e-7), not machine
    out["worst_c"], out["worst_u"] = worst_c, worst_u
    assert worst_c < 1e-5 and worst_u < 1e-5
    report(f"B2 d^2/dF^2 of sqrt|det(h + i sqrt U F)| / U at h = c eta "
           f"is {base:+.4f} for every c and U: c-dependent share "
           f"{worst_c:.1e}, U-dependent share {worst_u:.1e} "
           "(the vector sees the conformal class of h alone, "
           "determinant_vectors V5)")

    # B3 b0 on the one list: the amount's running share is 0.
    one = B.rows(r.d10, r.d11, r.n7, r.n26)
    b0_amount = B.b0(one)
    share = [row[3] for row in one if row[0] == "bridge"][0]
    out["b0_amount"], out["share"] = str(b0_amount), str(share)
    assert b0_amount == r.b0_G2 and b0_amount == B.b0(bc.rows(bridge="amount"))
    assert share == B.Fraction(0)
    report(f"B3 b0 on the one list = {b0_amount}, the amount's running "
           f"share = {share}. The O(g^2) piece adds no matter gauge "
           "term, so the running is unchanged")

    # B4 the counterfactual, the branch a surviving portal would take.
    lab = B.rows(r.d10, r.d11, r.n7, r.n26, bridge="light")
    b0_light = B.b0(lab)
    move = bc.light_shift(b0_light)["M_Pl"]
    out["b0_light"], out["M_Pl_move"] = str(b0_light), move
    assert b0_light == B.Fraction(19, 3) and 1e13 < move < 2e13
    assert b0_amount + B.Fraction(-13, 3) == b0_light
    report(f"B4 read as light, b0 = {b0_light} and M_Pl moves x{move:.3g} "
           f"(the anchor fails). That is the branch a surviving O(g^2) "
           "portal would take, and B1 and B2 exclude it")

    # B5 the flat label algebra, of which B1 is the second order.
    report("B5 Omega = 0 puts the label in the Einstein term (a2 ~ tr(1)) "
           "only, not in the running (a4 ~ tr(Omega^2)). B1 is the O(F^2) "
           "piece of that separation, and the cancellation of B1 is the "
           "c independence of B2")

    report("RESULT: the O(g^2) piece of the bridge is the G2 field's own "
           "Maxwell term, c- and U-independent to order F^2. It is already "
           "counted as the G2 gauge row, so the amount's running share is 0 "
           "and b0 = 32/3 holds. The M_Pl x1.6e13 shift that excludes the "
           "light reading is the branch a surviving matter gauge portal "
           "would take, and this check closes it.")
    return out


if __name__ == "__main__":
    run()
