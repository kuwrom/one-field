"""The boundary factor k from the coupling variable the engine already reads.

THE QUESTION. With M_b^2 = 1/(k G), the factors with which G enters
general relativity are k = 1 (force law), 4 pi (Gauss), 8 pi (field
equations), 16 pi (action), 32 pi (graviton coupling). The wired gates
pass k = 1 alone (tests/audit/g_path_trials.py G10). Does the engine's
own coupling variable fix k, with no measured G?

PREMISES. P1 one channel coupling definition for every channel:
"closure probability times the normalized dimension one current
response" (interference/wiring.py), read for the EM channel as 1/alpha
and for the common mode as G M_b^2. P2 the engine's EM output matches
the measured fine structure constant, which is the Coulomb coefficient:
alpha = V r/(hbar c) between two unit charges. P3 the common mode's
coefficient in that variable is the Newton coefficient,
alpha_G = V r/(hbar c) = G m1 m2/(hbar c) between two closures.

STATEMENTS (the numbers below check each one).

  N1 (the variable). The engine's 1/alpha(0) = 137.036 is the force law
     coefficient. In the rationalized variable e^2 = 4 pi alpha it would
     be compared with 10.905, and in 1/(4 pi alpha) with 1722.0. The
     measured value picks the force law variable, by a factor 4 pi in
     either direction, with no measured G.
  N2 (the factor). One definition reads the common mode in the same
     variable, so G M_b^2 is the force law coefficient: k = 1. Any other
     GR factor k would read the common mode in a variable rationalized
     by k, and the EM channel would then give 1/alpha off by the same
     factor.
  N3 (what G sees). With k = 1, G comes out at -2.43 standard
     deviations. k = 4 pi would put it a factor 12.6 low.

RESULT. k = 1 is the convention already built into the engine's
coupling variable. The measured fine structure constant fixes that
variable, by a factor 4 pi, with no measured G. The boundary factor is
therefore not a choice separate from the coupling definition. The
reading is an identification (IH). Reading the common mode in the EM
channel's variable is premise P1, one definition for every channel.
The wired gates exclude 4 pi through a threshold (g_path_trials.py
G10). This probe gives a reason for k = 1 that uses no threshold.

Scope: the normalization of the coupling variable. The unit charge of
the common mode, M_b, comes from the anchor inversion.

Usage: python3 tests/probes/coupling_normalization.py
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402

ALPHA_CS = (137.035999046, 27e-9)
G_CODATA = (6.67430e-11, 0.00015e-11)
FACTORS = {"1": 1.0, "4pi": 4 * math.pi, "8pi": 8 * math.pi, "16pi": 16 * math.pi, "32pi": 32 * math.pi}


def run(report=print):
    out = {}

    # N1: the variable the EM output is compared in.
    x = r.inv_alpha_phys
    readings = {"force law alpha": ALPHA_CS[0], "rationalized e^2 = 4 pi alpha": ALPHA_CS[0] / (4 * math.pi),
                "1/(4 pi alpha)": ALPHA_CS[0] * 4 * math.pi}
    pulls = {k: (x - v) / ALPHA_CS[1] for k, v in readings.items()}
    assert abs(pulls["force law alpha"]) < 1
    assert min(abs(p) for k, p in pulls.items() if k != "force law alpha") > 1e8
    out.update(readings=readings, pulls=pulls)
    report("N1 the engine's 1/alpha(0) against the measured value in each variable: "
           + ", ".join(f"{k} {v:.3f} ({p:+.1e} sigma)" for k, (v, p) in
                       zip(readings, zip(readings.values(), pulls.values()))))

    # N2, N3: the common mode in each GR factor's variable.
    G = {k: r.G_PRED / f for k, f in FACTORS.items()}
    g_pulls = {k: (g - G_CODATA[0]) / G_CODATA[1] for k, g in G.items()}
    assert abs(g_pulls["1"] + 2.43) < 0.01 and abs(G["1"] / G["4pi"] - 4 * math.pi) < 1e-12
    out.update(G=G, g_pulls=g_pulls)
    report("N2 G = common_mode hbar c/(k M_b^2): " + ", ".join(f"k = {k}: {g:.4e}" for k, g in G.items()))
    report(f"N3 k = 1 at {g_pulls['1']:+.2f} sigma; k = 4 pi a factor {4 * math.pi:.1f} low")
    return out


if __name__ == "__main__":
    run()
