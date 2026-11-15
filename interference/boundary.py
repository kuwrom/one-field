"""
boundary.py: the rows at the boundary M_b. One list serves both the running and the stiffness.

The transmutation (root.py) runs alpha_G2 = 1/(24 pi) down from the
boundary M_b, using the G2-charged rows that are light there. The
induced stiffness (gravity.py) sums the heat kernel weight of every row
that propagates at the regulator. Both use the list below, so a row
cannot contribute to the running without also venting.

How each row is read at M_b:
  light      propagates below M_b. It runs, and it vents with weight 1.
  amount     lies past the searchlight's horizon. It vents with weight 1,
             is minimally coupled, and does not run.
  boundary   has mass M_b. It leaves the running at M_b and vents
             with the heavy row weight E2(m^2 t0) of the proper time
             cutoff (tests/probes/one_list_kernel.py K1).
  composite  is bound below the regulator, so it is not a row (R6 of
             tests/probes/regulator_free_ratios.py).
  absent     has no spacetime row at all.

  row                states                  share of b0 if light    reading
  G2 gauge           2 dim(G2) = 28 vector   (11/3) h(G2) = 44/3     light
  electroweak gauge  2 (d11 + 1) = 8 vector  0                       light
  fermions           90 Weyl                 -(4/3) d11 T(7) = -4    light
  G2 scalars         n7 = 7 real             -(1/6) T(7)             composite
  bridge             n7 n26 = 182 real       none                    amount

The G2 scalars are the field that breaks G2: the condensate at
v_break = v_EW/sqrt 2 (couplings.py). Like the Higgs, they are not a
row at the regulator.

THE BRIDGE.  "Given the norm multiplicative ledger gate, n_7 = 7 is the
maximum number of imaginary directions this octonionic searchlight
distinguishes. Whatever folds beyond it is read as amount without
identity" (interference/__init__.py, THE HORIZON).  That is the bridge,
and n7 n26 = 182 is the searchlight's count of it.  It does not run. A
share of b0 is a Dynkin index, which is a label, and the amount carries
no label. It vents with weight 1: each channel has the minimal real
scalar weight 1/6 and is light at the regulator.

THE READING.  The list can read the bridge in five ways.  Read as light,
with its label as a running charge, it moves b0 away from the value the
anchor requires and moves M_Pl by 1.6e13. Read as absent or composite,
it leaves the stiffness at -4.5 and gate S fails. Read at the boundary's
mass, the sign of the stiffness depends on the spread of the regulator.
Read as amount, it keeps b0 and passes S and M at every contact value
(tests/probes/horizon_amount.py A5).  The amount is the only reading
that closes the running and both gates together.

THE AMOUNT'S COUPLING.  One covariant operator gives a row both its
share of b0 and its stiffness (Vassilevich arXiv:hep-th/0306138
Eqs. 4.27, 4.28).  The amount does not run, so its operator carries no
G2 connection. The coupling xi = alpha_G2 E[v^2] h_bridge = 1/(48 pi)
is a G2 self-coupling, so it has no vertex to come from. The amount is
therefore minimally coupled, with weight 1/6 (xi_of).  No gate depends
on this value: S and M hold at every contact for -0.178 < xi < 0.107
(tests/probes/horizon_amount.py A2, A4).  The readings bridge="light"
and bridge="boundary" of rows() keep the labelled bridge available for
comparison (tests/probes/boundary_content.py).

With G anchored, the regulator is t0 = Sigma/(2 pi) M_b^-2 (root.py), so
a list's stiffness is solved together with its regulator. The regulator
lies inside the boundary, where the list holds, exactly when Sigma >= 2 pi.
"""

import math
from fractions import Fraction

from scipy.optimize import brentq
from scipy.special import expn

T_7 = Fraction(1)                         # Dynkin index of the G2 fundamental
# Vector contact coefficient, the Ricci coefficient of the gauge operator's
# Weitzenbock term: +1 for the Hodge Laplacian (Bochner), the gauge field's
# Feynman gauge operator (tests/probes/gravity_premises.py Q2). The values
# 0 and -1 are also evaluated, to show that no verdict depends on this choice.
CONTACTS = (1.0, 0.0, -1.0)


def rows(d10, d11, n7, n26, bridge="amount", scalars="composite"):
    """The list at M_b: (name, states, kind, share of b0 when light, reading)."""
    label = -Fraction(1, 6) * n26 * T_7   # the bridge's searchlight label, taken as a running charge
    return [
        ("G2 gauge", 2 * (2 * n7), "vector", Fraction(11, 3) * d10**2, "light"),
        ("electroweak gauge", 2 * (d11 + 1), "vector", Fraction(0), "light"),
        ("fermions", d10 * d11 * (d11 * d10**2 + d10 + 1), "weyl",
         -Fraction(4, 3) * d11 * T_7, "light"),
        ("G2 scalars", n7, "scalar", -Fraction(1, 6) * T_7, scalars),
        ("bridge", n7 * n26, "bridge", Fraction(0) if bridge == "amount" else label, bridge),
    ]


def b0(table):
    """One loop coefficient of the G2 running: the light rows' shares."""
    return sum((b for _, _, _, b, reading in table if reading == "light"), Fraction(0))


def a1(kind, contact, xi):
    """Per state heat kernel weight (ledger_kernel.py R2)."""
    return {"vector": 1 / 6 - contact / 2, "weyl": 1 / 12, "scalar": 1 / 6, "bridge": 1 / 6 - xi}[kind]


def xi_of(reading, xi):
    """The curvature coupling of a bridge row: the label's xi when labelled, zero as amount.

    One covariant operator gives a row both its share of b0 and its
    stiffness (THE AMOUNT'S COUPLING). A row with no share has no G2
    vertex to carry xi = alpha_G2 E[v^2] h_bridge, so the amount is
    minimally coupled. The labelled readings keep the label's xi.
    """
    return 0.0 if reading == "amount" else xi


def sigma(table, t0, contact, xi, m_over_Mb=1.0, heavy=None):
    """Stiffness with the regulator at t0 (in M_b^-2). Row weights: light and amount 1, boundary E2(m^2 t0), others 0."""
    wb = float(expn(2, m_over_Mb**2 * t0)) if heavy is None else heavy
    w = {"light": 1.0, "amount": 1.0, "boundary": wb, "composite": 0.0, "absent": 0.0}
    return sum(n * a1(kind, contact, xi_of(reading, xi)) * w[reading] for _, n, kind, _, reading in table)


def _attractive_boundary_rows(table, contact, xi):
    """Every boundary row adds positive weight, so sigma falls as t0 grows."""
    assert all(a1(kind, contact, xi_of(reading, xi)) > 0
               for _, _, kind, _, reading in table if reading == "boundary")


def anchored(table, contact, xi, m_over_Mb=1.0, extra=0.0):
    """Positive solutions of Sigma = sigma(t0 = Sigma/(2 pi)): the stiffness with its anchored regulator.

    sigma - Sigma falls strictly, so there is at most one solution. It
    lies below the stiffness at t0 = 0, where every row has weight 1.
    """
    _attractive_boundary_rows(table, contact, xi)
    f = lambda s: sigma(table, s / (2 * math.pi), contact, xi, m_over_Mb) + extra - s
    lo, hi = 1e-9, sigma(table, 0.0, contact, xi, m_over_Mb) + extra + 1.0
    return [brentq(f, lo, hi)] if f(lo) > 0 else []


def inside(table, contact, xi):
    """An attractive stiffness whose anchored regulator lies inside the boundary exists.

    As the regulator moves below the boundary, sigma falls and 2 pi t0
    rises. So a solution with t0 >= M_b^-2 exists exactly when the
    stiffness at t0 = M_b^-2 reaches 2 pi.
    """
    _attractive_boundary_rows(table, contact, xi)
    return sigma(table, 1.0, contact, xi) >= 2 * math.pi
