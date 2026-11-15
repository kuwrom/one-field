"""The logarithmic class of the induced action, compared with data.

THE QUESTION. regulator_free_ratios.py treats the vacuum class and the
dimension two class of the induced action. The logarithmic class holds
the dimension four operators: C^2, the Euler density, R^2, and the
cross term m^2 R that makes the Einstein stiffness depend on the masses
of the rows. Its coefficients are the a2 heat kernel coefficients. They
are the one part of the induced action that no regulator affects. What do they say
about the stiffness when matter accumulates, and about the extra modes
of the induced metric?

PREMISES. The rows of interference/gravity.py as counted in
ledger_kernel.py (189 real scalars, minimally coupled,
18 massless vectors, 45 Weyl fermions), one loop, the proper time
cutoff t0 = d11 l_Pl^2 for the Einstein term, the diagnostic's cutoff
(sakharov_normalisation.py FACT 1, FACT 7). The logarithmic moment
M_1 = ln(1/(mu^2 t0)) at the scale mu where an operator is evaluated.

COEFFICIENTS (units 1/(16 pi^2), a2 ⊃ c C^2 - a E + d R^2 + k m^2 R per
field, fermions with their loop sign once; Vassilevich hep-th/0306138
Table 1 and Eq. 4.28, Visser gr-qc/0204062 Eqs. 14-16, Duff
hep-th/9308075 Eqs. 30-31, Herranen et al. arXiv:1407.3141 Table I):
    real scalar      c = 1/120  a = 1/360     d = (xi - 1/6)^2/2  k = xi - 1/6
    Weyl fermion     c = 1/40   a = 11/720    d = 0               k = 1/12 per state
    gauge vector     c = 1/10   a = 31/180    d = 0
    Proca            c = 13/120 a = 7/40      d = 1/72            k = 1/2

STATEMENTS (the numbers below check each one).

  L1 (the sums). Over the ledger's rows c = 9/2, a = 69/16,
     d = 21/8. Massless vectors and fermions are conformal: R^2 comes
     from the scalars alone.
  L2 (accumulation). A row of mass m contributes a1 (M_2 - m^2 M_1) to
     the stiffness, so a change of the rows' masses moves G by
     dG/G = -sum n a1 dm^2 M_1 t0 / Sigma. Switching every Standard
     Model mass off, the largest change matter can make, moves G by
     1.2e-33. Matter that accumulates sources the metric. It
     does not change the stiffness. Lunar laser ranging bounds the
     drift at 7.6e-14 per year (Hofmann and Mueller, cited in
     arXiv:1812.11181).
  L3 (the scalaron). With L ⊃ (M^2/2) R + (M_1/32 pi^2)(c C^2 + d R^2),
     the R^2 term carries a spin-0 mode of mass
     m0^2 = (8 pi^2/3) M^2/(d M_1), M the reduced Planck mass, healthy
     sign. For M_1 from 1 to 100 it lies between 0.32 M and 3.2 M: no
     fifth force at any distance tested (Eot-Wash excludes
     gravitational strength Yukawa ranges above 39 microns,
     arXiv:2002.11761). Inflation driven by R^2 needs
     m0 = 1.3e-5 M (Planck 2018 normalization): the induced R^2 is
     short by a factor d M_1 / 1.6e11, at least 5.9e8. The induced action
     does not carry Starobinsky inflation.
  L4 (the spin-2 mode). The C^2 term has the sign of every matter row:
     the massive spin-2 mode is a tachyonic ghost,
     m2^2 = -8 pi^2 M^2/(c M_1), the semiclassical instability of
     Horowitz and Wald. In cutoff units |m2^2| t0 = Sigma/(2 c M_1) and
     m0^2 t0 = Sigma/(6 d M_1). Both are O(1)/M_1, at the cutoff. There
     the local expansion stops and only the full kernel applies, and
     the full kernel has no zero on the Euclidean axis from q^2 t0 = 1e-3
     to 1e3 (ledger_kernel.py R3). The
     ratio m0^2/|m2^2| = c/(3d) is free of the regulator and of M_1.

RESULT. The logarithmic class is the regulator free part of the
induced action. It meets data only through structural statements: the
stiffness does not move when matter accumulates, and both extra modes
of the induced metric sit at the cutoff. The class carries no closure
condition of its own. The effect of accumulated matter on gravity is
carried by the sources and by the vacuum term.

Usage: python3 tests/probes/log_class.py     (about half a minute
standalone: the engine runs once)
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root  # noqa: E402
import composite_metric as cm  # noqa: E402
import ledger_kernel as lk  # noqa: E402

PER_FIELD = {                                   # (c, a, d(xi)) per field, units 1/(16 pi^2)
    "scalar": (1 / 120, 1 / 360, lambda xi: 0.5 * (xi - 1 / 6)**2),
    "dirac": (1 / 20, 11 / 360, lambda xi: 0.0),
    "vector": (1 / 10, 31 / 180, lambda xi: 0.0),
}
SM_MASSES = {"top": (173.0, 4 / 12), "W": (80.4, 2 * 3 / 6),     # (GeV, |n a1| of the row)
             "Z": (91.2, 3 / 6)}
LLR_GDOT = 7.6e-14                              # per year, one sigma
EOT_WASH_M = 0.1973e-6 / 39e-6                  # eV: hbar c / 39 microns
STAROBINSKY = 1.3e-5                            # m0 / reduced Planck mass


def sums():
    c = a = d = 0.0
    for kind, n, xi in lk.ledger_rows():
        ci, ai, di = PER_FIELD[kind]
        c, a, d = c + n * ci, a + n * ai, d + n * di(xi)
    return c, a, d


def run(report=print):
    c, a, d = sums()
    assert abs(c - 9 / 2) < 1e-12 and abs(d - 2.625) < 1e-12
    sigma = cm.ledger()["sigma"]
    M_Pl = root.M_Pl_GeV
    t0 = root.d11 / M_Pl**2                     # GeV^-2
    report(f"L1 sums over the ledger: c = {c:.4f}, a = {a:.4f}, d = {d:.4f}; Sigma = {sigma:.4f}.")

    shift = sum(w * m * m * math.log(1 / (m * m * t0)) for m, w in SM_MASSES.values()) * t0 / sigma
    assert 1e-33 < shift < 2e-33
    report(f"L2 accumulation: switching every Standard Model mass off moves G by {shift:.1e}; "
           f"lunar ranging allows {LLR_GDOT:.1e} per year.")

    M = M_Pl / math.sqrt(8 * math.pi)
    m0 = {M1: math.sqrt(8 * math.pi**2 / (3 * d * M1)) for M1 in (1.0, 100.0)}
    short = 8 * math.pi**2 / (3 * STAROBINSKY**2) / (d * 100.0)
    assert min(m0.values()) * M * 1e9 > EOT_WASH_M and short > 5.9e8
    report(f"L3 scalaron: m0/M = {m0[1.0]:.2f} (M_1 = 1) to {m0[100.0]:.2f} (M_1 = 100), "
           f"M = {M:.3e} GeV; Eot-Wash needs m0 > {EOT_WASH_M:.1e} eV. Starobinsky needs "
           f"m0/M = {STAROBINSKY:.1e}: the induced R^2 is short by {short:.1e} or more.")

    ratio = c / (3 * d)
    report(f"L4 spin-2 mode: tachyonic ghost; |m2^2| t0 = {sigma / (2 * c):.3f}/M_1, "
           f"m0^2 t0 = {sigma / (6 * d):.3f}/M_1, both at the cutoff; "
           f"m0^2/|m2^2| = c/(3d) = {ratio:.4f}, free of the regulator.")
    report("RESULT: the logarithmic class is regulator free and carries no closure condition "
           "of its own: the stiffness does not move when matter accumulates, and both extra "
           "modes of the induced metric sit at the cutoff.")
    return dict(c=c, a=a, d=d, sigma=sigma, accumulation=shift, m0_over_M=m0,
                starobinsky_shortfall=short, m0_over_m2_sq=ratio)


if __name__ == "__main__":
    run()
