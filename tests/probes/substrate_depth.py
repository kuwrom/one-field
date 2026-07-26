"""
substrate_depth.py: the substrate depth, measured through the
matter-stability channel.

SCOPE.  The one quantity the physical layer adds is the
interference weight of the medium at its density: the depth
rho xi^3, the coherent amplitude sharing one healing cell.  The code
derives the SHAPE of the self-interference (the coupling ratio
g_1/g_0 = 1/sqrt(2), nls_soliton.py) but not its strength; node
events, the only channel that unwinds a knot (knot_charge.py,
measured), are finite-action excursions off the fixed-norm manifold,
suppressed as exp(-c rho xi^3) with c an O(1) action coefficient
(the confined-phase entry test, confined_phase.py).  That weight is a measurement, not an
assumption.

THE MEASUREMENT, anchor-chained.  m_e, the one anchor,
gives M_Pl and with it the Planck attempt clock; the OBSERVED
STABILITY OF MATTER then bounds the medium's interference weight at
its density.  Nothing here is a second input to the engine: the
floor is read through the anchor's own chain, no output depends on
the height above it, and the height is visible through exactly one
channel, proton decay, which the web bets never occurs
(root.EchoTerm: no edge kind changes a knot's species).  A second
number could enter the framework only
by falsifying it.  The
attempt clock of a node event is the healing time, Planckian, since
the healing length is xi_0 = l_Pl/d10 (the Planck
identification, root.py C_EFF; this bound is the one
consumer sensitive to the O(1) factor, through xi^3), so a proton
that survives tau_p unwinds not, and

    c rho xi^3  >  ln(tau_p / t_Pl).

With the Super-K partial-lifetime bound tau_p > 1.6e34 yr (the same
bound the framework bets is never beaten) and t_Pl computed from the
framework's own derived M_Pl:

    per-proton reading:   c rho xi^3  >  ~ 1.96e2
    ensemble reading:     c rho xi^3  >  ~ 3.2e2
                          (no unwinding among ~1e80 baryons over a
                          Hubble time)

The deep-condensate requirement rho xi^3 >> 1 is therefore not a
premise: it is measured, with two decades of margin over unity, and
every future improvement of the proton-lifetime bound RAISES the
floor.  What remains is the bet that no derivation forces
the opposite (rho xi^3 ~ 1), which is the kill condition of
the confined-phase entry test (confined_phase.py).

STANDING OF THE PHYSICAL LAYER AFTER THIS READING.  It holds
no open premise: #1 depth = measured lower bound + its kill;
#4 winding protection = conditional derivation (the stabilizing
Faddeev-Skyrme quartic is constructive given the gap,
effective_action.py, and the gap is measured, knot_spectrum.py);
#2 (Lorentz), #5 (locked clocks) = untested predictions with stated
kills; #3 = the conjecture's scope, whose falsifiers are already in
PREDICTIONS (dark-matter nulls, Bell, EP).  Every entry is a
derivation, a measurement, or a bet.

Caveat, stated plainly: c is an O(1) action coefficient the
effective theory has not computed; the measured floor is on the
product c rho xi^3.  A computed c rescales the floor, not the
conclusion (two decades of margin).

ENDGAME (the grounded-strength reading, the gauge-fields-as-substrate-
excitations conjecture):
the algebra already fixes an inherent Planck-scale interference
strength, alpha_G2(M_Pl) = 1/(24 pi) (root.py, topological
bootstrap).  If the gauge sector is derived as substrate
excitations, that coupling becomes a function of substrate
parameters and 1/(24 pi) inverts into a DERIVED depth: the medium's
interference weight would be inherent to the minimal form, and the
contingency surface would collapse to existence itself plus the
names of units.

Usage: python3 substrate_depth.py
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.dirname(_here),
          os.path.join(os.path.dirname(os.path.dirname(_here)),
                       "interference"),
          os.path.join(_here, "interference")):
    if p not in sys.path:
        sys.path.insert(0, p)
import _bootstrap  # noqa: E402,F401

# observational constants single-homed in audit/reference.py
from root import SUPER_K

TAU_P_YR = SUPER_K['tau_p_yr']       # Super-K partial-lifetime bound
YR_S = SUPER_K['yr_s']               # Julian year, s
N_BARYON = SUPER_K['N_baryon']       # baryons in the observable universe
T_HUBBLE_S = SUPER_K['t_Hubble_s']   # ~13.8 Gyr


def run(report=print):
    from root import M_Pl_GeV, HBAR_SI, GEV_PER_J

    report("SUBSTRATE DEPTH: measured through the matter-stability "
           "channel")
    report("=" * 66)
    # Planck time from the framework's own derived M_Pl (m_e anchor)
    E_Pl_J = M_Pl_GeV / GEV_PER_J
    t_Pl = HBAR_SI / E_Pl_J
    assert 5.0e-44 < t_Pl < 5.8e-44
    tau_p = TAU_P_YR * YR_S
    per_proton = math.log(tau_p / t_Pl)
    ensemble = math.log(N_BARYON * T_HUBBLE_S / t_Pl)
    report(f"  attempt clock: t_Pl = hbar/M_Pl = {t_Pl:.3e} s "
           f"(M_Pl derived from the m_e anchor)")
    report(f"  per-proton (Super-K tau_p > 1.6e34 yr): "
           f"c rho xi^3 > ln(tau_p/t_Pl) = {per_proton:.1f}")
    report(f"  ensemble (~1e80 baryons, Hubble time):  "
           f"c rho xi^3 > {ensemble:.1f}")
    assert 190.0 < per_proton < 200.0
    assert 315.0 < ensemble < 330.0
    report("-" * 66)
    report("  the deep-condensate requirement rho xi^3 >> 1 is "
           "MEASURED (two decades over unity),")
    report("  not assumed; every improved proton-lifetime bound "
           "raises the floor.  The kill:")
    report("  a derivation forcing rho xi^3 ~ 1 "
           "falsifies the conjecture layer.")
    return dict(t_Pl=t_Pl, per_proton=per_proton, ensemble=ensemble)


if __name__ == "__main__":
    run()
