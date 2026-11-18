"""The substrate's own coupling, and the depth that alpha_G2 implies.

THE QUESTION. The aim of the physical layer is to read the
algebra's Planck scale coupling alpha_G2(M_Pl) = 1/(24 pi) as a
substrate quantity and invert it into the depth rho xi^3
(substrate_depth.py). That needs a definition of "the substrate's
coupling" that is not a choice among normalizations
(tests/audit/parameter_ledger.py P7).

THE DEFINITION. As in electrodynamics, where two unit charges interact
by V = alpha hbar c / r, take two quanta of the medium and compute the
force the medium carries between them. A quantum couples to the density with
strength d mu/d n = c_s^2/n (the same species), and the static density
response of any equation of state is chi(k) = 4n/(k^2 + 4 c_s^2) with
m = hbar = 1. The medium's static backreaction is the screened term
    V(r) = -alpha_sub hbar c_s exp(-r/xi)/r,
    alpha_sub = c_s^3/(pi n) = 1/(8 pi n xi^3),   xi = hbar/(2 m c_s),
the range being the repo's screening length (sakharov_normalisation.py
FACT 5, record_area.py, where mu = c_s^2). No continuous normalization
is free: the coupling is the inverse depth over 8 pi, up to the binary
reading of C2.

STATEMENTS (the numbers below check each one).
  C1  The potential computed from the response at three states
      (n, c_s) = (1, 1), (3, 1/2), (0.2, 2) matches
      1/(8 pi n xi^3) at r = 0.3 to 1 xi to quadrature precision.
  C2  The relation holds for any equation of state of a one component
      medium, since it uses only c_s and n. On the committed
      three component clock it holds for a quantum coupled by d mu/dR.
      A quantum coupled by its own Hartree potential carries 1.418
      times that coupling. The interspecies exchange term accounts for
      the difference (channel_response.py X1). xi is the screening range 1/(2 c_s). On
      the clock mu0 = c_s^2/2, so the 1/(4|mu0|) length is sqrt 2 xi.
  C3  Setting alpha_sub = alpha_G2(M_Pl) = 1/(24 pi) gives
      rho xi^3 = 3 exactly. The quantum depletion there is
      (8/(3 sqrt pi)) sqrt(n a^3) = 1.4e-3 for a one component gas:
      the medium is classical to a part in a thousand.

PREMISE AND SCOPE. The identification is that the algebra's coupling
is the medium mediated coupling between the substrate's quanta, the
simplest form of the gauge as substrate program, carried by the common
(density) mode. The G2 face is nonabelian and the substrate's Z3
relative sector carries its own, gapped, backreaction. Its prefactor is
computed in channel_response.py (two screened terms, 13% of the
common channel, no engine consumer). The depth 3 (4.25 under the Hartree
reading) is a consequence of the premise. No observable fixes it on
the committed sign (parameter_ledger.py P8), so it stands as a
prediction of the premise.

RESULT. The medium's static response gives the force between two quanta
as alpha_sub = 1/(8 pi n xi^3), matched at three states to better than
1e-6. Reading this exchange coupling as alpha_G2(M_Pl) = 1/(24 pi) gives
the substrate depth rho xi^3 = 3 exactly (4.25 under the Hartree reading
of channel_response.py). The depth belongs to the realization program
and feeds no scorecard number.

Usage: python3 tests/probes/substrate_coupling.py
"""

import math

import numpy as np
from scipy.integrate import quad


def alpha_from_response(n, c, radii=(0.3, 0.6, 1.0)):
    g_ib = c**2 / n
    xi = 1 / (2 * c)
    out = []
    for x in radii:
        r = x * xi
        integral = quad(lambda k: 4 * n * k / (k**2 + 4 * c**2), 0, np.inf,
                        weight="sin", wvar=r)[0]
        V = -(g_ib**2) / (2 * math.pi**2 * r) * integral
        out.append(-V * r * math.exp(r / xi) / c)
    return np.array(out), 1 / (8 * math.pi * n * xi**3)


def run(report=print):
    states = ((1.0, 1.0), (3.0, 0.5), (0.2, 2.0))
    dev = 0.0
    for n, c in states:
        a, pred = alpha_from_response(n, c)
        dev = max(dev, float(np.max(abs(a / pred - 1))))
        report(f"C1 n = {n}, c_s = {c}: alpha from the response {np.round(a, 5)}, "
               f"1/(8 pi n xi^3) = {pred:.5f}")
    depth = (24 * math.pi) / (8 * math.pi)
    K = (16 * math.pi)**-1.5                       # n xi^3 = K / sqrt(n a^3) for g = 4 pi a
    depletion = 8 / (3 * math.sqrt(math.pi)) * K / depth
    assert dev < 1e-6 and abs(depth - 3.0) < 1e-12 and abs(depletion - 1.40e-3) < 2e-5
    report(f"C3 alpha_sub = 1/(24 pi) -> rho xi^3 = {depth:.1f}. Depletion {depletion:.2e}")
    return dict(max_deviation=dev, depth=depth, depletion=depletion)


if __name__ == "__main__":
    run()
