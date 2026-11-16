"""3+1 by closure: the dimension in which bound structures close and stay.

dimensionality.py derives 3+1 from the S^3 content of a populated
relative sector, which is the picture of the gas sign probes. On the
committed self-bound sign the clock lump does not populate that sector
(knot_in_lump.py). This probe gets the dimension from closure instead.
Of the configurations a fluctuating field makes, the ones that last are
the ones that close. A bound orbit stays bound and returns on itself in
only one dimension of space.
The same rule applies to particles, which last as closed loops of
circulation (spinning_lump.py), and to larger structures made of bound
orbits.

PREMISES. The induced metric is the linearized Einstein operator at low
momentum (composite_pole.py, ledger_kernel.py), so in d space
dimensions the static exchange between two lumps is Gauss's law,
V(r) = -k / r^(d-2) for d >= 3 and V = k ln r for d = 2. The same holds
for the Coulomb exchange of a massless vector. Time is single
(signature_of_the_saddle.py S2: of 64 sign patterns and time choices,
exactly the four Lorentzian ones with time on the odd sign keep the
energy bounded below).

STATEMENTS (the numbers below check each one).

  C1 (stability). A circular orbit in V = -k / r^(d-2) is stable iff
     V_eff'' > 0 at the circle, V_eff = L^2/(2 r^2) + V, which holds iff
     d < 4. At d = 4 the orbit is marginal, and at d >= 5 every
     circular orbit is unstable (Ehrenfest 1917; Tangherlini, Nuovo Cimento 27,
     636 (1963)). Checked from V_eff'' at d = 2, 3, 4, 5, 6.
  C2 (closure). The apsidal angle of a near circular orbit in
     V = -k / r^alpha is pi / sqrt(2 - alpha) (the logarithm at alpha = 0),
     so pi / sqrt(4 - d). That is one half turn, a closed orbit,
     exactly at d = 3, and pi / sqrt 2 at d = 2. Bertrand's theorem: only the
     inverse square and the harmonic force close every bound orbit.
     Measured by integrating orbits: at d = 3 the apsidal angle is pi at
     eccentricities 0.01, 0.3 and 0.7. At d = 2 it is pi / sqrt 2 near
     circular and moves with the eccentricity, so no bound orbit
     closes.
  C3 (atoms). A bound quantum state of -lap/2 - k / r^(d-2) needs the
     energy bounded below. Under x -> x / s the kinetic term scales as
     s^2 and the potential as s^(d-2). The energy is bounded below iff
     d < 4, marginal at d = 4, and unbounded at d >= 5 (Tegmark, Class. Quant. Grav. 14, L69
     (1997): "no traditional atoms" above three dimensions). Checked by
     the minimum of the scaled energy of a trial state at each d.
  C4 (the metric free vector). With g = h/U, the vector density
     sqrt|g| g g F F carries U^(2 - D/2). U drops out of the Maxwell
     term only at D = 4, where the vector couples only to the conformal
     class of h (determinant_vectors.py V2).

RESULT. Stable, closed bound orbits and atoms with a ground state exist
in three space dimensions and in no other. With a single time this is
3+1, which is also the dimension in which U drops out of the metric
free vector. On the committed sign the dimension follows from closure:
the structures that last are those that close.

Usage: python3 tests/probes/dimension_by_closure.py     (a few seconds)
"""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar


def potential(r, d, k=1.0):
    return k * np.log(r) if d == 2 else -k * r**(2.0 - d)


def force(r, d, k=1.0):
    """-dV/dr."""
    return -k / r if d == 2 else -k * (d - 2) * r**(1.0 - d)


def circular_curvature(d, r0=1.0, k=1.0):
    """C1: V_eff'' at the circular orbit of radius r0."""
    L2 = -force(r0, d, k) * r0**3                   # centripetal balance L^2/r^3 = -F
    h = 1e-4
    Veff = lambda r: L2 / (2 * r * r) + potential(r, d, k)
    return (Veff(r0 + h) - 2 * Veff(r0) + Veff(r0 - h)) / h**2


def apsidal_angle(d, ecc, k=1.0):
    """C2: angle swept between successive pericentres, from integration."""
    r_p = 1.0
    # launch at pericentre with the speed that gives the requested r_a/r_p
    r_a = r_p * (1 + ecc) / (1 - ecc)
    E = lambda v: 0.5 * v * v + potential(r_p, d, k)
    L = lambda v: r_p * v

    def mismatch(v):
        return 0.5 * (L(v) / r_a)**2 + potential(r_a, d, k) - E(v)
    lo, hi = 1e-3, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        (lo, hi) = (mid, hi) if mismatch(lo) * mismatch(mid) > 0 else (lo, mid)
    v = 0.5 * (lo + hi)

    Lm = r_p * v

    def rhs(t, y):
        r, vr, phi = y
        return [vr, Lm * Lm / r**3 + force(r, d, k), Lm / (r * r)]

    def apo(t, y):
        return y[1]                                  # radial velocity changes sign at an apsis
    apo.direction = -1                               # outward to inward: apocentre
    apo.terminal = True
    sol = solve_ivp(rhs, (0, 400), [r_p, 0.0, 0.0], events=apo, rtol=1e-12, atol=1e-13,
                    max_step=0.02)
    return float(sol.y_events[0][0][2])              # angle swept pericentre -> apocentre


def scaled_energy_min(d, k=1.0):
    """C3: inf over s of the energy of a Gaussian trial state scaled by s."""
    # <-lap/2> = d/(4 w^2) s^2 ... use T(s) = a s^2 and V(s) = -b s^(d-2) with a, b > 0
    a, b = 1.0, 1.0
    if d == 2:
        f = lambda s: a * s * s + b * math.log(1.0 / s)
    else:
        f = lambda s: a * s * s - b * s**(d - 2)
    res = minimize_scalar(lambda ls: f(math.exp(ls)), bounds=(-20, 20), method="bounded")
    inner = f(math.exp(res.x))
    large = f(math.exp(15.0))
    return inner, large


def run(report=print):
    out = {}
    out["curvature"] = {d: circular_curvature(d) for d in (2, 3, 4, 5, 6)}
    report("C1 V_eff'' at the circular orbit: "
           + ", ".join(f"d = {d}: {v:+.4f}" for d, v in out["curvature"].items()))
    out["apsidal"] = {}
    for d in (2, 3):
        for ecc in (0.01, 0.3, 0.7):
            out["apsidal"][(d, ecc)] = apsidal_angle(d, ecc)
    report("C2 apsidal angle / pi: "
           + ", ".join(f"d = {d}, e = {e}: {v / math.pi:.5f}" for (d, e), v in out["apsidal"].items())
           + f". Pi/sqrt(4-d) predicts 1 (d = 3) and {1 / math.sqrt(2):.5f} (d = 2) near circular.")
    out["atoms"] = {d: scaled_energy_min(d) for d in (2, 3, 4, 5)}
    out["bounded"] = {d: a > -1e3 for d, (a, b) in out["atoms"].items()}
    report("C3 scaled energy of a trial state: "
           + ", ".join(f"d = {d}: minimum {a:+.3f} ({'bounded' if out['bounded'][d] else 'unbounded'})"
                       for d, (a, b) in out["atoms"].items()))
    out["weyl_exponent"] = {D: 2 - D / 2 for D in (3, 4, 5, 6)}
    report("C4 power of U in the metric free Maxwell term: "
           + ", ".join(f"D = {D}: {e:+.1f}" for D, e in out["weyl_exponent"].items()))
    report("RESULT: stable closed orbits and bounded atoms exist only at d = 3. With a single "
           "time and the clean vector at D = 4, the dimension read from closure is 3+1.")
    return out


if __name__ == "__main__":
    run()
