"""Radiation from a source through the kernel's pole: the quadrupole power, with Newton's G.

THE QUESTION.  The induced kernel gives Newton's law between two lumps
(two_lumps.py), and "the spin-2 pole that carries Newton's law carries
radiation with the same G" (post_newtonian.py N3).  What does a moving
source radiate through that pole, in which polarizations, and with
what coefficient?

PREMISES.  P1 the exchange between conserved sources through the
kernel K = F2 x (spin 2) + F0 x (spin 0) (two_lumps.py P3): up to one
overall constant c it is
    B(T1, T2) / k^2,   B = T1.T2 - a T1 T2,   a = 1/3 - (F2/F0)/3,
with T the trace.  The kernel's F0/F2 = -2 gives a = 1/2.  P2 the
constant is fixed by the static limit of the same exchange: two
masses at rest give B = (1 - a) m1 m2, and Newton's law is
c (1 - a) = 4 pi G, so c = 8 pi G.  P3 the pole's absorptive part is
the energy carried off: a source oscillating at frequency w loses
    dP/dOmega = (c/8 pi) (w^2/pi) B(T*, T)
per solid angle, B on the cone k = w (1, n) (for c = 8 pi G the
linearized result of Weinberg, Gravitation and Cosmology, Sec. 10.4).
P4 the source is conserved, w T^0i = k_j T^ij and w T^00 = k_i T^0i,
and slow: T^ij(w) = -(w^2/2) Q^ij(w), Q the second moment of the
energy density.

STATEMENTS (the numbers below check each one).

  Q1 (two polarizations).  The kernel's own ratio F0/F2, taken from
     composite_pole.py at q = 1e-3 m, is -2 to 2e-6, and its departure
     from -2 falls as q^2.  At a = 1/2, for every direction and every
     conserved source, B equals the squared transverse traceless part
     of T^ij, to 1e-12. So there are two polarizations, the power is
     positive, and there is no scalar and no vector wave.
  Q2 (the angular integral).  For a slow source the integral of B
     over directions is (8 pi/5) (T*^ij T^ij - |T^kk|^2/3): only the
     traceless quadrupole radiates.  A breathing sphere, Q
     proportional to the unit matrix, radiates nothing.
  Q3 (the coefficient).  With P2's constant the power is
     (G/5) <Q'''_ij Q'''_ij>, traceless, with the G of two_lumps.py.
     For two masses on a circular orbit of separation d and angular
     velocity W it is (32/5) G mu^2 d^4 W^6 (Peters and Mathews,
     Phys. Rev. 131, 435).
  Q4 (the ratio is not a choice).  The kernel's q^2 term is local
     only at F0/F2 = -2 (composite_pole.py K9), so no kernel of the
     lane carries a third, scalar polarization.  The orbital decay of
     PSR B1913+16 is 0.9983 +- 0.0016 of the quadrupole formula's
     (Weisberg and Huang, arXiv:1606.02744), 1.1 standard deviations
     from the kernel's own power.
  Q5 (what fails: the source not conserved).  For a source whose time
     components are not tied to its stress by P4, B takes both signs
     over random draws, so no power can be defined.

RESULT.  Through the pole that carries Newton's law a conserved slow
source radiates two transverse traceless polarizations with the power
(G/5) <Q'''_ij Q'''_ij>. This is the quadrupole formula with Newton's
G, and the binary pulsar's orbital decay tests it.

Scope: linear order in the kernel, slow sources, weak fields. Strong
fields and a horizon's temperature are outside the lane.

Usage: python3 tests/probes/quadrupole_radiation.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import composite_pole as cp  # noqa: E402

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
PULSAR = (0.9983, 0.0016)                  # PSR B1913+16, observed over predicted decay (arXiv:1606.02744)


def a_of(f0_over_f2):
    """The trace coefficient of the exchange for a kernel ratio F0/F2."""
    return 1.0 / 3.0 - 1.0 / (3.0 * f0_over_f2)


def kernel_ratio():
    """F0/F2 of the induced kernel at small momentum (composite_pole.py)."""
    J, _, _ = cp.proper_time(1.0, 1.0)
    F2, F0, _ = cp.form_factors(np.array([1e-3, 0.0, 0.0, 0.0]), J)
    return float(np.real(F0 / F2))


def conserved(Tij, n):
    """The full T^mu nu on the cone k = w (1, n) from its stress, by conservation."""
    T = np.zeros((4, 4), dtype=complex)
    T[1:, 1:] = Tij
    T[0, 1:] = T[1:, 0] = Tij @ n
    T[0, 0] = n @ Tij @ n
    return T


def bracket(T, a):
    """B(T*, T) = T*^mn T_mn - a |T^m_m|^2."""
    low = ETA @ T @ ETA
    return float(np.real(np.sum(np.conj(T) * low) - a * abs(np.trace(ETA @ T))**2))


def tt_square(Tij, n):
    """Squared transverse traceless part of T^ij for the direction n."""
    P = np.eye(3) - np.outer(n, n)
    S = P @ Tij @ P
    S = S - 0.5 * P * np.trace(S)
    return float(np.real(np.sum(np.conj(S) * S)))


def directions(n_theta=24, n_phi=48):
    """Gauss-Legendre in cos theta, uniform in phi: unit vectors and weights summing to 4 pi."""
    x, w = np.polynomial.legendre.leggauss(n_theta)
    phi = 2 * math.pi * np.arange(n_phi) / n_phi
    out = []
    for c, wc in zip(x, w):
        s = math.sqrt(1 - c * c)
        out += [(np.array([s * math.cos(p), s * math.sin(p), c]), wc * 2 * math.pi / n_phi) for p in phi]
    return out


def angular(Tij, a, grid):
    """Integral of B over directions for a slow conserved source of stress Tij."""
    return sum(w * bracket(conserved(Tij, n), a) for n, w in grid)


def binary(mu=1.0, d=1.0, W=1.0):
    """Circular orbit: the frequency w = 2 W and the amplitude Q_w of Q(t) = Q_w exp(-i w t) + c.c."""
    return 2 * W, (mu * d * d / 4) * np.array([[1, 1j, 0], [1j, -1, 0], [0, 0, 0]])


def power(w, Qw, a, c_over_8piG, grid):
    """Power in units of G: (c/8 pi G) (w^2/pi) times the angular integral of B, with T^ij = -(w^2/2) Q^ij."""
    return c_over_8piG * (w * w / math.pi) * angular(-(w * w / 2) * Qw, a, grid)


def run(report=print):
    out = {}
    rng = np.random.default_rng(7)
    grid = directions()
    assert abs(sum(w for _, w in grid) - 4 * math.pi) < 1e-12

    def stress():
        m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        return m + m.T

    def direction():
        v = rng.normal(size=3)
        return v / np.linalg.norm(v)

    # Q1: the kernel's ratio, and two polarizations.
    ratio = kernel_ratio()
    a = a_of(ratio)
    assert abs(ratio + 2) < 3e-6 and abs(a - 0.5) < 1e-6
    a = 0.5
    static = np.zeros((4, 4), dtype=complex)
    static[0, 0] = 1.0
    c_over_8piG = 0.5 / bracket(static, a)                      # P2: c (1 - a) = 4 pi G
    err, least = 0.0, math.inf
    for _ in range(200):
        Tij, n = stress(), direction()
        b, tt = bracket(conserved(Tij, n), a), tt_square(Tij, n)
        err, least = max(err, abs(b - tt) / tt), min(least, b)
    assert err < 1e-12 and least > 0 and abs(c_over_8piG - 1) < 1e-12
    out.update(kernel_ratio=ratio, tt_error=err)
    report(f"Q1 two polarizations: the kernel's F0/F2 = {ratio:+.6f}; B equals the squared transverse traceless "
           f"part to {err:.0e} over 200 conserved sources and directions, positive in each")

    # Q2: the angular integral.
    err2 = 0.0
    for _ in range(20):
        Tij = stress()
        exact = (8 * math.pi / 5) * (np.sum(abs(Tij)**2) - abs(np.trace(Tij))**2 / 3)
        err2 = max(err2, abs(angular(Tij, a, grid) / exact - 1))
    breathing = angular(np.eye(3, dtype=complex), a, grid)
    assert err2 < 1e-12 and abs(breathing) < 1e-12
    out.update(angular_error=err2, breathing=breathing)
    report(f"Q2 angular integral: (8 pi/5)(T*T - |tr T|^2/3) to {err2:.0e}; a breathing sphere gives "
           f"{abs(breathing):.0e}")

    # Q3: the coefficient.
    err3 = 0.0
    for _ in range(20):
        Qw, w = stress(), float(rng.uniform(0.5, 2.0))
        third = 2 * w**6 * (np.sum(abs(Qw)**2) - abs(np.trace(Qw))**2 / 3)   # <Q''' Q'''>, traceless
        err3 = max(err3, abs(power(w, Qw, a, c_over_8piG, grid) / (third / 5) - 1))
    w, Qw = binary()
    peters = power(w, Qw, a, c_over_8piG, grid)
    assert err3 < 1e-12 and abs(peters - 32 / 5) < 1e-12
    out.update(quadrupole_error=err3, binary=peters)
    report(f"Q3 coefficient: the power is (G/5) <Q''' Q'''> to {err3:.0e}. The circular binary gives "
           f"{peters:.6f} G mu^2 d^4 W^6, the 32/5 of Peters and Mathews")

    # Q4: the ratio is local only at -2 (composite_pole.py K9). The kernel's own power against the binary pulsar.
    assert abs(1.0 - PULSAR[0]) < 2 * PULSAR[1]
    out.update(einstein_pull=(1.0 - PULSAR[0]) / PULSAR[1])
    report(f"Q4 the kernel's own power sits {out['einstein_pull']:+.1f} standard deviations from PSR B1913+16's "
           f"{PULSAR[0]}({1e4 * PULSAR[1]:.0f})")

    # Q5: a source that is not conserved.
    signs = []
    for _ in range(200):
        m = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        signs.append(bracket(m + m.T, a))
    negative = sum(s < 0 for s in signs)
    assert 0 < negative < len(signs)
    out.update(negative_share=negative / len(signs))
    report(f"Q5 source not conserved: B is negative in {negative} of {len(signs)} random draws")

    report("RESULT: a conserved slow source radiates two transverse traceless polarizations through the pole, "
           "with the power (G/5) <Q''' Q'''> and Newton's G. The binary pulsar's decay tests it.")
    return out


if __name__ == "__main__":
    run()
