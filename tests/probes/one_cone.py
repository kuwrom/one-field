"""One cone: does the knot carrying branch share the phonon's light cone?

Gate C of the gravity lane (interference/gravity.py) asks whether matter
and metric propagate on one cone. It closes at the carrier level in the
empty field (tests/probes/lump_in_empty_field.py). This probe computes
the internal speeds of a uniform clock medium, the ones a clock filling
space would carry. At Gaussian order about the uniform clock (the q=1
condensate of z3_vacuum.py, phi_a = sqrt(rho) w^a) the three branches
are closed form (tests/probes/gravity_tensor.py):

    omega_0     = sqrt(epsilon (epsilon + 2 mu)),          epsilon = p^2/2
    omega_+/-   = sqrt((epsilon + 5r/4)^2 - r^2/4) +/- 3r/4,
    r = g1 rho,  mu = (3 g0 - g1) rho,  hbar = bare mass = 1.

THE TEST. A relativistic particle at low momentum reads
    omega = Delta + p^2/(2 m*) - p^4/(8 m*^3 c^2) + ...
with the Einstein relation m* c^2 = Delta. A Galilean band
omega = Delta + p^2/2m* is that expansion only if m* c^2 = Delta, with
c the metric's speed. Were the phonon the metric, that speed would be
c_ph^2 = mu. For the relative branches, expanding in epsilon,

    Delta_+/- = sqrt(21) r/4 +/- 3r/4,
    1/m*      = (5r/4)/(sqrt(21) r/4) = 5/sqrt(21)   (both branches),
    p^4 term  = -(r^2/4)/(8 E0^3),  E0 = sqrt(21) r/4.

WHAT IS MEASURED (no parameter is chosen here; g1/g0 = 1/sqrt(2) is
the repo's standing value, nls_soliton.G1_OVER_G0):

  1. Einstein ratio m* c_ph^2 / Delta for each relative branch.
  2. The same quantity with the branch's own best speed: the shifted
     frequency Omega = omega -/+ 3r/4 obeys
         Omega^2 = E0^2 + (5r/4) p^2 + p^4/4
     exactly, so the relative pair is a charged Klein-Gordon doublet
     with speed c_rel^2 = 5r/4 seen at a charge chemical potential 3r/4,
     plus the Bogoliubov p^4/4 term. m* c_rel^2 = E0 exactly, but
     c_rel != c_ph, and the physical gaps are E0 -/+ 3r/4, not E0.
  3. A c-free test: the p^4 coefficient of a relativistic particle is
     fixed by Delta and m*. The branch misses it, because the gap lies
     at the healing scale where the Bogoliubov term is order one.
  4. Pair no go: m*_+ = m*_- while Delta_+ - Delta_- = 3r/2, so no
     single c satisfies the Einstein relation for both members of the
     doublet at any r > 0.
  5. No limiting speed on any branch, the phonon included: the group
     velocity grows as p at large p (Galilean free particle).

All analytic numbers are cross checked against eigenvalues of the
six component Nambu matrix of gravity_tensor.py, fitted at small p.
The symmetric background (1,1,1) is reported for completeness: it is
not the vacuum (z3_vacuum.py) and its relative band curves downward.

SCOPE. This is the uniform vacuum quasiparticle band. A bound knot's
own dispersion E(P) is not measured anywhere in the repository
(effective_action.py calls the in knot dispersion a registered
target). This probe measures only the branch the knot is built from.

Usage: python3 tests/probes/one_cone.py     (under a second)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from gravity_tensor import (bdg_matrix, clock_background,  # noqa: E402
                            clock_frequencies)
from nls_soliton import G1_OVER_G0  # noqa: E402

SQRT21 = math.sqrt(21.0)


def clock_cone_constants(mu, r):
    """Closed form low momentum constants of the clock vacuum."""
    e0 = SQRT21 * r / 4
    inv_m = (5 * r / 4) / e0
    return dict(
        c_ph=math.sqrt(mu),
        e0=e0,
        gap_lower=e0 - 3 * r / 4,
        gap_upper=e0 + 3 * r / 4,
        m_star=1.0 / inv_m,
        quartic=-(r * r / 4) / (8 * e0**3),
        c_rel=math.sqrt(5 * r / 4),
        charge_shift=3 * r / 4,
    )


def fitted_bands(A, B, momenta=(0.002, 0.004, 0.006, 0.008, 0.01)):
    """Quadratic fit in p^2 of the relative branches, phonon slope.

    Uses only eigenvalues of the Nambu matrix, not the analytic formula.
    """
    p = np.asarray(momenta)
    rows = []
    for value in p:
        eig = np.linalg.eigvals(bdg_matrix(A, B, value * np.array([0.0, 0.0, 1.0])))
        assert max(abs(eig.imag)) < 1e-10
        rows.append(np.sort(eig.real)[3:])
    rows = np.array(rows)
    phonon_slope = rows[0, 0] / p[0]
    lower = np.polyfit(p**2, rows[:, 1], 2)
    upper = np.polyfit(p**2, rows[:, 2], 2)
    return dict(c_ph=phonon_slope,
                gap_lower=lower[2], half_inv_m_lower=lower[1], quartic_lower=lower[0],
                gap_upper=upper[2], half_inv_m_upper=upper[1], quartic_upper=upper[0])


def group_velocity_ratio(k, mu, r, h=1e-5):
    """v_g(k)/k for the three branches; a limiting speed needs this -> 0."""
    lo, hi = clock_frequencies(k - h, mu, r), clock_frequencies(k + h, mu, r)
    return (hi - lo) / (2 * h) / k


def symmetric_background(g0=1.0, g1=G1_OVER_G0, rho=1.0):
    """Relative band on the (1,1,1) state: omega^2 = (eps + a1)^2 - b1^2."""
    a1, b1 = -2.5 * g1 * rho, -g1 * rho
    mu = (3 * g0 + 2 * g1) * rho
    gap = math.sqrt(a1 * a1 - b1 * b1)
    m_star = gap / a1                       # d omega/d eps = a1/gap < 0
    return dict(c_s=math.sqrt(mu), gap=gap, m_star=m_star,
                ratio=m_star * mu / gap)


def run(report=print):
    phi, mu, r, A, B = clock_background()
    k = clock_cone_constants(mu, r)
    fit = fitted_bands(A, B)

    # 0. The fitted matrix spectrum reproduces the closed forms.
    assert abs(fit["c_ph"] - k["c_ph"]) < 1e-5
    assert abs(fit["gap_lower"] - k["gap_lower"]) < 1e-10
    assert abs(fit["gap_upper"] - k["gap_upper"]) < 1e-10
    for half in (fit["half_inv_m_lower"], fit["half_inv_m_upper"]):
        assert abs(half - 1 / (2 * k["m_star"])) < 1e-7
    for quart in (fit["quartic_lower"], fit["quartic_upper"]):
        assert abs(quart - k["quartic"]) < 1e-5
    assert abs(k["m_star"] - SQRT21 / 5) < 1e-15
    report(f"Clock vacuum: r = g1 rho = {r:.12f}, mu = {mu:.12f}")
    report(f"  phonon  c_ph = sqrt(mu) = {k['c_ph']:.12f} "
           f"(matrix slope {fit['c_ph']:.9f})")
    report(f"  relative gaps Delta_-/+ = {k['gap_lower']:.12f} / "
           f"{k['gap_upper']:.12f}")
    report(f"  band mass m* = sqrt(21)/5 = {k['m_star']:.12f} on both "
           f"(matrix 1/2m* = {fit['half_inv_m_lower']:.9f})")

    # 1. Einstein ratio against the phonon cone.
    ratio_lower = k["m_star"] * k["c_ph"]**2 / k["gap_lower"]
    ratio_upper = k["m_star"] * k["c_ph"]**2 / k["gap_upper"]
    assert ratio_lower > 7 and ratio_upper > 1.5
    c_einstein_lower = math.sqrt(k["gap_lower"] / k["m_star"])
    c_einstein_upper = math.sqrt(k["gap_upper"] / k["m_star"])
    report(f"1. m* c_ph^2 / Delta = {ratio_lower:.9f} (lower), "
           f"{ratio_upper:.9f} (upper); Lorentz on the phonon cone needs 1")
    report(f"   speed the Einstein relation would need: "
           f"{c_einstein_lower:.9f} = {c_einstein_lower / k['c_ph']:.6f} c_ph "
           f"(lower), {c_einstein_upper:.9f} = "
           f"{c_einstein_upper / k['c_ph']:.6f} c_ph (upper)")

    # 2. Hidden Klein-Gordon form with its own speed and a charge shift.
    kg_error = 0.0
    for p in (0.0, 0.03, 0.4, 1.1, 3.0, 9.0):
        w = clock_frequencies(p, mu, r)
        for omega_shifted in (w[1] + k["charge_shift"], w[2] - k["charge_shift"]):
            target = k["e0"]**2 + k["c_rel"]**2 * p * p + p**4 / 4
            kg_error = max(kg_error, abs(omega_shifted**2 - target) / target)
    assert kg_error < 1e-12
    assert abs(k["m_star"] * k["c_rel"]**2 - k["e0"]) < 1e-14
    report(f"2. Omega = omega -/+ 3r/4 obeys Omega^2 = E0^2 + c_rel^2 p^2 + "
           f"p^4/4 (rel. error {kg_error:.1e})")
    report(f"   c_rel = sqrt(5r/4) = {k['c_rel']:.12f} = "
           f"{k['c_rel'] / k['c_ph']:.9f} c_ph; m* c_rel^2 = E0 = "
           f"{k['e0']:.12f} exactly. Charge shift 3r/4 = {k['charge_shift']:.12f}")

    # 3. c-free test at order p^4.
    quartic_rel_lower = -1 / (8 * k["m_star"]**2 * k["gap_lower"])
    quartic_rel_upper = -1 / (8 * k["m_star"]**2 * k["gap_upper"])
    q_lower = k["quartic"] / quartic_rel_lower
    q_upper = k["quartic"] / quartic_rel_upper
    q_kg = k["quartic"] / (-k["c_rel"]**4 / (8 * k["e0"]**3))
    assert abs(q_kg - 4 / 25) < 1e-12
    assert q_lower < 0.1
    report(f"3. p^4 coefficient {k['quartic']:.12f}. Relativistic value from "
           f"(Delta, m*): ratio {q_lower:.9f} (lower), {q_upper:.9f} (upper); "
           f"against the Klein-Gordon form (E0, c_rel): 4/25")

    # 4. The doublet cannot share one Einstein speed.
    split = k["gap_upper"] - k["gap_lower"]
    assert abs(split - 1.5 * r) < 1e-14
    assert abs(fit["half_inv_m_lower"] - fit["half_inv_m_upper"]) < 1e-9
    report(f"4. equal band masses, gaps split by 3r/2 = {split:.12f}: no single "
           "c gives m* c^2 = Delta for both branches at any r > 0")

    # 5. No limiting speed anywhere.
    vg = {kk: group_velocity_ratio(kk, mu, r) for kk in (10.0, 40.0)}
    assert all(abs(v - 1).max() < 2e-3 for v in vg.values())
    report(f"5. v_g/p at p = 10: {vg[10.0]}. At p = 40: {vg[40.0]}. "
           "Every branch, phonon included, ends as a free particle. The "
           "phonon cone is a p << 1/xi statement only")

    sym = symmetric_background()
    assert sym["m_star"] < 0
    report(f"Symmetric (1,1,1) state (not the vacuum): gap {sym['gap']:.9f}, "
           f"m* = {sym['m_star']:.9f} < 0, m* c_s^2/Delta = {sym['ratio']:.9f}")
    report("RESULT: measured against the phonon cone, the relative branch of a "
           "uniform clock medium is Galilean. The medium has two internal "
           "speeds, not a single one. Gate C closes at the carrier level in the empty field, where "
           "neither is the cone (tests/probes/lump_in_empty_field.py).")
    return dict(mu=mu, r=r, constants=k, fit=fit,
                ratio_lower=ratio_lower, ratio_upper=ratio_upper,
                c_einstein_lower=c_einstein_lower,
                c_einstein_upper=c_einstein_upper,
                kg_error=kg_error, quartic_ratio_lower=q_lower,
                quartic_ratio_upper=q_upper, quartic_ratio_kg=q_kg,
                gap_split=split, group_velocity_ratio=vg, symmetric=sym)


if __name__ == "__main__":
    run()
