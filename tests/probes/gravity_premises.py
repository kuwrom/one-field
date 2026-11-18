"""The gravity lane's three readings, each read without G: the amount's share, the vector contact, the gates' thresholds.

THE QUESTION. The chain from the four integers to G reads three things
besides the integers: the bridge read as amount takes no share of b0
(horizon_amount.py A6), the vector rows enter the stiffness at the
covariant contact c = +1 (tests/audit/parameter_ledger.py), and gates
S and M are read at a threshold and in a band (g_path_trials.py G10).
What does each rest on without a measured G?

PREMISES. P1 the horizon: the amount past it carries no names
(interference/__init__.py, THE HORIZON). P2 "One covariant operator
gives a row both its share of b0 and its stiffness" (interference/
boundary.py, THE AMOUNT'S COUPLING). P3 the engine's transmutation:
1/alpha_G2 = 24 pi at M_b, S0 = 2 pi/(b0 alpha_G2), the lever arm to
mu* = M_b e^-(S0 - C2(26)), and the EM base 1/alpha_alg =
singh x 1/alpha_G2(mu*) (root.py). P4 the gauge field's kinetic
operator in Feynman gauge, the Hodge Laplacian d delta + delta d on
1-forms, with two scalar ghosts. P5 the measured 1/alpha(0) (Cs
recoil) as the test. No measured G is read.

STATEMENTS (the numbers below check each one).

  Q1 (the share). A row's share of b0 is its Dynkin index, a name
     (P2), and the amount carries none (P1), so its share is zero. The
     EM base depends on b0 directly, 1/alpha_alg = 48 b0/pi, so
     1/alpha(0) tests the share with no measured G.  One real 7 of the
     bridge in the running moves 1/alpha(0) by 1.6 percent, eighty
     million standard deviations.  The bridge read with its label
     (b0 = 19/3) puts 1/alpha(0) at 81.4.
  Q2 (the contact). The vector weight is 1/6 - c/2 with c the Ricci
     coefficient of the operator's Weitzenbock term. The Hodge
     Laplacian on 1-forms is the rough Laplacian plus Ricci
     (Bochner), so c = 1: on the unit four sphere, from the spectra of
     coexact and exact 1-forms less two ghosts, the heat trace's a1 is
     -8 = 2 x (1/6 - 1/2) x R, and the rough Laplacian gives +4 =
     2 x (1/6) x R. The same Weitzenbock term carries the field
     strength, the paramagnetic 4 of the vector's 11/3: without it b0
     is -16/3, S0 is negative and no scale transmutes below M_b.
  Q3 (the thresholds). On the one list every row is light or amount
     and vents with weight one, so the stiffness is 25.83 / 43.83 /
     61.83 at every regulator.  The sign needs no threshold. On the
     lattice of shares the boundary reading c is the whole
     (boundary_lattice.py), and it uses no gate. G depends on the ruler
     and the common mode, which root.py computes before gravity.py
     runs.  The threshold of S and the band of M decide no number for
     the sign or the stiffness. They state that the anchored regulator
     lies inside the boundary, t0 M_b^2 = 4.11 / 6.98 / 9.84. The
     verdicts on the boundary reading do depend on the threshold.  Three
     exclusions rest on the wired threshold: the GR factor 4 pi in k,
     which S passes under two of its six conventions (g_path_trials.py
     G10); c = 2 pi on the menu of shares (boundary_whole.py B3). And
     the 71 lattice values (boundary_lattice.py L3). The ground for
     k = 1 that uses no threshold is the coupling variable,
     coupling_normalization.py, an identification.

RESULT. The amount's zero share is the horizon reading on the operator that gives
a share, and 1/alpha(0) tests it without G. The covariant contact is
the Weitzenbock identity of the gauge field's own operator, the one
whose field strength term gives b0 its 11/3. No verdict on the sign
or the stiffness depends on the gates' thresholds. The exclusion of
4 pi from k does depend on the wired threshold.  coupling_normalization.py
grounds k = 1 without it, as an identification.

Scope: one loop heat kernel weights on the one list. The four sphere is
a check of the local coefficient, which holds on any background.

Usage: python3 tests/probes/gravity_premises.py
"""

import math
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

ALPHA_CS = (137.035999046, 27e-9)


def inv_alpha_of(b0):
    """1/alpha(0) from the engine's transmutation with the given b0 (P3), solved as root.py solves it."""
    inv_aG2 = 1 / r.alpha_G2_Pl
    S0 = 2 * math.pi / (b0 * r.alpha_G2_Pl)
    base = float(r.singh_ratio) * (inv_aG2 - b0 / (2 * math.pi) * (S0 - float(r.C2_26)))
    x = base
    for _ in range(200):
        x = base * (1 - 1 / (2 * math.pi) - 2 / (2 * math.pi * x)**2)
    return base, x, S0


def heat_a1(kappa, lmax=4000):
    """a1 of the gauge vector less two ghosts on the unit S^4, for Weitzenbock Ricci coefficient kappa."""
    vol, l = 8 * math.pi**2 / 3, np.arange(0, lmax, dtype=float)
    sc_l, sc_d = l * (l + 3), (2 * l + 3) * (l + 1) * (l + 2) / 6
    k = l[1:]
    co_l, co_d = (k + 1) * (k + 2), k * (k + 3) * (2 * k + 3) / 2      # coexact 1-forms, Hodge
    ex_l, ex_d = sc_l[1:], sc_d[1:]                                   # exact 1-forms
    shift = 3.0 * (1 - kappa)                                         # Ric = 3 on the unit S^4
    trace = lambda t: (np.sum(co_d * np.exp(-(co_l - shift) * t)) + np.sum(ex_d * np.exp(-(ex_l - shift) * t))
                       - 2 * np.sum(sc_d * np.exp(-sc_l * t)))
    ts = np.linspace(0.002, 0.02, 12)
    coeffs = np.polyfit(ts, [trace(t) * (4 * math.pi * t)**2 / vol for t in ts], 4)[::-1]
    return coeffs[0], coeffs[1]


def run(report=print):
    out = {}

    # Q1: the share.
    one = B.rows(r.d10, r.d11, r.n7, r.n26)
    b0 = float(B.b0(one))
    labelled = float(B.b0(B.rows(r.d10, r.d11, r.n7, r.n26, bridge="light")))
    seven = b0 - 1 / 6
    base, x, _ = inv_alpha_of(b0)
    assert abs(base - 48 * b0 / math.pi) < 1e-9 and abs(base - r._inv0) < 1e-9
    assert abs(x - r.inv_alpha_phys) < 1e-9
    shares = {}
    for name, b in (("amount, b0 = 32/3", b0), ("one real 7 running", seven), ("labelled, b0 = 19/3", labelled)):
        _, xi, _ = inv_alpha_of(b)
        shares[name] = (b, xi, (xi - ALPHA_CS[0]) / ALPHA_CS[1])
    assert abs(labelled - 19 / 3) < 1e-12
    assert abs(shares["one real 7 running"][1] / x - 1 + 0.0156) < 1e-3
    assert abs(shares["one real 7 running"][2]) > 5e5
    out.update(b0=b0, base_over_b0=base * math.pi / b0, shares=shares)
    report("Q1 1/alpha_alg = 48 b0/pi; 1/alpha(0) by reading: "
           + "; ".join(f"{k}: {v[1]:.6f} ({v[2]:+.2e} sigma)" for k, v in shares.items()))

    # Q2: the contact.
    hodge, rough = heat_a1(1.0), heat_a1(0.0)
    R = 12.0
    assert abs(hodge[0] - 2) < 1e-6 and abs(rough[0] - 2) < 1e-6
    assert abs(hodge[1] - 2 * (1 / 6 - 1 / 2) * R) < 1e-4 and abs(rough[1] - 2 * (1 / 6) * R) < 1e-4
    assert abs(B.a1("vector", 1.0, 0.0) - hodge[1] / (2 * R)) < 1e-5
    hv = r.hv_G2
    b0_no_spin = -Fraction(1, 3) * hv - Fraction(4, 3) * r.n_dirac_G2 * r.T_7
    assert b0_no_spin == Fraction(-16, 3)
    _, _, S_no_spin = inv_alpha_of(float(b0_no_spin))
    assert S_no_spin < 0
    out.update(a1_hodge=hodge[1], a1_rough=rough[1], b0_no_spin=float(b0_no_spin), S_no_spin=S_no_spin)
    report(f"Q2 S^4 heat trace, vector less two ghosts: a1 = {hodge[1]:+.4f} (Hodge, per state "
           f"{hodge[1] / (2 * R):+.4f}), {rough[1]:+.4f} (rough, {rough[1] / (2 * R):+.4f}). Without the "
           f"Weitzenbock term b0 = {float(b0_no_spin):.3f} and S0 = {S_no_spin:.2f}")

    # Q3: the thresholds.
    sig = {ct: [B.sigma(one, t, ct, float(r.xi_bridge)) for t in np.logspace(-3, 3, 13)] for ct in B.CONTACTS}
    assert all(max(v) - min(v) < 1e-12 and min(v) > 0 for v in sig.values())
    t0 = {ct: B.anchored(one, ct, float(r.xi_bridge))[0] / (2 * math.pi) for ct in B.CONTACTS}
    assert all(t > 1 for t in t0.values())
    out.update(sigma={ct: v[0] for ct, v in sig.items()}, t0=t0)
    report("Q3 stiffness at every regulator from 1e-3 to 1e3 M_b^-2: "
           + " / ".join(f"{v[0]:.2f}" for v in sig.values()) + ". Anchored t0 M_b^2 = "
           + " / ".join(f"{t:.2f}" for t in t0.values()))

    report("RESULT: the zero share is the horizon reading on the operator that gives a share, tested by 1/alpha(0) without G. The "
           "contact is the gauge operator's Weitzenbock term. No verdict on the sign or the stiffness depends on the thresholds. "
           "k = 1 against 4 pi depends on the wired threshold, and the coupling variable grounds it without the threshold "
           "(an identification).")
    return out


if __name__ == "__main__":
    run()
