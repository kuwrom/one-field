"""Scale ties: which of the engine's scales the 125 GeV state moves.

THE QUESTION. The engine fixes every light scale by closure: an action
exponent below M_Pl or a ratio to another scale. A closure fixes a
value. It does not say which field the value follows. By the
low energy theorem, every mass that follows the 125 GeV field h
couples h to the gauge fields whose running that mass changes, so a
scale tied to h through a field relation is seen in Higgs data even
when it is not a particle (regulator_free_ratios.py R7). Which ties
are relations between values, and which are relations between fields?

PREMISES. The engine's scale ties as committed (root.py, masses.py,
couplings.py). The engine's alpha_s chain (couplings.py), evaluated at
125 GeV. The low energy theorem d(1/alpha_s)/d ln v as the gluon
coupling of h, normalised to the top's. The theorem is read at beta
order 1, the order of the threshold matching. The two loop chain moves
the top's slope by its b1 step but carries no two loop decoupling term,
so it mixes orders. It is kept as a cross check and gives the same
verdicts.

THE TIES.
  T1  M_Pl from m_e (anchor inversion, root.py): a closure. The
      stiffness follows no light field (log_class.py L2: 1.2e-33).
  T2  v_EW = M_Pl exp(-S_quark) and T3 mu* = M_Pl exp(-(S_lepton -
      C2(26))): closure exponents; v_EW/mu* is a fixed ratio.
  T4  lepton and quark masses: closure values (masses.py), but a chiral
      mass term needs the electroweak condensate by gauge invariance,
      so each follows h through its Yukawa operator: kappa_F = 1.
  T5  M_W, M_Z = g v/2: gauge, so they follow h: kappa_V = 1.
  T6  alpha_s(mu*) = pi/32 (couplings.py): a closure at mu*.
  T7  M_V = g v_break sqrt(C2(7)/6), v_break^2 = E[v^2] v_EW^2
      (couplings.py): a closure value. The field it follows is the
      exposure stated in couplings.py, THE GLUON CONDITION.

STATEMENTS (the numbers below check each one).

  S1 (reading A). h is the radial mode of the electroweak condensate
     and G2 is broken by a field of its own whose vacuum value the
     connector fixes: h moves T4 and T5 only, and kappa_g = 1 by
     construction. The cost is a colour singlet G2-face radial mode
     with a gluon coupling 21/4 times the top's per unit ln v_break.
  S2 (reading B). h also moves M_V (T7): kappa_g = -17/4 (R7).
  S3 (reading C). h is a dilaton of the transmuted scale, so it
     rescales every scale proportional to it: mu*, M_V and m_t together.
     With alpha_s pinned at mu*, kappa_g = +23/2.
     A dilaton is the fluctuation of the scale, and a relation
     between scales follows it by definition. THE GLUON CONDITION
     (the connector transports the vacuum value, not the field)
     therefore cannot hold for a dilaton: under the engine's ties,
     the measured gluon coupling admits reading A alone.

RESULT. Every tie in the engine is a closure between values. The
measured couplings decide which field each value follows. T4 and T5
follow h by gauge invariance, and T7 must not. That fixes h as the
radial mode of the electroweak condensate, not a dilaton of the
master scale. The reading carries a G2-face radial mode of its own
(couplings.py, tests/probes/g2_radial_mode.py).

Usage: python3 tests/probes/scale_ties.py     (about half a minute
standalone: the engine runs once)
"""

import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import couplings as cpl  # noqa: E402
from interference import masses, root  # noqa: E402

KAPPA_G = (0.92, 0.08)            # CMS, Nature 607, 60 (2022)


def _run_order1(a0, mu0, mu1, nf):
    """QCD transport at beta order 1: 1/alpha moves by b0 ln(mu1/mu0)."""
    return 1.0 / (1.0 / a0 + cpl._b0(nf) * math.log(mu1 / mu0))


def chain(order=2):
    """1/alpha_s(125 GeV) as a function of which scales follow v.

    order 2 is the engine's alpha_s chain, order 1 the low energy
    theorem's order."""
    with contextlib.redirect_stdout(io.StringIO()):
        m = masses.derive(root.derive())
    M_Pl = root.M_Pl_GeV
    mu0 = M_Pl * math.exp(-(root.S_lepton - float(root.C2_26)))
    a_star = math.pi / 32
    g_G2 = math.sqrt(4 * math.pi * a_star)
    v0, mt0 = m["v_EW_pred_GeV"], m["m_t"] / 1e3

    def inv_alpha(v, MV_follows, mt_follows, mu_follows, mu=125.0):
        r = v / v0
        mu_star = mu0 * (r if mu_follows else 1.0)
        MV = g_G2 * (v if MV_follows else v0) * math.sqrt(float(root.E_v2) / 3)
        mt = mt0 * (r if mt_follows else 1.0)
        lam3 = 1.0 - root.d11 * root.n7 * math.log(MV / mu_star)
        a = 1.0 / (1.0 / a_star - lam3 / (12 * math.pi))
        run = cpl._run_SM_2loop if order == 2 else _run_order1
        return 1.0 / run(run(a, mu_star, mt, 6), mt, mu, 5)

    def slope(*follows, h=1e-4):
        return (inv_alpha(v0 * math.exp(h), *follows)
                - inv_alpha(v0 * math.exp(-h), *follows)) / (2 * h)
    return slope, inv_alpha, v0


def readings(order):
    """Top slope, coset slope and kappa_g of readings A, B, C at one beta order."""
    slope, inv_alpha, v0 = chain(order)
    top = slope(False, True, False)
    kappa = {k: slope(*f) / top for k, f in
             (("A", (False, True, False)), ("B", (True, True, False)), ("C", (True, True, True)))}
    return top, slope(True, False, False), kappa, inv_alpha, v0


def run(report=print):
    top, coset_only, kappa, _, _ = readings(1)
    _, _, kappa_chain, inv_alpha, v0 = readings(2)
    lo, hi = KAPPA_G[0] - 2 * KAPPA_G[1], KAPPA_G[0] + 2 * KAPPA_G[1]
    passes = {k: lo <= v <= hi for k, v in kappa.items()}
    alpha_MZ = 1.0 / inv_alpha(v0, True, True, False, mu=91.1876)
    assert abs(alpha_MZ - 0.11839) < 5e-5
    assert passes == {"A": True, "B": False, "C": False}
    assert abs(coset_only / top + 21 / 4) < 1e-6
    assert abs(kappa["B"] + 17 / 4) < 1e-6 and abs(kappa["C"] - 23 / 2) < 1e-6
    assert {k: lo <= v <= hi for k, v in kappa_chain.items()} == passes
    report(f"chain alpha_s(M_Z) = {alpha_MZ:.5f}. d(1/alpha_s(125))/d ln v at beta order 1: "
           f"top {top:+.4f}, coset threshold alone {coset_only:+.4f} ({coset_only / top:+.2f} x top).")
    for k, name in (("A", "radial mode, G2 broken by its own field"),
                    ("B", "radial mode that also moves M_V"),
                    ("C", "dilaton of the transmuted scale (mu*, M_V, m_t)")):
        report(f"  reading {k} ({name}): kappa_g = {kappa[k]:+.2f}, "
               f"{'inside' if passes[k] else 'outside'} [{lo:.2f}, {hi:.2f}]")
    report(f"  two loop chain cross check: B {kappa_chain['B']:+.2f}, C {kappa_chain['C']:+.2f}, "
           "same verdicts.")
    report("RESULT: the measured gluon coupling admits reading A alone: h is the radial mode "
           "of the electroweak condensate, not a dilaton, and the G2 face carries a radial "
           "mode of its own.")
    return dict(kappa=kappa, passes=passes, top=top, coset_only=coset_only,
                kappa_chain=kappa_chain, alpha_s_MZ=alpha_MZ)


if __name__ == "__main__":
    run()
