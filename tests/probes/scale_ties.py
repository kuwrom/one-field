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
couplings.py). The engine's two loop alpha_s chain (couplings.py),
evaluated at 125 GeV. The low energy theorem
d(1/alpha_s)/d ln v as the gluon coupling of h, normalised to the top's.

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
     with a gluon coupling 4.5 times the top's per unit ln v_break.
  S2 (reading B). h also moves M_V (T7): kappa_g = -3.52 (R7).
  S3 (reading C). h is a dilaton of the transmuted scale, so it
     rescales every scale proportional to it: mu*, M_V and m_t together.
     With alpha_s pinned at mu*, kappa_g = +10.3 on the two loop chain.
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
(couplings.py).

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


def chain():
    """1/alpha_s(125 GeV) as a function of which scales follow v."""
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
        return 1.0 / cpl._run_SM_2loop(cpl._run_SM_2loop(a, mu_star, mt, 6), mt, mu, 5)

    def slope(*follows, h=1e-4):
        return (inv_alpha(v0 * math.exp(h), *follows)
                - inv_alpha(v0 * math.exp(-h), *follows)) / (2 * h)
    return slope, inv_alpha, v0


def run(report=print):
    slope, inv_alpha, v0 = chain()
    top = slope(False, True, False)
    readings = {"A": top, "B": slope(True, True, False), "C": slope(True, True, True)}
    kappa = {k: s / top for k, s in readings.items()}
    coset_only = slope(True, False, False)
    lo, hi = KAPPA_G[0] - 2 * KAPPA_G[1], KAPPA_G[0] + 2 * KAPPA_G[1]
    passes = {k: lo <= v <= hi for k, v in kappa.items()}
    alpha_MZ = 1.0 / inv_alpha(v0, True, True, False, mu=91.1876)
    assert abs(alpha_MZ - 0.11839) < 5e-5
    assert passes == {"A": True, "B": False, "C": False}
    assert abs(kappa["B"] + 3.52) < 0.02 and abs(kappa["C"] - 10.31) < 0.05
    report(f"chain alpha_s(M_Z) = {alpha_MZ:.5f}; d(1/alpha_s(125))/d ln v: top {top:+.4f}, "
           f"coset threshold alone {coset_only:+.4f} ({coset_only / top:+.2f} x top).")
    for k, name in (("A", "radial mode, G2 broken by its own field"),
                    ("B", "radial mode that also moves M_V"),
                    ("C", "dilaton of the transmuted scale (mu*, M_V, m_t)")):
        report(f"  reading {k} ({name}): kappa_g = {kappa[k]:+.2f}, "
               f"{'inside' if passes[k] else 'outside'} [{lo:.2f}, {hi:.2f}]")
    report("RESULT: the measured gluon coupling admits reading A alone: h is the radial mode "
           "of the electroweak condensate, not a dilaton, and the G2 face carries a radial "
           "mode of its own.")
    return dict(kappa=kappa, passes=passes, top=top, coset_only=coset_only,
                alpha_s_MZ=alpha_MZ)


if __name__ == "__main__":
    run()
