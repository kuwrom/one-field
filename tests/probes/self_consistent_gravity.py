"""Self-consistent graviton correction and induced effective potential.

The composite metric g = h/U is made of the field. Its fluctuations
propagate on themselves. The Sakharov ledger at one loop counts every
matter row but not the metric's own loop. This probe computes the
self-consistent correction: the metric induces a stiffness, that
stiffness governs the metric's own fluctuations, and those fluctuations
contribute to the stiffness.

PREMISES. P1 the proper time cutoff rho = theta(t - t0) of FACT 1.
P2 the locality lemma (composite_pole.py K9): F0/F2 = -2 at every
order, so the metric's q^2 term is always the linearized EH operator.
P3 the matter ledger Sigma_matter from gravity.py (the rows that are
NOT the composite metric itself).

THE SELF-CONSISTENT EQUATION.  At one loop of the matter fields,
1/G_1 = Sigma_matter / (2 pi t0). The composite metric (spin 2) then
propagates with this stiffness. Its own one loop contribution to the
stiffness is the graviton row.

The graviton row per the heat kernel expansion (Bastianelli, Bonezzi,
Melis 2206.13287, Christensen-Duff 1979):
    on shell:      a1 = -8/3  (gauge invariant, per spin-2 field)
    de Donder:     a1 = -23/3  (off shell, includes gauge fixing)
    type Ia:       a1 = -11/3
    type Ib:       a1 = -35/8
In units where a real scalar = +1/6.

STATEMENTS.

  SC1 (the graviton row is additive). The a1 coefficient depends on the
      field's Lorentz structure (spin 2), not on the coupling. It enters
      Sigma as an additive constant, just like every other row: one
      graviton propagator ~ 1/Sigma, one vertex ~ Sigma, the factors
      cancel under the loop integral.

  SC2 (the equation is additive). Because a1 is universal (it depends
      on spin, not coupling), the corrected stiffness is simply
      Sigma_eff = Sigma_matter + c_grav. No iteration, no fixed point.
      Higher loops change the graviton row by O(c/Sigma) ~ 10%.

  SC3 (CW sextic). The metric loop on a background of density R induces
      a CW effective potential. The R^3 (sextic) term is suppressed by
      1/Sigma^3, far below the needed saturating floor.

RESULT. The on shell graviton row is -8/3, gauge invariant at leading
order in 1/N (subleading correction O(10%)). The literature range
[-7.67, -2.67] is gauge scheme variation, not physical uncertainty.
The graviton moves the covariant convention ratio from 0.698 to 0.775.
The substrate's p^4 dispersion does not enter: it belongs to a
condensate, which gate C keeps out of the vacuum.
Sigma != 6pi under the covariant convention: this is a structural
statement about the field content, not a correctable gap. The CW
sextic is killed (200x too small). The Sakharov ratio remains a
regulator dependent diagnostic (FACT 7). The framework's G prediction
is the anchor inversion (M_Pl at +22 ppm, root.py).

Usage: python3 tests/probes/self_consistent_gravity.py
"""

import math
import os
import sys


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root  # noqa: E402

C16 = 1.0 / (16 * math.pi**2)


def sigma_from_rows():
    """The matter ledger's Sigma (no graviton loop)."""
    d10, d11, n7, n26 = root.d10, root.d11, root.n7, root.n26
    N_bridge = n7 * n26
    a1 = {"scalar": 1/6, "weyl": 1/12, "gauge": -1/3}
    n_gauge = 2 * d11 * (d11 + 1)
    n_fermion = d10 * d11 * (d11 * d10**2 + d10 + 1)
    n_coset = d10 * (2 * n7 - (d11**2 - 1))
    n_G2_scal = n7
    sigma = (n_gauge * a1["gauge"]
             + n_fermion * a1["weyl"]
             + n_coset * a1["gauge"]
             + n_G2_scal * a1["scalar"]
             + N_bridge * a1["scalar"])               # the amount, minimally coupled
    return sigma


GRAVITON_ROW = {
    "on_shell": -8/3,
    "de_Donder": -23/3,
    "type_Ia": -11/3,
    "type_Ib": -35/8,
}


def higher_order_estimate(sigma_matter, c1):
    """Estimate the size of the 2-loop correction to the graviton row.

    The graviton row a1 is universal at leading order. At next order,
    the background curvature shifts it by O(c1/Sigma_matter): the
    graviton propagating in the induced background sees corrections
    proportional to G * R ~ 1/Sigma. Returns the fractional shift.
    """
    return abs(c1) / sigma_matter


def cw_sextic_coefficient(n_dof, sigma, t0):
    """CW induced sextic from the metric's one loop effective potential.

    The metric couples to matter density with strength C = 2pi/Sigma.
    The third derivative of the CW potential gives the sextic term
    kappa_CW = n pi^2 / (2 Sigma^3) * t0, suppressed by 1/Sigma^3.
    """
    return n_dof * math.pi**2 / (2 * sigma**3) * t0


def run(report=print):
    sigma_matter = sigma_from_rows()
    target = 6 * math.pi
    bare_ratio = target / sigma_matter

    report(f"\nSC0: matter ledger Sigma = {sigma_matter:.4f}, "
           f"target 6pi = {target:.6f}, "
           f"bare ratio = {bare_ratio:.6f}")

    # SC1: the graviton row, four schemes.
    report("\nSC1: graviton a1 coefficient (literature, scalar = 1/6 units):")
    for name, val in GRAVITON_ROW.items():
        sigma_eff = sigma_matter + val
        ratio = target / sigma_eff
        report(f"  {name:15s}: a1 = {val:+.4f}, "
               f"Sigma_eff = {sigma_eff:.4f}, "
               f"6pi/Sigma_eff = {ratio:.6f}")

    # SC2: relative importance (1/N_eff expansion parameter).
    report("\nSC2: relative importance (1/N_eff):")
    for name, val in GRAVITON_ROW.items():
        r = abs(val) / sigma_matter
        report(f"  {name:15s}: |c|/Sigma = {r:.4f} = {1/r:.1f}^{{-1}}")

    # SC3: the value needed for exact closure at 6pi.
    c_needed = target - sigma_matter
    lo, hi = min(GRAVITON_ROW.values()), max(GRAVITON_ROW.values())
    report(f"\nSC3: for 6pi/Sigma_eff = 1, need c = {c_needed:.4f}")
    report(f"  on_shell gives  {GRAVITON_ROW['on_shell']:+.4f} -> "
           f"ratio {target/(sigma_matter + GRAVITON_ROW['on_shell']):.6f}")
    report(f"  de_Donder gives {GRAVITON_ROW['de_Donder']:+.4f} -> "
           f"ratio {target/(sigma_matter + GRAVITON_ROW['de_Donder']):.6f}")
    report(f"  needed          {c_needed:+.4f}")
    report(f"  literature range: [{lo:.4f}, {hi:.4f}]")
    in_range = lo <= c_needed <= hi
    report(f"  the needed value is {'inside' if in_range else 'outside'} the "
           f"scheme spread. The spread is gauge artifact (SC6), so this is "
           f"not a closure route.")

    # SC4: the on shell (gauge invariant) answer.
    # At leading order in 1/N, the composite graviton is indistinguishable
    # from a fundamental one (Carone-Claringbold-Vaman 1710.09367).
    # The bilinear structure constrains the fluctuations just as gauge
    # invariance does: h is not an independent field.
    c_phys = GRAVITON_ROW["on_shell"]
    sigma_1 = sigma_matter + c_phys
    ratio_1 = target / sigma_1
    report(f"\nSC4: on shell (gauge invariant) graviton row:")
    report(f"  Sigma_eff = {sigma_1:.4f}, "
           f"6pi/Sigma_eff = {ratio_1:.6f}, "
           f"deficit = {(ratio_1 - 1)*100:+.2f}%")

    # SC5: higher order estimate.
    frac = higher_order_estimate(sigma_matter, c_phys)
    report(f"\nSC5: 2-loop fractional shift to the graviton row: "
           f"|c|/Sigma = {frac:.4f} = {1/frac:.0f}^{{-1}}")
    report(f"  the row is perturbatively stable at 10%.")

    # SC6: the on shell value IS the physical answer.
    # The literature range [-23/3, -8/3] is gauge scheme variation, not
    # physical uncertainty. At leading order in 1/N the composite
    # graviton matches the fundamental one (Carone et al. 1710.09367),
    # and the gauge invariant answer is the on shell value: -8/3.
    # The de Donder value includes gauge fixing artifacts that shift the
    # OPERATOR (not just the ghost subtraction): the gauge fixed h-sector
    # operator has a1 = -97/6, while the pure Lichnerowicz on 10
    # symmetric tensor components has a1 = +13/6 (POSITIVE, wrong sign).
    # The physical operator is the linearized Einstein tensor, which is
    # degenerate. Its nondegenerate projection onto TT modes gives -8/3.
    report(f"\nSC6: the on shell value (-8/3) is the physical answer.")
    report(f"  gauge fixed h-operator a1 = {-97/6:.4f} (NOT the Lichnerowicz)")
    report(f"  pure Lichnerowicz a1      = {+13/6:.4f} (positive, wrong sign)")
    report(f"  on shell (TT modes)       = {c_phys:+.4f} (gauge invariant)")
    report(f"  The literature range is scheme variation, not physics.")

    # SC7: the substrate's p^4 dispersion belongs to a condensate, not to the vacuum.
    report("\nSC7: no p^4 correction to the vacuum term.")
    report("  The Bogoliubov dispersion belongs to a condensate, and a condensate")
    report("  filling space needs xi > 2.7e6 l_Pl to keep c_T = c")
    report("  (carrier_sakharov.py), so it does not regulate the empty field.")

    # SC8: CW induced sextic from the metric's loop.
    n_grav_dof = 5
    t0_val = float(root.d11)
    kappa = cw_sextic_coefficient(n_grav_dof, sigma_matter, t0_val)
    report(f"\nSC8: CW sextic from the graviton loop:")
    report(f"  kappa_CW = {kappa:.2e}, floor = O(1)")
    report(f"  ratio = {kappa:.2e}, suppressed by {1/kappa:.0e}")
    report(f"  KILL: the metric loop does not generate the saturating term.")

    report("\n" + "=" * 60)
    report("SUMMARY")
    report("=" * 60)
    report(f"Sigma_matter = {sigma_matter:.4f}, target 6pi = {target:.4f}")
    report(f"bare ratio              = {bare_ratio:.4f}")
    report(f"+ on shell graviton     = {target/sigma_1:.4f} "
           f"(deficit {(ratio_1-1)*100:+.1f}%)")
    report(f"CW sextic               = KILLED ({1/kappa:.0e}x too small)")
    report(f"\nThe on shell graviton row is a 10% correction in the right")
    report(f"direction.  It shrinks the gap but does not close it.")
    report(f"The Sakharov ratio remains a diagnostic (FACT 7), not 1.")
    report(f"The framework's G prediction is the anchor (+22 ppm).")

    return dict(
        sigma_ledger=sigma_matter, target=target,
        bare_ratio=bare_ratio, literature=GRAVITON_ROW,
        c_needed=c_needed, physical_ratio=ratio_1,
        kappa_CW=kappa,
    )


if __name__ == "__main__":
    run()
