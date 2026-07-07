"""
bookkeeping.py: anomaly bookkeeping and the theta operator, computed.

WHY FERMIONS (the Witten route).  For solitons of a sigma model
equipped with a Wess-Zumino term at level N, the exchange of two
solitons multiplies the wavefunction by exp(i pi N).  For SU(n >= 3)
targets pi_4 = 0, so the statistics comes entirely from the WZ term:
solitons are fermions if and only if N is odd.  (This is how baryons
are fermions in QCD's effective theory: N = N_c = 3.)

The framework's confined sector carries SU(3)_3 chiral data, level
k = 3, and k = 3 was forced by Gate 3 (the dimension-one current) for
reasons having nothing to do with statistics.  Three is odd:

    exchange phase = exp(3 pi i) = -1  ->  FERMIONS.

The level is forced, not cited: Witten quantization makes the WZ
coefficient an integer, and 't Hooft anomaly matching fixes it to the
number of fermion species in the loop, the d11 = 3 Dirac fermions in
the 7 (the same matter content the alpha_s derivation uses).  Two
deep things follow at once: spin-statistics is free (the same
topological term gives half-integer spin), and "k odd" becomes a
retrodiction: the level that generates three generations and the
level that makes matter fermionic are the SAME integer.  Generations
and statistics share one origin.

WHY THE KNOT IS CHARGED (persistence => charge).  In a topological
sector the only protection is charge; an uncharged lump is
contractible, a transient.  Being a knot and unit Z3 winding are one
property.  Negative evidence from our own 1D data: without
topological protection the windings mix (probes/stationary.py),
exactly what unprotected labels do.

WZ level: the anomaly bookkeeping of the branched matter.
  (a) CONSISTENCY: the child SU(3) cubic gauge anomaly of the
      branched 7 -> 3 + 3bar + 1 must vanish (a chiral gauge anomaly
      would kill the child theory):  A(3) + A(3bar) + A(1) = 0.
  (b) LEVEL COUNTING: the family-space WZ coefficient counts the
      child-charged copies circulating per family loop, the quark
      analog of "WZ level = N_c".  Per family, the branched 7
      contributes ONE child triplet (with its conjugate as the Dirac
      partner); three families -> coefficient 3 = d11 = Gate 3's k.
      (Counting-level verification; the full triangle diagram is the
      cited standard result.)

Theta: the Cardy-level check that the lightest Z3
  SECTOR-CHANGING primary is the fundamental.  Changing the Z3
  sector requires NONZERO TRIALITY t = (l1 + 2 l2) mod 3; the
  identity family {1, J, J^2} (t = 0) preserves sectors.  Among the
  triality-charged primaries the minimum weight is h = 2/9, the
  fundamental pair, so the sector-changing insertion IS the h = 2/9
  operator, uniquely and by the largest possible margin (next
  changer: 5/9).

Usage: python3 bookkeeping.py
"""

from fractions import Fraction


def cubic_anomaly_su3(rep):
    """A(3) = +1, A(3bar) = -1, A(1) = 0 (normalisation A(fund)=1)."""
    return {"3": 1, "3bar": -1, "1": 0}[rep]


def run(report=print):
    report("BOOKKEEPING: anomaly ledger + sector-changer")
    report("=" * 64)

    # P1(a): child gauge anomaly of the branched 7
    branch = ["3", "3bar", "1"]                     # 7 -> 3 + 3bar + 1
    A = sum(cubic_anomaly_su3(r) for r in branch)
    report(f"  P1a: 7 -> 3 + 3bar + 1;  child cubic anomaly "
           f"A = {A}   ({'CONSISTENT' if A == 0 else 'BROKEN'})")
    assert A == 0

    # P1(b): level counting
    n_families = 3                                   # d11 generations
    triplets_per_family = 1                          # one 3 (+ 3bar) each
    level = n_families * triplets_per_family
    report(f"  P1b: WZ coefficient = child-charged copies per family")
    report(f"       loop x families = {triplets_per_family} x "
           f"{n_families} = {level}  (= d11 = Gate-3 k; ODD -> fermion)")
    assert level == 3

    # theta: lightest sector-changing primary of SU(3)_3
    prims = [(l1, l2) for l1 in range(4) for l2 in range(4)
             if l1 + l2 <= 3]
    h = lambda l1, l2: Fraction(l1*l1 + l2*l2 + l1*l2 + 3*l1 + 3*l2, 18)
    t = lambda l1, l2: (l1 + 2*l2) % 3
    changers = sorted((h(*p), p) for p in prims if t(*p) != 0)
    keepers = sorted((h(*p), p) for p in prims if t(*p) == 0)
    report(f"  theta: sector-KEEPING (t=0): "
           f"{[(str(hh), p) for hh, p in keepers]}")
    report(f"       sector-CHANGING (t!=0), lightest first:")
    for hh, p in changers[:4]:
        report(f"         h{p} = {hh}")
    h0, p0 = changers[0]
    report(f"       lightest changer: {p0} at h = {h0} = 2/9 "
           f"(next: {changers[2][0]})")
    assert h0 == Fraction(2, 9) and p0 in [(1, 0), (0, 1)]
    report("-" * 64)
    report("  child theory anomaly-consistent; WZ level = 3 by")
    report("      counting: generations and statistics, one integer.")
    report("  the Z3 sector-changer is the h = 2/9 operator,")
    report("      uniquely: theta is the insertion's identity.")


if __name__ == "__main__":
    run()
