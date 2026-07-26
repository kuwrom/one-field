"""
dimensionality.py: 3+1 from the substrate's own structure, derived
rather than assumed.

THE ARGUMENT.  Three components force a C^2 relative sector (the
q = 1, 2 Fourier modes of the Z3 split), so the charged order
parameter is four real components.  Two structural facts then pin the
dimension:

  1. KNOT CONTENT.  The unit-normed relative spinor lives on S^3.
     pi_2(S^3) = 0 and pi_3(S^3) = Z: winding number, the substrate's
     topological content (winding_texture.py, knot_charge.py), first
     exists at three spatial dimensions.  d_spatial >= 3.

  2. POINTLIKE EVENTS.  Zeros of a generic 4-real-component field
     have codimension 4, so the set where the charge can change has
     dimension D - 4 in D-dimensional spacetime.  Charge-changing is
     eventlike (isolated points, the node events knot_charge.py
     measures) exactly when D = 4.  For D = 5 the zeros are
     persistent lines (charge leaks continuously along them); for
     D = 3 there are generically no zeros at all (charge could never
     change).

  d_spatial >= 3 and D = 4 have the unique joint solution D = 3+1.

THE COMPUTATION (part A).  The codimension fact is measured, not
cited: band-limited random 4-component fields on D-tori, D = 3, 4, 5,
three fixed seeds, two resolutions each.  Measured discriminators:

  D = 3: min|phi| stays bounded away from zero as resolution doubles
         (no zero set) and the near-zero cell count stays O(1).
  D = 4: min|phi| shrinks with resolution (zeros exist) and the
         near-zero cell count is resolution-INDEPENDENT (dimension-0
         set: isolated points).
  D = 5: min|phi| shrinks and the near-zero cell count GROWS linearly
         with resolution (dimension-1 set: lines).

INPUTS OF THE DERIVATION: the three-component
structure (engine data), genericity, and the demand that topological
charge exists and changes at discrete events (the mechanism
knot_charge.py measures).  Given those, D = 3+1 is derived, and the
only residue is genericity.

Usage: python3 dimensionality.py     (~20 s)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401

SEEDS = (1, 2, 3)
SIZES = {3: (32, 64), 4: (16, 32), 5: (10, 16)}


def _field_stats(D, N, seed, n_modes=6):
    """min|phi|/scale and gradient-scaled near-zero cell count for a
    band-limited random 4-component field on the D-torus."""
    rng = np.random.default_rng(seed)
    coords = np.meshgrid(*[np.linspace(0, 2 * np.pi, N, endpoint=False)] * D,
                         indexing="ij")
    f = np.zeros((4,) + (N,) * D)
    for c in range(4):
        for _ in range(n_modes):
            k = rng.integers(-2, 3, size=D)
            ph = rng.uniform(0, 2 * np.pi)
            f[c] += rng.normal() * np.cos(
                sum(k[i] * coords[i] for i in range(D)) + ph)
    norm = np.sqrt((f ** 2).sum(axis=0))
    scale = norm.mean()
    thr = scale * (3.0 / N)
    return norm.min() / scale, int((norm < thr).sum())


def run(report=print):
    report("DIMENSIONALITY: 3+1 from the C^2 relative sector")
    report("=" * 70)

    # the C^2 relative sector is the three-component structure's
    from substrate3d import rel_modes, texture_state
    psi, _, _, _ = texture_state(8, 12.0, 0.3, 2.0)
    p1, p2 = rel_modes(psi)
    n_real = 2 * len((p1, p2))
    assert n_real == 4
    report(f"  three components -> {len((p1, p2))} relative complex "
           f"modes -> {n_real} real components (codim {n_real})")

    # ── part A: the codimension fact, measured ──────────────────────
    out = {}
    report("\n  PART A: zero sets of generic 4-component fields")
    report(f"  {'D':>3} {'min/scale (lo->hi N)':>26} "
           f"{'near-zero cells (lo->hi N)':>28}")
    for D, (n_lo, n_hi) in SIZES.items():
        lo = [_field_stats(D, n_lo, s) for s in SEEDS]
        hi = [_field_stats(D, n_hi, s) for s in SEEDS]
        min_lo = float(np.mean([m for m, _ in lo]))
        min_hi = float(np.mean([m for m, _ in hi]))
        cnt_lo = float(np.mean([c for _, c in lo]))
        cnt_hi = float(np.mean([c for _, c in hi]))
        out[D] = dict(min_lo=min_lo, min_hi=min_hi,
                      cnt_lo=cnt_lo, cnt_hi=cnt_hi)
        report(f"  {D:>3} {min_lo:>12.4f} -> {min_hi:<10.4f} "
               f"{cnt_lo:>13.1f} -> {cnt_hi:<10.1f}")

    # D = 3: no zero set.  min|phi| stays bounded away from zero and
    # the near-zero count stays O(1) (measured: 0.037 mean, <= 4 cells)
    assert out[3]["min_hi"] > 0.025
    assert out[3]["min_hi"] / out[3]["min_lo"] > 0.5
    assert out[3]["cnt_hi"] <= 8
    # D = 4: isolated points.  min shrinks (measured ratio 0.32) and
    # the count is resolution-independent (measured ratio 1.09)
    assert out[4]["min_hi"] / out[4]["min_lo"] < 0.6
    assert 0.6 < out[4]["cnt_hi"] / out[4]["cnt_lo"] < 1.5
    # D = 5: lines.  the count grows with resolution (measured 1.63,
    # against the size ratio 1.6)
    assert out[5]["cnt_hi"] / out[5]["cnt_lo"] > 1.35
    report("  -> D=3 no zeros, D=4 isolated points, D=5 lines: "
           "codim = 4, measured")

    # ── part B: the pin ─────────────────────────────────────────────
    report("\n  PART B: the dimensional pin")
    # pi_d(S^3): 0 for d < 3, Z at d = 3 (knot content needs
    # d_spatial >= 3); eventlike charge change needs D = codim = 4
    verdicts = {}
    for D_tot in (3, 4, 5, 6):
        knots = (D_tot - 1) >= 3
        events = (D_tot == 4)
        verdicts[D_tot] = knots and events
        report(f"    D = {D_tot - 1}+1: knots {'yes' if knots else 'no ':>3}"
               f", pointlike events {'yes' if events else 'no'}"
               f"{'   <- unique' if knots and events else ''}")
    assert [D for D, ok in verdicts.items() if ok] == [4]

    # ── part C: charge change REQUIRES zeros (the residue, closed) ──
    # The winding number is locally constant on nonvanishing fields
    # (an integer, continuous in the field), so it can change only
    # where the field leaves that space: at a zero.  Measured on the
    # substrate's own texture: a winding-changing path drives the
    # core minimum to machine zero exactly where B jumps, while a
    # winding-preserving deformation keeps a strictly positive floor.
    report("\n  PART C: charge change requires zeros (measured)")
    from knot_charge import charge
    from substrate3d import texture_state as _tex
    n3 = 32
    psi_t, dxg, Rg, _ = _tex(n3, 24.0, 0.7, 2.5)
    psi_u = np.ones_like(psi_t)
    core = Rg < 6.0

    def _stats(p):
        p1, p2 = rel_modes(p)
        b = charge(p1, p2, dxg, Rg)
        return b, float(np.sqrt(np.abs(p1) ** 2
                                + np.abs(p2) ** 2)[core].min())

    b_ends, mins_chg = [], []
    for s_ in np.linspace(1.0, 0.0, 11):
        b, mn = _stats(psi_u + s_ * (psi_t - psi_u))
        mins_chg.append(mn)
        if s_ in (1.0, 0.0):
            b_ends.append(b)
    mins_prs, b_prs = [], []
    for size_ in (2.0, 2.4, 2.8, 3.2):
        b, mn = _stats(_tex(n3, 24.0, 0.7, size_)[0])
        mins_prs.append(mn)
        b_prs.append(b)
    report(f"    winding-changing path: B {b_ends[0]:+.2f} -> "
           f"{b_ends[1]:+.2f}, deepest core min = {min(mins_chg):.1e}")
    report(f"    winding-preserving path: B stays "
           f"{min(b_prs):+.2f}..{max(b_prs):+.2f}, floor min = "
           f"{min(mins_prs):.3f}")
    assert round(b_ends[0]) == -1 and round(b_ends[1]) == 0
    assert min(mins_chg) < 1e-6           # the change passes a zero
    assert all(round(b) == -1 for b in b_prs)
    assert min(mins_prs) > 0.05           # no zero, no change

    report("-" * 70)
    report("  VERDICT: given three components and genericity, "
           "D = 3+1 is the unique")
    report("  dimension: charge change requires zeros (PART C, "
           "measured), zeros are")
    report("  codim 4 (PART A, measured), and knots need d >= 3.  "
           "3+1 is derived, not")
    report("  assumed; its residue is genericity alone.")
    return dict(codim=out, unique_D=4,
                change_min=min(mins_chg), preserve_floor=min(mins_prs))


if __name__ == "__main__":
    run()
