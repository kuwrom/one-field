"""The scope of "coherence carries one flow": the initial state, not every time.

THE QUESTION. The dust reading (dark_sector.py, T10) reads the amount
as a count with no pressure: "kinetic pressure needs constituents
moving against the common frame, and coherence carries one flow, the
frame of the state it was left in". Halos are held up by velocity
dispersion. Is "one flow" a statement about the state the amount is
left in, or about every later time? Which reading does the count
support?

PREMISES. P1 the amount is a count conserved per comoving volume that
couples only through the common mode (dark_sector.py): no pressure, no
collisions, every element moving on its own. P2 the plane wave
growing mode of a cold count in one dimension, the Zel'dovich map
x(q, s) = q + s sin q, with s the growth of the perturbation. The map
is exact for one dimensional pressureless collisionless motion before
and after the first crossing. P3 the alternative reading, one
velocity at every place at every time: the pressureless flow that
stays single valued is the adhesion limit, where crossing elements
stick and the flow carries no dispersion.

STATEMENTS (the numbers below check each one).

  O1 (the state left in). For s < 1 every place holds one element of
     the count: one stream, one velocity, local dispersion at the cell
     width (9e-4). The count is one flow.
  O2 (the first fold). At s = 1 the map's Jacobian 1 + s cos q
     vanishes at q = pi and the single valued velocity gradient
     diverges. At that point, a finite time in, the count stops being a
     flow field.
  O3 (the count folds). At s = 3 the count holds three streams over
     0.51 of the period, carrying 0.87 of the mass, with a local
     velocity dispersion of 0.49 of the perturbation's velocity
     amplitude there (single stream places: the cell width, 2e-4), and
     no pressure: the dispersion is the phase mixing of one cold sheet.
  O4 (one flow at every time fails). Held single valued at every time,
     the same count must stick where it crosses (P3), and the folded
     region carries no dispersion: one velocity at a place has none,
     against O3's 0.49. A pressureless count then cannot supply a
     halo's dispersion. A field level coherent state could supply it
     only through quantum pressure, below a de Broglie length set by a
     mass the framework does not name (dust_gate_state.py S1).

RESULT. "Coherence carries one flow" is a statement about the state
the amount is left in: a cold sheet, one stream at each place. The
count's own dynamics folds that sheet, and the folded count carries
the dispersion halos show, with no pressure, as cold dark matter does.
Read as one flow at every time, it fails at the first fold. The
statement is a reading whose scope is the initial state.

Usage: python3 tests/probes/one_flow_scope.py
"""

import math

import numpy as np

N = 200_000
Q = np.linspace(0.0, 2 * math.pi, N, endpoint=False)


def streams(s, cells=2000):
    """Per place: the number of streams, the fraction of places and of mass with three, and the local dispersion there."""
    x = Q + s * np.sin(Q)                          # unwrapped Lagrangian map
    v = np.sin(Q)                                  # velocity per unit growth rate
    jac = 1 + s * np.cos(Q)
    count = np.zeros(cells, int)
    cuts = np.flatnonzero(np.diff(np.sign(jac)) != 0) + 1
    for seg in np.split(np.arange(N), cuts):       # each monotone branch of x(q) is one stream
        lo, hi = x[seg].min(), x[seg].max()
        cells_hit = np.arange(int(math.floor(lo / (2 * math.pi) * cells)), int(math.ceil(hi / (2 * math.pi) * cells)))
        np.add.at(count, np.mod(cells_hit, cells), 1)
    cell = (np.mod(x, 2 * math.pi) / (2 * math.pi) * cells).astype(int)
    multi = count >= 3
    in_multi = multi[cell]
    # local dispersion: the spread of velocities within each place, mass weighted over the folded places
    sums = np.bincount(cell, minlength=cells).astype(float)
    mean = np.bincount(cell, v, cells) / np.maximum(sums, 1)
    var = np.bincount(cell, v * v, cells) / np.maximum(sums, 1) - mean**2
    disp = float(math.sqrt(max(0.0, (var * sums)[multi].sum() / max(sums[multi].sum(), 1)))) if multi.any() else 0.0
    single = float(math.sqrt(max(0.0, (var * sums)[~multi].sum() / sums[~multi].sum())))
    return count, float(multi.mean()), float(in_multi.mean()), disp, single, jac


def run(report=print):
    out = {}

    # O1: before the first fold, one stream everywhere.
    count, frac, mass, disp, single, jac = streams(0.5)
    assert jac.min() > 0 and count.max() == 1 and single < 2e-3
    report(f"O1 s = 0.5: Jacobian above {jac.min():.2f}, one stream at every place, local dispersion {single:.0e} (the cell width)")

    # O2: the first fold.
    s1 = 1.0
    gradient = lambda s: np.max(np.abs(np.cos(Q) / (1 + s * np.cos(Q) + 1e-300)))   # dv/dx = cos q / J
    g_before, g_near = gradient(0.9), gradient(0.999)
    assert abs((1 + s1 * np.cos(Q)).min()) < 1e-9 and g_near > 50 * g_before
    report(f"O2 s = 1: the Jacobian vanishes at q = pi; |dv/dx| grows from {g_before:.0f} at s = 0.9 to {g_near:.0f} at 0.999")

    # O3: the count folds into three streams with dispersion and no pressure.
    count, frac, mass, disp, single, _ = streams(3.0)
    assert count.max() == 3 and abs(frac - 0.51) < 0.01 and abs(mass - 0.87) < 0.01
    assert abs(disp - 0.49) < 0.01 and single < 1e-3
    out.update(fold_fraction=frac, fold_mass=mass, dispersion=disp)
    report(f"O3 s = 3: three streams over {frac:.2f} of the period carrying {mass:.2f} of the mass, "
           f"velocity dispersion {disp:.2f}, pressure 0")

    # O4: held single valued, a place carries one velocity, so its local dispersion is zero by construction.
    out.update(single_stream_dispersion=single)
    report(f"O4 one velocity at every place: local dispersion 0 by construction, against the folded count's {disp:.2f}")
    return out


if __name__ == "__main__":
    run()
