"""A knot inside the self-bound clock lump.

In the gravity lane, the closure's empty field realization is a
self-trapped lump with the clock inside it, on the self-bound sign
(lump_in_empty_field.py). In the matter probes, the knot is a unit pi_3
winding of the relative spinor. It is set up over the uniform symmetric
background on the gas sign (knot_charge.py, winding_texture.py). There
the relative instability keeps the relative sector populated. This probe
asks whether the two pictures describe one object: does a knotted state
exist inside the clock lump?

PREMISES. The nonrelativistic carrier on the self-bound sign,
g0 = -1, g1 = 1/sqrt 2, kappa = 1, in three dimensions. The 3D clock
lump at mu = mu0/2 (N = 121.34). The clock occupies the Z3 Fourier mode
q = 1. The relative spinor is the pair of the other two modes
(psi_hat_0, psi_hat_2), and its winding is read by the estimator of
knot_charge.py.

STATEMENTS (the numbers below check each one).

  K1 (no vacuum topology). With g1 != 0 the symmetry of the potential is
     the common U(1) alone (arena_rulers.py), so the lump's internal
     vacuum manifold is a circle times the two clock orientations
     q = 1, 2, and pi_3 of it is trivial. A knot needs a nonzero
     relative spinor over its core. The clock lump does not carry one.
  K2 (the seeded knot drains). Write a unit hedgehog texture into the
     relative spinor, at amplitude 0.5 of the lump's (winding estimate
     -0.84 at 48^3, -0.90 at 64^3), and at amplitude 0.9 with width 2.5
     (-0.91). At fixed N each flows to the plain clock lump. The
     relative population falls below 1e-3 of its start, the winding
     estimate goes to zero, and the energy ends at the clock lump's
     value, -7.8873, starting from +4.4 and +28.3. The winding is lost
     while the relative population is still finite, so the texture
     unwinds before it drains. No knotted stationary state exists
     inside the lump.
  K3 (what the Faddeev-Skyrme term cannot do). The induced quartic of
     effective_action.py stiffens the director of a relative spinor of
     fixed norm. Here the downhill direction is the norm itself. At
     g1 != 0 the populations of the Z3 modes do not conserve the norm,
     so no term in the director can block this direction.

RESULT. On the self-bound sign the clock lump holds no knot: a
texture written into it flows to the plain clock lump. The knot based
statements of the gas sign probes read a populated relative sector
over a uniform background, which the particle of the gravity lane does
not carry. Its persistence is its circulation (spinning_lump.py).

Usage: python3 tests/probes/knot_in_lump.py     (about a minute and a half)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from knot_charge import charge  # noqa: E402
import lump_formation as fm  # noqa: E402
import lump_in_empty_field as lf  # noqa: E402

G0, G1, KAPPA = -1.0, lf.G1, 1.0
W = np.exp(2j * np.pi / 3)
AXES = (1, 2, 3)


def seeded(n, amp, size, L=12.0, frac=0.5):
    """The 3D clock lump with a unit hedgehog in its relative spinor."""
    d = lf.lump(frac * lf.MU0)
    dx = 2 * L / n
    xs = (np.arange(n) - n // 2) * dx
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12
    F = d["profile"](R.ravel()).reshape(R.shape)
    f = np.pi * np.exp(-R / size)
    zeta = (np.cos(f) + 1j * np.sin(f) * Z / R, 1j * np.sin(f) * (X + 1j * Y) / R)
    hat = {1: F * math.sqrt(1 - amp**2), 0: amp * F * zeta[0], 2: amp * F * zeta[1]}
    psi = np.stack([sum(W**(q * a) * hat[q] for q in (0, 1, 2)) for a in range(3)]) / math.sqrt(3)
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    return psi, dx, R, KX**2 + KY**2 + KZ**2, d


def fourier(psi):
    return [sum(np.conj(W)**(q * a) * psi[a] for a in range(3)) / 3 for q in range(3)]


def energy(psi, dx, K2):
    rho = np.abs(psi)**2
    Rt = rho.sum(0)
    kin = 0.5 * np.sum(np.abs(np.fft.fftn(psi, axes=AXES))**2 * K2) / psi[0].size * dx**3
    ex = sum(rho[j] * np.real(np.conj(psi[(j + 1) % 3]) * psi[(j + 2) % 3]) for j in range(3))
    return float(kin + np.sum(0.5 * G0 * Rt**2 + KAPPA * Rt**3 / 3 + G1 * ex) * dx**3)


def flow(psi, dx, R, K2, steps=1600, dtau=0.05, every=400):
    """Semi implicit gradient flow at fixed total charge."""
    N = np.sum(np.abs(psi)**2) * dx**3
    rows = []
    for i in range(steps + 1):
        if i % every == 0:
            h = fourier(psi)
            rel = 3 * (np.sum(np.abs(h[0])**2) + np.sum(np.abs(h[2])**2)) * dx**3
            rows.append(dict(tau=i * dtau, relative=float(rel),
                             winding=charge(h[0], h[2], dx, R, r_core=6.0),
                             E=energy(psi, dx, K2)))
        if i == steps:
            break
        nl = fm.force(psi, G0, G1, KAPPA)
        psi = np.fft.ifftn((np.fft.fftn(psi, axes=AXES) - dtau * np.fft.fftn(nl, axes=AXES))
                           / (1 + 0.5 * dtau * K2), axes=AXES)
        psi *= math.sqrt(N / (np.sum(np.abs(psi)**2) * dx**3))
    return rows


def run(report=print, cases=((48, 0.5, 1.8), (48, 0.9, 2.5))):
    out = {}
    for n, amp, size in cases:
        psi, dx, R, K2, d = seeded(n, amp, size)
        rows = flow(psi, dx, R, K2)
        out[(n, amp, size)] = dict(rows=rows, E_clock=d["E_nr"])
        first, last = rows[0], rows[-1]
        report(f"K2 grid {n}^3, amplitude {amp}, width {size}: winding {first['winding']:+.3f} -> "
               f"{last['winding']:+.3f} (zero by tau = "
               f"{next(r['tau'] for r in rows if abs(r['winding']) < 0.05):.0f}). Relative charge "
               f"{first['relative']:.2f} -> {last['relative']:.4f}; E {first['E']:+.3f} -> "
               f"{last['E']:+.4f} against the clock lump {d['E_nr']:+.4f}.")
    report("RESULT: the seeded knot unwinds and its relative population drains into the clock. "
           "The self-bound clock lump holds no knotted stationary state.")
    return out


if __name__ == "__main__":
    run()
