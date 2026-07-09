"""
dynamics.py: the territory-side probe, Eq. 2.1 evolved directly.

First dynamical test of the narrative's two core mechanisms, run on
the framework's own equation (the Z3-symmetric coupled NLS with the
G2 constraint g1 = g0/sqrt(2)), 1D, split-step Fourier:

  TEST 1 (forgetting as erosion + the shadow):
    a localized v-biased relative excitation (the "newborn loop",
    v = (1, -1/2, -1/2)) sheds its bias by radiation while a
    persistent common-mode DEPLETION forms at the site.
    The bias/relative column is a NULL CONTROL: the initial state is
    purely v-patterned, so the ratio is trivially constant (= |v|
    norm), its constancy certifies the regime is linear (no mode
    mixing), and it does NOT yet test selectivity.  A genuine
    selectivity test needs a superposition of independent relative
    patterns (v and w = (0, 1, -1)) with their decay rates compared.
    That is the next experiment.

  TEST 2 (a knot and its shadow):
    a dark-bright vector structure stays localized and its
    total-density dip persists while the trial profile relaxes.

STATUS (honest): bias erosion and shadow formation are OBSERVED
(qualitative).  Strict stability of a gapped knot is NOT yet
demonstrated, the trial ansatz broadens. The open task is the true
stationary profile (imaginary-time solver) and a 2D run where
windings are topological.

Runtime ~20 s.  Usage: python3 dynamics.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'interference')))
from root import d10

G1_OVER_G0 = 1.0 / np.sqrt(d10)          # the G2 constraint


def make_grid(N=1024, L=160.0):
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    k = 2*np.pi*np.fft.fftfreq(N, d=L/N)
    return x, k, L/N


def evolve(psi, k, G, dt, steps):
    Kin = np.exp(-0.5j * k**2 * dt)
    for _ in range(steps):
        dens = np.abs(psi)**2
        psi = psi * np.exp(-1j * dt * (G @ dens))
        psi = np.fft.ifft(np.fft.fft(psi, axis=1) * Kin, axis=1)
    return psi


def run(report=print):
    x, k, dx = make_grid()
    g0 = 1.0
    g1 = g0 * G1_OVER_G0
    G = np.array([[g0, g1, g1], [g1, g0, g1], [g1, g1, g0]])
    dt, quarter = 0.002, 2000            # 4 quarters -> t = 16
    v = np.array([1.0, -0.5, -0.5])
    c = np.abs(x) < 6.0
    rho0 = 1.0 / 3.0

    # ── TEST 1: biased excitation, erosion, selectivity, shadow ────
    env = 0.4 * np.exp(-x**2 / 8.0)
    psi = np.sqrt(rho0) * (1.0 + v[:, None] * env[None, :]) + 0j
    report("TEST 1: v-biased excitation (newborn loop)")
    report(f"  {'t':>4} {'bias(center)':>13} {'bias/relative':>14} "
           f"{'shadow dn':>11}")
    out = []
    for i, t in enumerate([0, 4, 8, 16]):
        if i:
            psi = evolve(psi, k, G, dt, quarter if t <= 8 else 2*quarter)
        dens = np.abs(psi)**2
        n = dens.sum(0)
        d = dens - n / 3.0                        # relative content
        b = (v[:, None] * dens).sum(0)            # v-weighted bias
        bias_c = float(np.sqrt((b[c]**2).mean()))
        rel_c = float(np.sqrt((d[:, c]**2).sum(0).mean()))
        shadow = float((n[c] - n.mean()).mean())
        ratio = bias_c / rel_c if rel_c else float('nan')
        out.append((t, bias_c, ratio, shadow))
        report(f"  {t:4d} {bias_c:13.5f} {ratio:14.4f} {shadow:+11.6f}")
    eroded = out[-1][1] < 0.6 * out[0][1]
    shadowed = out[-1][3] < 0.0
    report(f"  bias erosion: {'OBSERVED' if eroded else 'NOT OBSERVED'}."
           f"  persistent shadow: {'OBSERVED' if shadowed else 'NOT OBSERVED'}")

    # ── TEST 2: dark-bright trial knot ──────────────────────────────
    w = 1.5
    psi = np.zeros((3, len(x)), dtype=complex)
    psi[0] = np.sqrt(rho0) * np.tanh(x / w)
    psi[1] = np.sqrt(0.5 * rho0) / np.cosh(x / w)
    psi[2] = np.sqrt(rho0) * np.ones(len(x))
    d0 = np.abs(psi)**2
    dip0 = float((d0.sum(0)[c] - d0.sum(0).mean()).mean())
    psi = evolve(psi, k, G, dt, 4 * quarter)
    d1 = np.abs(psi)**2
    dip1 = float((d1.sum(0)[c] - d1.sum(0).mean()).mean())
    localized = d1[1].max() > 3 * np.median(d1[1])
    report("\nTEST 2: dark-bright trial knot")
    report(f"  total-density dip: {dip0:+.5f} (t=0) -> {dip1:+.5f} (t=16)")
    report(f"  bright component localized at t=16: "
           f"{'YES' if localized else 'NO'} "
           f"(profile broadened: stationary state NOT yet demonstrated)")
    return {'erosion': eroded, 'shadow': shadowed,
            'dip_persist': abs(dip1) > 0.5 * abs(dip0),
            'localized': bool(localized)}


if __name__ == "__main__":
    run()
