"""
substrate3d.py: the one copy of the 3D substrate machinery.

Shared by the 3D dynamics probes knot_charge.py (charge
protection), confined_phase.py (Z3 disorder operator),
knot_spectrum.py (ring-down tower), and skyrmion_3d.py (formation),
so the grid, the Z3 Fourier projections, and the split-step
evolution exist exactly once.  The saturated-knot seed
(texture_state) is shared by the first three. skyrmion_3d builds
its own small-amplitude formation seed, and each charge-estimator
variant lives with its probe.  The equation is the repo's exact
Z3-coupled NLS (interference/nls_soliton.py), g1 = g0/sqrt(2).  An
absorbing sponge on the relative modes at r > 9 is the saturated-knot
probes' far-field treatment (all reported physics happens at r < 6).
The formation probe runs sponge-free and norm-conserving.
"""

import os
import sys

import numpy as np
import scipy.fft as sfft

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "interference"))
from nls_soliton import G1_OVER_G0  # noqa: E402

g0, rho0 = 1.0, 1.0
g1 = g0 * G1_OVER_G0
omega = np.exp(2j * np.pi / 3)
DT = 0.005


def grid(N, L=24.0):
    dx = L / N
    xs = (np.arange(N) - N // 2) * dx
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12
    k1 = 2 * np.pi * sfft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    return dx, X, Y, Z, R, KX**2 + KY**2 + KZ**2


def texture_state(N, L, amp, size):
    """Unit-winding texture (winding_texture profile) on the
    condensate."""
    dx, X, Y, Z, R, K2 = grid(N, L)
    f = np.pi * np.exp(-R / size)
    nx, ny, nz = X / R, Y / R, Z / R
    c, s = np.cos(f), np.sin(f)
    chi1, chi2 = c + 1j * s * nz, s * (1j * nx - ny)
    a = amp / np.cosh(R / size)
    psi = np.empty((3, N, N, N), dtype=complex)
    for k in range(3):
        psi[k] = np.sqrt(rho0) * (1.0 + omega**k * a * chi1
                                  + omega**(2 * k) * a * chi2)
    return psi, dx, R, K2


def rel_modes(psi):
    p1 = (psi[0] + omega**(-1) * psi[1] + omega**(-2) * psi[2]) / 3.0
    p2 = (psi[0] + omega**(-2) * psi[1] + omega**(-4) * psi[2]) / 3.0
    return p1, p2


def rel_mode(psi):
    return rel_modes(psi)[0]


def evolve(psi, K2, R, n_steps, g1_eff=g1, callback=None, every=0,
           sponge=True):
    """Split-step evolution. Optional relative-mode sponge at r > 9.

    sponge=True (saturated-knot probes) absorbs outgoing relative
    radiation at the far field. sponge=False (formation probe,
    matching the review's sim3d run) keeps the evolution exactly
    norm-conserving, with conservation trackable as a diagnostic.
    callback(psi), if given, is called every `every` steps (dense
    observables: node detection, ring-down time series).  The
    evolution itself is identical with or without a callback.
    """
    kin = np.exp(-0.5j * K2 * DT)
    lam = 0.5 * (1 + np.tanh((R - 9.0) / 1.5)) if sponge else None
    for n in range(1, n_steps + 1):
        psi = sfft.ifftn(sfft.fftn(psi, axes=(1, 2, 3)) * kin,
                         axes=(1, 2, 3))
        rho_t = (np.abs(psi)**2).sum(axis=0)
        for c in range(3):
            V = g0 * rho_t + g1_eff * np.real(np.conj(psi[(c + 1) % 3])
                                              * psi[(c + 2) % 3])
            psi[c] *= np.exp(-1j * V * DT)
        psi = sfft.ifftn(sfft.fftn(psi, axes=(1, 2, 3)) * kin,
                         axes=(1, 2, 3))
        if sponge:
            psi -= (lam * DT) * (psi - psi.mean(axis=0)[None])
        if every and callback is not None and n % every == 0:
            callback(psi)
    return psi
