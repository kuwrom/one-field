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
Z3-coupled NLS (tests/substrate/nls_soliton.py), g1 = g0/sqrt(2).  An
absorbing sponge on the relative modes at r > 9 is the saturated-knot
probes' far-field treatment (all reported physics happens at r < 6).
The formation probe runs sponge-free and norm-conserving.
"""

import os
import sys

import numpy as np
import scipy.fft as sfft

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
from nls_soliton import (G1_OVER_G0,  # noqa: E402
                         interaction_matrix as _im_1d,
                         _apply_unitary)

g0, rho0 = 1.0, 1.0
g1 = g0 * G1_OVER_G0
omega = np.exp(2j * np.pi / 3)
DT = 0.005


def _interaction_matrix(psi, g0_, g1_, alpha=0.0):
    """The hermitian 3x3 interaction matrix on the 3D grid.

    Thin wrapper over nls_soliton.interaction_matrix, which is written
    for a (3, ...) field of any spatial rank, so the 1D and 3D probes
    integrate literally the same equation.
    """
    return _im_1d(psi, g0_, g1_, alpha)


def energy(psi, dx, g0_=g0, g1_=g1, alpha=0.0):
    """Total energy on the 3D grid (same functional as nls_soliton).

    GUARD: assumes a CUBIC grid with uniform spacing dx (K2 is
    rebuilt from psi.shape[1] alone); a non-cubic field would be
    silently mis-integrated, so it is rejected loudly instead."""
    assert psi.shape[1] == psi.shape[2] == psi.shape[3], (
        f"energy() requires a cubic grid, got {psi.shape[1:]}")
    ph = sfft.fftn(psi, axes=(1, 2, 3))
    N = psi.shape[1] * psi.shape[2] * psi.shape[3]
    # K2 is rebuilt by the caller's grid(); recompute here from shape
    n = psi.shape[1]
    k1 = 2 * np.pi * sfft.fftfreq(n, d=dx)
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    K2 = KX**2 + KY**2 + KZ**2
    kin = 0.5 * float(np.sum(K2 * np.abs(ph)**2)) / N
    rho = np.abs(psi)**2
    rho_t = rho.sum(axis=0)
    inter = 0.5 * g0_ * float(np.sum(rho_t**2))
    phase = np.exp(1j * alpha)
    for m in range(3):
        m1, m2 = (m + 1) % 3, (m + 2) % 3
        inter += g1_ * float(np.sum(rho[m] * np.real(
            phase * np.conj(psi[m1]) * psi[m2])))
    return (kin + inter) * dx**3


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


def isolated_seed(N, L, amp, size):
    """ONE localized relative-sector bump, spherically symmetric.

    texture_state seeds a WINDING texture; this seeds a nodeless
    spherically symmetric relative-sector structure instead, so that
    angular momentum stays a good quantum number of the relaxed
    background and a multiplet's degeneracy can be counted.
    """
    dx, X, Y, Z, R, K2 = grid(N, L)
    envelope = amp * np.exp(-(R / size) ** 2)
    psi = np.empty((3, N, N, N), dtype=complex)
    for k in range(3):
        psi[k] = np.sqrt(rho0) * (1.0 + omega**k * envelope)
    return psi, dx, R, K2


def h_apply(psi, K2, g0_=g0, g1_=g1):
    """H psi for the stationarity residual: kinetic + interaction."""
    lap = sfft.ifftn(-K2 * sfft.fftn(psi, axes=(1, 2, 3)), axes=(1, 2, 3))
    M = _interaction_matrix(psi, g0_, g1_)
    return -0.5 * lap + np.einsum('...ij,j...->i...', M, psi)


def relax(psi, K2, dx, n_iter, dtau, n0=None, g0_=g0, g1_=g1):
    """3D imaginary-time relaxation with norm projection.

    The 1D probes have carried this for a while (knot_circulant._relax);
    the 3D module had `evolve` but no relaxation, so no certified
    stationary 3D state could be built.  Same scheme: a half-step
    kinetic propagator in Fourier space, the exact exp(-dtau M)
    interaction substep by batched eigendecomposition, then rescale to
    the target norm.
    """
    if n0 is None:
        n0 = float(np.sum(np.abs(psi) ** 2)) * dx**3
    kin = np.exp(-0.5 * dtau * K2)
    for _ in range(n_iter):
        psi = sfft.ifftn(kin * sfft.fftn(psi, axes=(1, 2, 3)),
                         axes=(1, 2, 3))
        M = _interaction_matrix(psi, g0_, g1_)
        w, V = np.linalg.eigh(M)
        comp = np.moveaxis(psi, 0, -1)
        amp = np.einsum('...ji,...j->...i', V.conj(), comp) * np.exp(-dtau * w)
        psi = np.moveaxis(np.einsum('...ij,...j->...i', V, amp), -1, 0)
        psi *= np.sqrt(n0 / (float(np.sum(np.abs(psi) ** 2)) * dx**3))
    return psi


def stationarity(psi, K2, dx, g0_=g0, g1_=g1):
    """(mu, relative residual) of H psi = mu psi."""
    hp = h_apply(psi, K2, g0_, g1_)
    n0 = float(np.sum(np.abs(psi) ** 2)) * dx**3
    mu = float(np.real(np.sum(np.conj(psi) * hp)) * dx**3) / n0
    res = float(np.linalg.norm(hp - mu * psi) / np.linalg.norm(psi))
    return mu, res


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

    TWO STATED LIMITS OF THIS INTEGRATOR:
      * the sponge is NON-UNITARY and NON-VARIATIONAL, applied
        outside the Strang splitting: any sponge=True run conserves
        neither norm nor energy by construction.  Its use is
        confined by the module docstring (reported physics at
        r < 6, absorption at r > 9), nothing in the code enforces
        that radius discipline, so a probe reporting far-field
        physics under sponge=True would be wrong silently.
      * DT is a module-level constant (0.005) rather than an evolve()
        parameter, and no convergence study lives in this file; the
        1D module carries the splitting-order check
        (nls_soliton.derive), the 3D probes inherit DT on trust.
    """
    # HALF-step kinetic propagator for -1/2 grad^2 (the full-step
    # factor exp(-0.5j*K2*DT) applied twice per step would run at
    # m = 1/2).
    kin = np.exp(-0.25j * K2 * DT)
    lam = 0.5 * (1 + np.tanh((R - 9.0) / 1.5)) if sponge else None
    for n in range(1, n_steps + 1):
        psi = sfft.ifftn(sfft.fftn(psi, axes=(1, 2, 3)) * kin,
                         axes=(1, 2, 3))
        # Hermitian nonlinear substep: exp(-i M DT) with the exchange
        # terms (nls_soliton.interaction_matrix).  A diagonal-V
        # substep integrates one half of Hamilton's equations and
        # conserves no energy.
        M = _interaction_matrix(psi, g0, g1_eff)
        psi_pred = _apply_unitary(M, psi, DT)
        M = 0.5 * (M + _interaction_matrix(psi_pred, g0, g1_eff))
        psi = _apply_unitary(M, psi, DT)
        psi = sfft.ifftn(sfft.fftn(psi, axes=(1, 2, 3)) * kin,
                         axes=(1, 2, 3))
        if sponge:
            psi -= (lam * DT) * (psi - psi.mean(axis=0)[None])
        if every and callback is not None and n % every == 0:
            callback(psi)
    return psi
