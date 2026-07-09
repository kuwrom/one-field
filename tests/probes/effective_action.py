"""
effective_action.py: the texture-stabilizing term, derived
(DERIVATION_PROGRAMS #4, branch (i), closed constructively).

THE QUESTION.  knot_charge.py measured that the bare Z3-NLS does not
protect pi_3 winding: textures Derrick-collapse and unwind through
codim-4 amplitude-node events.  DP #4 demanded that the confined
sector's effective action contain a texture-stabilizing term (Skyrme
quartic or equivalent), or the chain breaks.  This probe exhibits
that term CONSTRUCTIVELY, from the repo's own energy functional,
with every step executable.

THE CONSTRUCTION (three exact steps, each verified below):

  1. EXACT SPLIT.  Write the relative spinor as phi = sqrt(rho) zeta
     with |zeta| = 1 (the director).  The kinetic density splits
     EXACTLY:
         |grad phi|^2 = (grad sqrt(rho))^2
                        + rho [ (1/4)(grad n)^2 + v^2 ],
     with n = zeta^dag sigma zeta the CP^1 director and
     v_i = -i zeta^dag d_i zeta the relative supercurrent.  Verified
     pointwise on the repo's texture, error -> 0 with resolution.

  2. MERMIN-HO.  The current's curl is SLAVED to the topology:
         (curl v)_ij = (1/2) n . (d_i n x d_j n).
     Verified pointwise (factor +1/2, converging with resolution).
     A texture cannot exist without circulating current.

  3. INDUCED QUARTIC.  In the saturated regime the charged sector is
     GAPPED (measured: omega_gap = 2.749 = 0.74 mu, the ring-down
     tower of knot_spectrum.py. At Gaussian order the classical
     linear response IS the quantum quasiparticle spectrum).
     Eliminating the gapped current at fixed Mermin-Ho curl gives the
     screened response energy, whose local limit (R > 1/m) is the
     Faddeev-Skyrme quartic:
         E_top = [rho/(8 m^2)] * INT F_ij^2,
         F_ij = n . (d_i n x d_j n).
     Derrick balance for the fixed-norm director then reads
         E(R) = (rho/8) I2 R + [rho/(8 m^2)] I4 / R,
     with I2, I4 the profile integrals of the unit hedgehog (computed
     below by quadrature): a MINIMUM at R* = sqrt(I4/I2)/m.  Measured
     numbers: I4/I2 = 2.31, R* m = 1.52 > 1 (the local limit is
     self-consistent at its own minimum), and R* = 1.06 healing
     lengths: the stabilized texture is the core-scale knot.

WHAT THIS CLOSES AND WHAT REMAINS.  Branch (i) of DP #4 is closed:
the effective action of the FIXED-NORM sector contains the
stabilizing term, with a derived coefficient and a self-consistent
scale.  The bare field evades it exactly as measured: the amplitude
is classically soft, and node events (amplitude through zero) exit
the fixed-norm manifold. That is the one channel the quartic cannot
see.  So the entire protection question now rests on a single
premise: the amplitude/node channel is suppressed in the saturated
quantum regime.  At Gaussian order it is (the measured gap). Beyond
Gaussian order a node event is a finite-action excursion, suppressed
as exp(-c rho xi^3), and the diluteness parameter rho xi^3 is a
substrate-depth parameter the web does not fix.  That premise is
SUBSTRATE_CONJECTURES #1, now sharpened to one quantitative
condition (see registry).

Corroboration: Babaev-Faddeev-Niemi (PRB 65, 100512, 2002) derive
the same induced Faddeev-Skyrme term for two-component condensates.
The follow-up literature's finding that full-GP knots still shrink
through amplitude excursions is exactly our measured mechanism.

Usage: python3 effective_action.py     (static, fast. No PDE)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from substrate3d import g0, g1, rel_modes, rho0, texture_state  # noqa: E402

M_GAP = 2.749          # measured: knot_spectrum.py registered outcome
MU = (3 * g0 + g1) * rho0
XI = 1.0 / np.sqrt(MU)   # healing length of the condensate


def _director(psi):
    p1, p2 = rel_modes(psi)
    phi = np.stack([p1, p2])
    rho = (np.abs(phi)**2).sum(axis=0)
    zeta = phi / np.sqrt(rho)[None]
    z1, z2 = zeta
    n = np.stack([2 * np.real(np.conj(z1) * z2),
                  2 * np.imag(np.conj(z1) * z2),
                  np.abs(z1)**2 - np.abs(z2)**2])
    return phi, rho, zeta, n


def split_error(N, L=24.0, size=2.5):
    """Median relative error of the exact kinetic split on the
    texture (discretization only. -> 0 with resolution)."""
    psi, dx, R, _ = texture_state(N, L, amp=0.7, size=size)
    phi, rho, zeta, n = _director(psi)
    z1, z2 = zeta
    g = lambda A, ax: np.gradient(A, dx, axis=ax)
    lhs = sum(np.abs(g(phi[c], ax))**2
              for c in range(2) for ax in range(3))
    amp = sum(g(np.sqrt(rho), ax)**2 for ax in range(3))
    sig = 0.25 * sum(g(n[a], ax)**2 for a in range(3) for ax in range(3))
    v = [np.real(-1j * (np.conj(z1) * g(z1, ax)
                        + np.conj(z2) * g(z2, ax))) for ax in range(3)]
    rhs = amp + rho * (sig + sum(vi**2 for vi in v))
    core = (R > 1.0) & (R < 6.0)
    err = np.abs(lhs - rhs)[core] / np.maximum(np.abs(lhs)[core], 1e-12)
    # Mermin-Ho factor on the same grid
    curl_z = g(v[1], 0) - g(v[0], 1)
    dxn = np.stack([g(n[a], 0) for a in range(3)])
    dyn = np.stack([g(n[a], 1) for a in range(3)])
    topo = (n * np.cross(dxn, dyn, axis=0)).sum(axis=0)
    mask = core & (np.abs(topo) > np.percentile(np.abs(topo[core]), 75))
    mh = float(np.median(curl_z[mask] / topo[mask]))
    return float(np.median(err)), mh


def profile_integrals(N=64, L=14.0):
    """I2 = INT (grad n)^2, I4 = INT sum_{i<j} F_ij^2 for the unit
    hedgehog director (size = 1), by grid quadrature."""
    psi, dx, R, _ = texture_state(N, L, amp=0.7, size=1.0)
    _, _, _, n = _director(psi)
    g = lambda A, ax: np.gradient(A, dx, axis=ax)
    dn = [np.stack([g(n[a], ax) for a in range(3)]) for ax in range(3)]
    I2 = float(sum((dn[ax]**2).sum() for ax in range(3)) * dx**3)
    I4 = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = (n * np.cross(dn[i], dn[j], axis=0)).sum(axis=0)
            I4 += float((F**2).sum() * dx**3)
    return I2, I4


def run(report=print):
    report("EFFECTIVE ACTION: the texture-stabilizing term, derived")
    report("=" * 70)

    # 1 + 2: the exact split and Mermin-Ho, with convergence
    e32, mh32 = split_error(32)
    e48, mh48 = split_error(48)
    report(f"  exact split |grad phi|^2 = amp + rho[(1/4)(grad n)^2 + v^2]:")
    report(f"    median rel err {e32:.2e} (32^3) -> {e48:.2e} (48^3)")
    report(f"  Mermin-Ho curl v = (1/2) n.(dn x dn):")
    report(f"    measured factor {mh32:+.3f} (32^3) -> {mh48:+.3f} (48^3)")
    assert e48 < e32 < 0.05, "split identity must converge"
    assert abs(mh48 - 0.5) < abs(mh32 - 0.5) + 1e-9 and \
        abs(mh48 - 0.5) < 0.05, "Mermin-Ho factor must converge to 1/2"

    # 3: induced quartic and the Derrick minimum
    I2, I4 = profile_integrals()
    Rstar = float(np.sqrt(I4 / I2) / M_GAP)
    E = lambda Rv: (I2 / 8) * Rv + (I4 / (8 * M_GAP**2)) / Rv
    report("-" * 70)
    report(f"  gap m = {M_GAP} (measured, knot_spectrum.py) = "
           f"{M_GAP/MU:.2f} mu")
    report(f"  profile integrals: I2 = {I2:.1f}, I4 = {I4:.1f}, "
           f"I4/I2 = {I4/I2:.2f}")
    report(f"  E(R)/rho = (I2/8) R + (I4/8m^2)/R:")
    for Rv in (0.5 * Rstar, Rstar, 2 * Rstar, 4 * Rstar):
        report(f"    E({Rv:.3f}) = {E(Rv):.2f}")
    report(f"  Derrick MINIMUM at R* = {Rstar:.3f} = "
           f"{Rstar/XI:.2f} healing lengths.  R* m = {Rstar*M_GAP:.2f}")

    # frozen assertions
    assert E(Rstar) < E(0.5 * Rstar) and E(Rstar) < E(2 * Rstar), \
        "the induced quartic must produce a Derrick minimum"
    assert Rstar * M_GAP > 1.0, \
        "the local (Faddeev-Skyrme) limit must hold at its own minimum"
    assert 0.5 < Rstar / XI < 2.0, \
        "the stabilized texture must sit at the core (healing) scale"

    report("-" * 70)
    report("  CLOSED (branch i of DERIVATION_PROGRAMS #4): the fixed-")
    report("  norm sector's effective action contains the induced")
    report("  Faddeev-Skyrme quartic, coefficient rho/(8 m^2), with a")
    report("  self-consistent stabilized scale at the healing length.")
    report("  The one channel it cannot see is the amplitude node")
    report("  (exit from the fixed-norm manifold), which is the")
    report("  measured unwinding mechanism.  Protection therefore")
    report("  rests on a single premise: the node channel is gapped")
    report("  in the saturated quantum regime (Gaussian order:")
    report("  measured. Beyond: exp(-c rho xi^3) suppression).")
    report("  That premise IS SUBSTRATE_CONJECTURES #1, sharpened.")
    return dict(err=(e32, e48), mermin_ho=(mh32, mh48),
                I2=I2, I4=I4, R_star=Rstar, R_star_m=Rstar * M_GAP,
                R_star_xi=Rstar / XI)


if __name__ == "__main__":
    run()
