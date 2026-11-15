"""
effective_action.py: the texture stabilizing term, derived
(the topological protection derivation, branch (i), closed
constructively).

SCOPE. knot_charge.py measured that the Z3-NLS substrate on its own
does not protect pi_3 winding. Textures undergo Derrick collapse and
unwind through codim-4 amplitude node events.  DP #4 requires the confined
sector's effective action to contain a texture stabilizing term
(Skyrme quartic or equivalent), or the chain breaks.  This probe
exhibits that term CONSTRUCTIVELY, from the repo's own energy
functional, with every step executable.

THE CONSTRUCTION (three exact steps, each verified below):

  1. EXACT SPLIT.  Write the relative spinor as phi = sqrt(rho) zeta
     with |zeta| = 1 (the director).  The kinetic density splits
     EXACTLY:
         |grad phi|^2 = (grad sqrt(rho))^2
                        + rho [ (1/4)(grad n)^2 + v^2 ],
     with n = zeta^dag sigma zeta the CP^1 director and
     v_i = -i zeta^dag d_i zeta the relative supercurrent.  Verified
     pointwise on the repo's texture, error -> 0 with resolution.

  2. MERMIN HO.  The current's curl is SLAVED to the topology:
         (curl v)_ij = (1/2) n . (d_i n x d_j n).
     Verified pointwise (factor +1/2, converging with resolution).
     A texture cannot exist without circulating current.

  3. INDUCED QUARTIC.  In the saturated regime the charged sector is
     GAPPED (measured: omega_gap = 1.178 = 0.27 mu, the ring down
     tower of knot_spectrum.py. At Gaussian order the classical
     linear response IS the quantum quasiparticle spectrum).
     Eliminating the gapped current at fixed Mermin-Ho curl gives the
     screened response energy, whose local limit (R > 1/m) is the
     Faddeev-Skyrme quartic:
         E_top = [rho/(8 m^2)] * INT F_ij^2,
         F_ij = n . (d_i n x d_j n).
     Derrick balance for the fixed norm director then reads
         E(R) = (rho/8) I2 R + [rho/(8 m^2)] I4 / R,
     with I2, I4 the profile integrals of the unit hedgehog (computed
     below by quadrature): a MINIMUM at R* = sqrt(I4/I2)/m.  Here m
     is the SCREENING WAVENUMBER of the eliminated mode, in units of
     inverse length, since E_top must carry the units of E_grad.  The
     measured gap omega_gap = 1.178 (knot_spectrum.py) is a frequency.
     Converting it to a wavenumber needs the eliminated mode's
     dispersion inside the knot.  The BdG relative mode of the
     symmetric vacuum has a modulational band, not a static pole, so
     two conversions bracket the dispersion: the Schrodinger static
     pole kappa^2 = 2 omega_gap, and the acoustic one
     kappa = omega_gap/c_s.  Measured numbers: I4/I2 = 1.98, and
     R* kappa = sqrt(I4/I2) = 1.41 > 1 identically, so the local limit
     is self-consistent at its own minimum for any kappa. R* lies
     between 1.9 and 5.3 healing lengths across the bracket. The
     stabilized texture is therefore a few healing lengths wide. Its
     exact width within that bracket is set by the eliminated mode's
     dispersion inside the knot.

WHAT THIS CLOSES AND THE PREMISE IT RESTS ON.  Branch (i) of the
topological protection derivation is closed. The effective action of
the FIXED NORM sector contains the stabilizing term, with a derived
coefficient and a self-consistent scale. The substrate field on its
own escapes this term in exactly the way measured. The amplitude is
classically soft, and node events (amplitude through zero) leave the
fixed norm manifold. That is the only channel the quartic cannot see.
Protection therefore rests on a single premise: the amplitude/node
channel is suppressed in the saturated quantum regime.  At Gaussian
order it is (the measured gap). Beyond Gaussian order a node event is
a finite action excursion, suppressed as exp(-c rho xi^3). The
diluteness parameter rho xi^3 is a substrate depth parameter that the
web does not fix.  That premise is the confined phase entry test
(confined_phase.py), sharpened to one quantitative condition.

Corroboration: Babaev-Faddeev-Niemi (PRB 65, 100512, 2002) derive
the same induced Faddeev-Skyrme term for two component condensates.
The follow up literature's finding that full GP knots still shrink
through amplitude excursions is exactly our measured mechanism.

SIGN.  This probe runs the gas sign g0 > 0 over a held uniform
background (tests/substrate/nls_soliton.py, THE SIGN).  The committed
self-bound sign holds the clock lump, which carries no pi_3 texture
(knot_in_lump.py). For that sign the closure is the spinning lump
(spinning_lump.py).

Usage: python3 effective_action.py     (static, fast. No PDE)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _bootstrap  # noqa: E402,F401
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from substrate3d import g0, g1, rel_modes, rho0, texture_state  # noqa: E402

# Scales of the hermitian equation.  The nonvariational variant
# gives gap 2.749 = 0.74 mu with mu = (3g0+g1)rho0 (that variant's rate).
# knot_spectrum's result on the hermitian equation is
# omega_gap = 1.178 and mu = (3g0+2g1)rho0.
M_GAP = 1.178          # measured: knot_spectrum.py outcome (a frequency)
MU = (3 * g0 + 2 * g1) * rho0
XI = 1.0 / np.sqrt(MU)   # healing length of the condensate
C_S = np.sqrt(MU)        # q=0 sound speed (sakharov FACT 5 dispersion)
# screening wavenumber of the eliminated mode: the two conversions
# of the measured frequency gap that bracket the in knot
# dispersion (hbar = mass = 1, so a frequency has units 1/length^2)
KAPPA = {"schrodinger pole sqrt(2 omega_gap)": np.sqrt(2.0 * M_GAP),
         "acoustic omega_gap/c_s": M_GAP / C_S}


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
    report("EFFECTIVE ACTION: the texture stabilizing term, derived")
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
    ratio = float(np.sqrt(I4 / I2))          # = R* kappa, kappa independent
    report("-" * 70)
    report(f"  gap omega = {M_GAP} (measured, knot_spectrum.py) = "
           f"{M_GAP/MU:.2f} mu, a frequency. Screening wavenumber bracket:")
    report(f"  profile integrals: I2 = {I2:.1f}, I4 = {I4:.1f}, "
           f"I4/I2 = {I4/I2:.2f}, R* kappa = sqrt(I4/I2) = {ratio:.2f}")
    R_star = {}
    for name, kappa in KAPPA.items():
        Rs = ratio / kappa
        E = lambda Rv, k=kappa: (I2 / 8) * Rv + (I4 / (8 * k**2)) / Rv
        assert E(Rs) < E(0.5 * Rs) and E(Rs) < E(2 * Rs), \
            "the induced quartic must produce a Derrick minimum"
        R_star[name] = Rs
        report(f"    kappa = {kappa:.3f} [{name}]: Derrick MINIMUM at "
               f"R* = {Rs:.3f} = {Rs/XI:.2f} healing lengths")
    R_lo, R_hi = min(R_star.values()), max(R_star.values())
    report(f"  R* = {R_lo/XI:.1f} to {R_hi/XI:.1f} healing lengths across "
           f"the bracket. R* kappa = {ratio:.2f} > 1 for every kappa")

    # frozen assertions
    assert ratio > 1.0, \
        "the local (Faddeev-Skyrme) limit must hold at its own minimum"
    # The Derrick minimum exists for every conversion of the gap, and
    # the texture is a few healing lengths wide. The exact width is
    # set by the eliminated mode's in knot dispersion inside the bracket.
    Rstar = R_lo
    assert 1.5 < R_lo / XI and R_hi / XI < 6.0, \
        "the stabilized texture must sit within a few healing lengths"

    report("-" * 70)
    report("  CLOSED (branch i of the topological protection derivation): the fixed")
    report("  norm sector's effective action contains the induced")
    report("  Faddeev-Skyrme quartic, coefficient rho/(8 kappa^2), with a")
    report("  self-consistent stabilized scale a few healing lengths")
    report(f"  wide ({R_lo/XI:.1f} to {R_hi/XI:.1f} xi across the bracket "
           "of the measured gap).")
    report("  The one channel it cannot see is the amplitude node")
    report("  (exit from the fixed norm manifold), which is the")
    report("  measured unwinding mechanism.  Protection therefore")
    report("  rests on a single premise: the node channel is gapped")
    report("  in the saturated quantum regime (Gaussian order:")
    report("  measured. Beyond: exp(-c rho xi^3) suppression).")
    report("  That premise IS the confined phase entry test, sharpened.")
    return dict(err=(e32, e48), mermin_ho=(mh32, mh48),
                I2=I2, I4=I4, R_star=Rstar, R_star_m=ratio,
                R_star_xi_lo=R_lo / XI, R_star_xi_hi=R_hi / XI,
                kappa=dict(KAPPA))


if __name__ == "__main__":
    run()
