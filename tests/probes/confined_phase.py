"""
confined_phase.py: the layer-bridge test (substrate-conjecture ledger
entry #1. Derivation program D6 of the gaps audit).

THE QUESTION.  The web's statistics/chirality chain routes through a
CONFINED TOPOLOGICAL PHASE of the substrate (unbroken Z3 centre): the
WZ citation, winding quantization, and the substrate-web bridge all
presuppose it, and knot_charge.py showed the BARE substrate protects
nothing: protection must live in the confined phase if it lives
anywhere.  This probe defines an executable order parameter and
measures it in the saturated-knot regime.

THE ORDER PARAMETER (Z3 disorder operator, exact lattice form).  The
triality-charged field is the relative mode p1 (Z3 Fourier q = 1).
Its phase chi = arg(p1) carries the charged direction, and the
disorder operator W(C) = <exp(i * winding of chi around C)> obeys an
AREA law exactly when vortex lines of chi proliferate: in the dilute
vortex gas, ln W(a) ~ -c * rho_v * a^2, with rho_v the vortex-line
density.  The executable order parameter is therefore rho_v itself,
measured EXACTLY: the fraction of elementary plaquettes (all three
orientations) whose wrapped phase circulation is nonzero, the
standard lattice vortex detector. No sampling, every plaquette is
counted.  rho_v > 0 in the saturated regime -> area law -> the
disorder signature confinement needs. rho_v = 0 -> perimeter law
(W = 1 for every contractible loop) -> the classical configuration
is topologically TRIVIAL in the charged sector.

EPOCHS.  t = 0 (imprinted knot: the unit texture forces a charged
vortex ring, the detector's positive control, verified by an exact
winding count around the ring). t = 2 (mid-collapse, node-event
epoch). t = 8 (the saturated B = 0 endpoint). And a g1 = 0 control
at t = 8.

WHAT THIS PROBE IS AND IS NOT.  It is the entry test made executable:
a defined order parameter, measured exactly, with a positive control
and a coupling control, frozen in CI.  It is NOT a demonstration of
topological order (gap + degeneracy need the quantum spectrum on
distinct backgrounds, the remaining obligation of
SUBSTRATE_CONJECTURES #1).  The measured answer recorded below: the
classical endpoint is DEFECT-FREE (rho_v = 0: perimeter law, no
confinement signature at the classical level).  Whatever realizes the
confined phase is therefore not a property of the classical substrate
configuration, consistent with knot_charge.py, and it sharpens where
the obligation lives: the quantum theory of the confined sector.

Usage: python3 confined_phase.py     (~2 min, mark slow)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# the shared 3D machinery (one copy: substrate3d.py)
from substrate3d import DT, evolve, g1, rel_mode, texture_state  # noqa: E402


def _wrap(x):
    return np.angle(np.exp(1j * x))


def vortex_density(chi, R, r_max=8.0):
    """Exact lattice vortex detector: fraction of elementary
    plaquettes (three orientations, centres within r < r_max) whose
    wrapped phase circulation is nonzero (|sum| > pi)."""
    hits, total = 0, 0
    for ax1, ax2 in ((0, 1), (0, 2), (1, 2)):
        d1 = _wrap(np.diff(chi, axis=ax1))          # link ax1
        d2 = _wrap(np.diff(chi, axis=ax2))          # link ax2
        sl1 = [slice(None)] * 3
        sl1[ax2] = slice(0, -1)
        sl2 = [slice(None)] * 3
        sl2[ax1] = slice(1, None)
        sl3 = [slice(None)] * 3
        sl3[ax2] = slice(1, None)
        sl4 = [slice(None)] * 3
        sl4[ax1] = slice(0, -1)
        circ = (d1[tuple(sl1)] + d2[tuple(sl2)]
                - d1[tuple(sl3)] - d2[tuple(sl4)])
        # plaquette-centre distance from the origin
        slc = [slice(None)] * 3
        slc[ax1] = slice(0, -1)
        slc[ax2] = slice(0, -1)
        Rc = R[tuple(slc)]
        mask = Rc < r_max
        hits += int((np.abs(circ[mask]) > np.pi).sum())
        total += int(mask.sum())
    return hits / total, hits


def ring_winding(chi, N):
    """Exact winding of chi around a small loop threading the
    imprinted vortex ring (positive control for the detector)."""
    c = N // 2
    a = 4
    i0, k0 = c + 3 - a // 2, c - a // 2   # centred on the ring piercing
    patch = chi[i0:i0 + a + 1, c, k0:k0 + a + 1]
    seq = np.concatenate([patch[0, :-1], patch[:-1, -1],
                          patch[-1, :0:-1], patch[::-1, 0][:-1]])
    d = _wrap(np.diff(np.append(seq, seq[0])))
    return d.sum() / (2 * np.pi)


def run(report=print):
    report("CONFINED PHASE: Z3 disorder operator, exact lattice form "
           "(entry test)")
    report("=" * 70)
    N, L, size = 48, 24.0, 2.5

    psi0, dx, R, K2 = texture_state(N, L, amp=0.7, size=size)

    # positive control: the imprinted texture carries a charged vortex
    # ring (zeros of p1 on the f = pi/2, z = 0 circle). The detector
    # must see it, both as a loop winding and as a plaquette density
    chi0 = np.angle(rel_mode(psi0))
    w_ring = ring_winding(chi0, N)
    rho0_v, hits0 = vortex_density(chi0, R)
    report(f"  t = 0 (imprinted): ring winding = {w_ring:+.2f} "
           f"(exact -1 expected)")
    report(f"                     vortex-plaquette density = "
           f"{rho0_v:.2e} ({hits0} plaquettes)")

    # epochs
    psi2 = evolve(psi0.copy(), K2, R, int(2.0 / DT), g1)
    rho2_v, hits2 = vortex_density(np.angle(rel_mode(psi2)), R)
    psi8 = evolve(psi2, K2, R, int(6.0 / DT), g1)
    rho8_v, hits8 = vortex_density(np.angle(rel_mode(psi8)), R)
    psic = evolve(psi0.copy(), K2, R, int(8.0 / DT), 0.0)
    rhoc_v, hitsc = vortex_density(np.angle(rel_mode(psic)), R)
    report(f"  t = 2 (collapse):  density = {rho2_v:.2e} "
           f"({hits2} plaquettes)")
    report(f"  t = 8 (endpoint):  density = {rho8_v:.2e} "
           f"({hits8} plaquettes)")
    report(f"  t = 8 (g1 = 0):    density = {rhoc_v:.2e} "
           f"({hitsc} plaquettes)")

    report("-" * 70)
    # ── frozen assertions ────────────────────────────────────────────
    # 1. The detector works: the imprinted ring is seen exactly
    assert abs(abs(w_ring) - 1.0) < 0.05, "ring winding must be ±1"
    assert hits0 > 0, "the imprinted ring must register as plaquettes"
    # 2. REGISTERED OUTCOME (2026-07 measurement): the classical
    #    endpoint is defect-free: the disorder operator obeys a
    #    PERIMETER law (W = 1 on every contractible loop), i.e. NO
    #    confinement signature survives at the classical level.  The
    #    confined phase, if it exists, is a property of the QUANTUM
    #    confined sector, not of the classical substrate configuration
    #    (consistent with knot_charge.py. SUBSTRATE_CONJECTURES #1
    #    remains open, its obligation now sharply located).
    assert hits8 == 0, \
        f"classical endpoint must be defect-free (found {hits8})"
    report("  REGISTERED OUTCOME: charged-sector vortex content drains")
    report("  out with the unwinding. The classical endpoint is")
    report("  topologically TRIVIAL (perimeter law, rho_v = 0). The")
    report("  confined phase is not a classical-field property. The")
    report("  entry test's remaining obligation is the QUANTUM spectrum")
    report("  of the confined sector (gap + topological degeneracy).")
    return dict(w_ring=w_ring, rho=(rho0_v, rho2_v, rho8_v, rhoc_v),
                hits=(hits0, hits2, hits8, hitsc))


if __name__ == "__main__":
    run()
