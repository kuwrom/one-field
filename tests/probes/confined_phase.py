"""
confined_phase.py: the layer-bridge test (substrate verification
ledger entry #1. Derivation program D6 of the gaps audit).

SCOPE.  The web's statistics/Skyrmion chain routes through a
CONFINED TOPOLOGICAL PHASE of the substrate (unbroken Z3 centre): the
WZ citation, winding quantization, and the substrate-web bridge all
presuppose it, and knot_charge.py showed the substrate-only theory protects
nothing: protection lives in the confined phase if it lives
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

STATUS.  The entry test made executable: a defined order parameter,
measured exactly, with a positive control and a coupling control,
frozen in CI.  The Gaussian excitation gap and classical area law are
measured.  Quantum topological-order diagnostics use spectra on distinct
backgrounds and finite-size scaling.

SYNTHESIS (the four-probe chain that carries the entry test):
  * effective_action.py: the stabilizing quartic exists given the gap.
  * knot_spectrum.py: the gap is measured at Gaussian order (one gapped
    relative line: omega = 1.178 = 1.67 |b_1| = 0.267 mu. Bogoliubov
    linear response = the quantum quasiparticle spectrum at that order).
  * this probe: the classical disorder signature is measured (area law).
  * substrate_depth.py: the depth rho xi^3 is bounded from below.
Beyond Gaussian order a node event is a finite-action excursion,
suppressed as exp(-c rho xi^3). The diluteness rho xi^3 is a
substrate-depth parameter the web does not fix.  READ AS A PREDICTION:
stable fermionic matter exists, so through the closed implication chain
the framework requires the substrate to be a deep condensate
(rho xi^3 >> 1). Any independent handle on substrate depth must land
there.  Killed by: a derivation forcing rho xi^3 ~ 1, or the
charged-sector gap closing beyond Gaussian order in the saturated
regime.

MEASURED.  Under the non-variational variant of the equation, the
classical endpoint is DEFECT-FREE
(rho_v = 0, perimeter law), which reads as the confined phase not
being a property of the classical substrate at all.  That variant
drops the exchange terms, which ARE the
inter-component population current (see tests/substrate/nls_soliton.py).
With them present the same measurement inverts: the plaquette count
grows 36 (t = 0, the imprint) -> 52 (t = 2) -> 517 (t = 8),
rho_v ~ 1e-2, while the g1 = 0 control still ends at exactly 0.  A
candidate classical disorder signature exists after all.  And the AREA
LAW is now measured (area_law() below, mark slow): the Z3 disorder
operator with fractional test charge 1/3 falls as ln W = -c a^2
(c = 0.021, residual 8x better than the perimeter fit, W from 0.93 to
0.12 over a = 2..10 at the saturated endpoint).  The classical
disorder signature is complete: vortex proliferation plus area law.
Spectra on distinct backgrounds and box-size scaling are the corresponding
quantum topological-order diagnostics.

Usage: python3 confined_phase.py     (~7 min single-core, mark slow)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _bootstrap  # noqa: E402,F401
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
    ## frozen assertions ##
    # 1. The detector works: the imprinted ring is seen exactly
    assert abs(abs(w_ring) - 1.0) < 0.05, "ring winding must be ±1"
    assert hits0 > 0, "the imprinted ring must register as plaquettes"
    # 2. The energy-conserving equation produces measured plaquette
    #    counts 36 (t=0), 52 (t=2), and 517 (t=8), with rho_v ~ 1e-2
    #    at the endpoint.  The g1 = 0 control ends at exactly 0.
    #    The assertions freeze the direction of the effect because the
    #    exact count depends on resolution and evolution time.
    assert hits8 > 0, \
        f"true-coupling endpoint must carry defects (found {hits8})"
    assert hitsc == 0, \
        f"g1 = 0 control must stay defect-free (found {hitsc})"
    assert hits8 > hits0, \
        f"defects must proliferate, not drain ({hits0} -> {hits8})"
    report("  OUTCOME: the")
    report("  charged-sector vortex content PROLIFERATES under the true")
    report(f"  coupling ({hits0} -> {hits8} plaquettes, rho_v = "
           f"{rho8_v:.1e}) while the g1 = 0")
    report("  control ends at exactly 0.  The non-variational 'defect-free")
    report("  classical endpoint' was an artefact of the non-variational")
    report("  equation, which dropped the exchange terms carrying the")
    report("  inter-component current.  A classical disorder signature")
    report("  is therefore back on the table; an AREA law still needs")
    report("  the rho_v scaling with loop and box size, which this probe")
    report("  does not measure.")
    return dict(w_ring=w_ring, rho=(rho0_v, rho2_v, rho8_v, rhoc_v),
                hits=(hits0, hits2, hits8, hitsc))


if __name__ == "__main__":
    run()


def area_law(report=print):
    """Conjecture #1's second half: the disorder operator W(C) versus
    loop size at the saturated endpoint.  W(a) = |<exp(i circ_C)>|
    over square loops of side a in the three coordinate planes,
    centres in the core.  Dilute-vortex-gas prediction: area law
    ln W ~ -c a^2 when the charged vortices proliferate; perimeter
    law ln W ~ -c' a when they do not.  The measured exponent is the
    verdict; this function reports the fits and returns them."""
    N, L, size = 48, 24.0, 2.5
    psi0, dx, R, K2 = texture_state(N, L, amp=0.7, size=size)
    cache = "/tmp/one_field_chi8.npy"
    if os.path.exists(cache):
        chi = np.load(cache)
    else:
        psi8 = evolve(psi0.copy(), K2, R, int(8.0 / DT), g1)
        chi = np.angle(rel_mode(psi8))
        try:
            np.save(cache, chi)
        except OSError:
            pass

    sides = list(range(2, 11))
    W = []
    for a in sides:
        vals = []
        for ax1, ax2 in ((0, 1), (0, 2), (1, 2)):
            d1 = _wrap(np.diff(chi, axis=ax1))
            d2 = _wrap(np.diff(chi, axis=ax2))
            n = N - a
            # loop sums via cumulative differences along the perimeter
            for i in range(0, n, 2):
                for j in range(0, n, 2):
                    for k in range(0, N, 4):
                        idx = [0, 0, 0]
                        idx[ax1], idx[ax2] = i, j
                        rest = 3 - ax1 - ax2
                        idx[rest] = k
                        ci, cj = i + a // 2, j + a // 2
                        cc = [0, 0, 0]
                        cc[ax1], cc[ax2], cc[rest] = ci, cj, k
                        if R[tuple(cc)] > 5.0:
                            continue
                        circ = 0.0
                        for m in range(a):
                            p = list(idx); p[ax1] = i + m
                            circ += d1[tuple(p)]
                            p = list(idx); p[ax1] = i + a; p[ax2] = j + m
                            circ += d2[tuple(p)]
                            p = list(idx); p[ax1] = i + a - 1 - m + 1 - 1
                            p[ax1] = i + (a - 1 - m); p[ax2] = j + a
                            circ -= d1[tuple(p)]
                            p = list(idx); p[ax2] = j + (a - 1 - m)
                            circ -= d2[tuple(p)]
                        # closed-loop circulation is exactly
                        # 2 pi n_enc; the Z3 disorder operator uses
                        # the fractional test charge 1/3:
                        n_enc = int(round(circ / (2 * np.pi)))
                        vals.append(np.exp(2j * np.pi * n_enc / 3.0))
        W.append(abs(np.mean(vals)) if vals else 1.0)
    W = np.array(W)
    lnW = np.log(np.maximum(W, 1e-12))
    A2 = np.array([a * a for a in sides], float)
    A1 = np.array(sides, float)
    # least-squares slopes through the origin, and fit residuals
    c_area = float(-(A2 @ lnW) / (A2 @ A2))
    c_perim = float(-(A1 @ lnW) / (A1 @ A1))
    res_area = float(np.linalg.norm(lnW + c_area * A2))
    res_perim = float(np.linalg.norm(lnW + c_perim * A1))
    report("  AREA LAW at the saturated endpoint (t = 8):")
    for a, w in zip(sides, W):
        report(f"    a = {a:2d}: W = {w:.4f}, ln W = {float(np.log(max(w, 1e-12))):+.4f}")
    report(f"    area fit c = {c_area:.5f} (residual {res_area:.3f}); "
           f"perimeter fit c' = {c_perim:.5f} (residual {res_perim:.3f})")
    report(f"    verdict: {'AREA' if res_area < res_perim else 'PERIMETER'} "
           f"law fits better")
    # frozen assertions (measured: c_area = 0.0207, residual
    # 0.113 vs perimeter 0.890; W falls 0.93 -> 0.12 over a = 2..10)
    assert res_area < 0.5 * res_perim, "area law must fit decisively"
    assert 0.01 < c_area < 0.04
    assert all(W[i + 1] < W[i] + 1e-6 for i in range(len(W) - 1))
    assert W[-1] < 0.2
    return dict(sides=sides, W=[float(w) for w in W],
                c_area=c_area, c_perim=c_perim,
                res_area=res_area, res_perim=res_perim,
                area_wins=bool(res_area < res_perim))
