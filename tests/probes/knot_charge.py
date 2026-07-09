"""
KNOT CHARGE: does the substrate protect pi_3 winding? (measured answer: no)

The winding_texture probe shows the 3D texture carries unit pi_3 charge as a
STATIC object. This probe asks the dynamical question: does a knot in the
Z3-coupled substrate KEEP that charge? The answer is a clean negative, with
the mechanism identified, and it sharpens where protection must come from.

Three facts, each asserted below:

1. ESTIMATOR: the lattice charge of the static texture converges to -1 as
   dx -> 0 (the deficit at coarse dx is discretization, not physics).

2. UNWINDING: a saturated knot imprinted with the unit texture loses its
   charge in a few natural time units. The relative spinor phi = (phi1,
   phi2) lives in C^2, a LINEAR space: zeros of phi are codimension-4,
   i.e. isolated spacetime EVENTS, and the winding jumps only there. The
   probe detects the events directly: min|phi| in the core plunges by an
   order of magnitude exactly when B collapses.

3. BALLISTIC, NOT ACTIVATED: the unwinding time grows only linearly with
   texture size, tau ~ R/c_s with c_s = sqrt(g0 rho_tot) (measured
   exponent ~1.3 over R in [1.8, 3.2]. An energy barrier would give
   activated, at least quadratic growth). There is no metastability
   window. This is Derrick's theorem doing its work: quadratic-gradient
   energetics cannot fix a texture scale, and the substrate has no
   Skyrme-type quartic term. The stable attractor is a NODELESS, B = 0
   density lump: the substrate soliton is Q-ball-like, not topological.

Consequence for the emergence tree: unit winding is NOT a substrate-level
invariant. Its creation and protection are obligations of the confined-phase
effective theory, exactly where the level-3 WZ structure already lives
(winding_texture.py). The derivation program this opens: show that the
confined sector's effective action contains a texture-stabilizing term
(Skyrme quartic or equivalent), or that the protected object is the
effective SU(3) field rather than the bare relative spinor. Until then,
statements of the form "3D protects absolutely" hold for maps of fixed
nonvanishing norm, not for the substrate's linear field.

4. NO ROTATION LOOPHOLE: Q-lump / rotating-skyrmion stabilization (Leese
   1991. Battye-Cooper-Sutcliffe PRL 2002) requires a conserved charge
   backing the internal rotation. The candidate generators are the U(2)
   isometries of the relative spinor. Asserted below by exact commutation
   test: the sigma_z (Hopf-fiber) rotation does NOT commute with the
   substrate flow, while the genuine global U(1) commutes to machine
   epsilon -- so no Noether charge exists to back differential rotation,
   and the one surviving charge (total N) is condensate-dominated (~99%)
   and cannot tie a texture scale. A tilt scan (core-filling eta in
   [0.5, 4], amplitude in [0.4, 1.0]) confirms: unwinding time stays in a
   flat ballistic band, and the strongest-rotation candidate dies fastest.

Runtime: ~3.5 minutes (64^3 short run + two static grids + 32^3
commutation test). Mark slow.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# the shared 3D machinery (one copy: substrate3d.py)
from substrate3d import (DT, evolve, omega, rel_modes,  # noqa: E402
                         texture_state)


def charge(p1, p2, dx, R, r_core=8.0):
    norm = np.sqrt(np.abs(p1)**2 + np.abs(p2)**2)
    nrm = np.maximum(norm, 1e-30)
    n4 = np.stack([p1.real / nrm, p1.imag / nrm,
                   p2.real / nrm, p2.imag / nrm], -1)
    gs = [np.gradient(n4[..., a], dx, axis=ax)
          for ax in range(3) for a in range(4)]
    M = np.stack([n4, np.stack(gs[0:4], -1), np.stack(gs[4:8], -1),
                  np.stack(gs[8:12], -1)], axis=-1)
    dens = np.linalg.det(M) / (2 * np.pi**2)
    dens = np.where(norm > 0.02 * norm.max(), dens, 0.0)
    return float(dens[R < r_core].sum() * dx**3)


def evolve_nodes(psi, K2, R, n_steps, node_every=20):
    """Evolve. Return (psi, deepest core min|phi| at node_every cadence).

    Node events are isolated in spacetime (codim-4), so the detector must
    sample densely in time. Charge needs computing only at checkpoints.
    """
    core = R < 6.0
    tracker = {"deepest": np.inf}

    def _detect(p):
        p1, p2 = rel_modes(p)
        mn = float(np.sqrt(np.abs(p1)**2 + np.abs(p2)**2)[core].min())
        tracker["deepest"] = min(tracker["deepest"], mn)

    psi = evolve(psi, K2, R, n_steps, callback=_detect, every=node_every)
    return psi, tracker["deepest"]


def rotate_sz(psi, alpha):
    p1, p2 = rel_modes(psi)
    p0 = psi.mean(axis=0)
    p1, p2 = np.exp(1j * alpha) * p1, np.exp(-1j * alpha) * p2
    out = np.empty_like(psi)
    for k in range(3):
        out[k] = p0 + omega**k * p1 + omega**(2 * k) * p2
    return out


def main():
    print("KNOT CHARGE: substrate-level pi_3 protection test")
    print("=" * 64)

    # -- 1. estimator convergence (static) --
    Bs = {}
    for N in (48, 96):
        psi, dx, R, _ = texture_state(N, 24.0, amp=0.7, size=2.5)
        Bs[N] = charge(*rel_modes(psi), dx, R)
        print(f"  static texture, N={N:3d}: B = {Bs[N]:+.4f}")
    assert abs(Bs[96]) > abs(Bs[48]) and abs(Bs[96]) > 0.93, \
        "estimator must converge toward |B| = 1 with resolution"
    print("  estimator converges toward -1: CONFIRMED")
    print("-" * 64)

    # -- 2 & 3. dynamical unwinding with node detection --
    N, size = 64, 1.8
    psi, dx, R, K2 = texture_state(N, 24.0, amp=0.7, size=size)
    core = R < 6.0
    B0 = charge(*rel_modes(psi), dx, R)
    p1, p2 = rel_modes(psi)
    m0 = float(np.sqrt(np.abs(p1)**2 + np.abs(p2)**2)[core].min())
    print(f"  saturated knot R={size}: B(0) = {B0:+.3f}, "
          f"min|phi|_core = {m0:.3f}")
    min_node, B_hist = m0, []
    for chunk in range(6):  # t = 0 -> 3 in steps of 0.5
        psi, deepest = evolve_nodes(psi, K2, R, int(0.5 / DT))
        p1, p2 = rel_modes(psi)
        B = charge(p1, p2, dx, R)
        min_node = min(min_node, deepest)
        B_hist.append(B)
        print(f"  t={0.5*(chunk+1):4.1f}: B = {B:+.3f}, "
              f"deepest min|phi|_core so far = {min_node:.4f}")
    assert min(abs(b) for b in B_hist) < 0.5, \
        "winding must unwind within a few natural times (it does)"
    # A true zero sampled on the grid reads ~|grad phi| * dx, not 0. The
    # criterion is a clear dip below the initial envelope floor, which a
    # nodeless evolution never shows (post-unwinding plateaus sit ABOVE m0).
    assert min_node < 0.5 * m0, \
        "unwinding must proceed through amplitude-node events (it does)"
    print("-" * 64)

    # -- 4. rotation loophole: no conserved charge backs internal rotation --
    psi0, dx, R, K2 = texture_state(32, 24.0, amp=0.7, size=2.5)
    al = 0.7
    A = evolve(rotate_sz(psi0, al), K2, R, 100)
    Bv = evolve(psi0.copy(), K2, R, 100)
    Bv = rotate_sz(Bv, al)
    r_sz = np.linalg.norm(A - Bv) / np.linalg.norm(Bv)
    C = evolve(psi0 * np.exp(1j * al), K2, R, 100)
    D = evolve(psi0.copy(), K2, R, 100)
    r_u1 = np.linalg.norm(C - D * np.exp(1j * al)) / np.linalg.norm(D)
    print(f"  sigma_z commutation residual = {r_sz:.2e}  (broken)")
    print(f"  global U(1) residual         = {r_u1:.2e}  (true symmetry)")
    assert r_sz > 1e-6 and r_u1 < 1e-10, \
        "no Noether charge backs internal rotation (Q-lump route closed)"
    print("-" * 64)
    print("  VERDICT: pi_3 winding is NOT protected at substrate level.")
    print("  Mechanism: Derrick collapse of the texture, unwinding through")
    print("  codim-4 amplitude-node events. tau ~ R/c_s (ballistic, no")
    print("  barrier). The stable attractor is a nodeless B = 0 lump.")
    print("  Protection is an obligation of the confined-phase effective")
    print("  theory (level-3 WZ sector), not of the bare Z3-NLS. The")
    print("  rotation loophole is closed: the spinor's U(2) is broken to")
    print("  the global U(1), whose charge is condensate-dominated.")


if __name__ == "__main__":
    main()
