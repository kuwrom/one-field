"""
winding_texture.py (statistics): the explicit 3D knot and
its topological charge, computed.

THE CONSTRUCTION: the confined child is SU(3). 3D textures are
classified by pi_3(SU(3)) = Z.  Build the explicit unit hedgehog
(SU(2) hedgehog embedded in the upper-left block of SU(3)):

    U(x) = exp(i f(r) x_hat . sigma),  f(0) = pi, f(inf) = 0,

and compute its winding number by the standard integral

    B = (1/24 pi^2) INT eps_{ijk} Tr[(U+ dU_i)(U+ dU_j)(U+ dU_k)].

RESULT: B = 1 numerically.  This is the knot as a genuine
topological object, the thing 1D could not protect (we watched
windings mix there).  In 3D the charge is protected for maps of
fixed nonvanishing norm. The substrate's relative spinor is a
LINEAR C^2 field, not such a map, and does not protect it.
Measured: tests/probes/knot_charge.py (ballistic Derrick collapse,
unwinding through codim-4 amplitude-node events).  Protection is a
property claimed of the confined phase (registry:
DERIVATION_PROGRAMS #4).

THE STATISTICS STEP (Witten 1983. Conditional as labeled in the
registry's CLOSURE_RECORD and DERIVATION_PROGRAMS #4).  The honest
split: the WZ term's LEVEL is derived by counting the three Dirac
species ('t Hooft anomaly matching, tests/probes/bookkeeping.py),
and its PRESENCE is forced by the same matching GIVEN the confined
phase (the anomaly must be reproduced in the IR, and the WZ term is
the unique local functional that does it. Corroboration-grade
citation).  With the level-3 WZ term, a B = 1 texture acquires
exchange phase e^{i pi k B} = e^{3 pi i} = -1: FERMION.  What is
computed here: the texture EXISTS and carries unit charge.  The
single open premise is the confined phase itself
(SUBSTRATE_CONJECTURES #1).

Usage: python3 winding_texture.py   (~15 s)
"""

import numpy as np


def run(report=print):
    n = 42
    L = 6.0
    xs = np.linspace(-L, L, n)
    dx = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12
    f = np.pi * np.exp(-R / 1.5)          # f(0)=pi, f(inf)->0 (smooth)

    # SU(2) hedgehog embedded in SU(3): U = cos f + i sin f (x_hat.sigma)
    nx, ny, nz = X/R, Y/R, Z/R
    c, s = np.cos(f), np.sin(f)
    # 2x2 block entries
    U = np.zeros((n, n, n, 3, 3), dtype=complex)
    U[..., 0, 0] = c + 1j*s*nz
    U[..., 0, 1] = s*(1j*nx + ny)
    U[..., 1, 0] = s*(1j*nx - ny)
    U[..., 1, 1] = c - 1j*s*nz
    U[..., 2, 2] = 1.0

    def d(A, axis):
        return np.gradient(A, dx, axis=axis)

    Ud = np.conj(np.swapaxes(U, -1, -2))
    Ls = []
    for ax in range(3):
        dU = np.stack([d(U[..., i, j], ax) for i in range(3)
                       for j in range(3)], axis=-1).reshape(n, n, n, 3, 3)
        Ls.append(np.einsum("...ij,...jk->...ik", Ud, dU))
    L0, L1, L2 = Ls
    comm = (np.einsum("...ij,...jk,...kl->...il", L0, L1, L2)
            - np.einsum("...ij,...jk,...kl->...il", L0, L2, L1)
            + np.einsum("...ij,...jk,...kl->...il", L1, L2, L0)
            - np.einsum("...ij,...jk,...kl->...il", L1, L0, L2)
            + np.einsum("...ij,...jk,...kl->...il", L2, L0, L1)
            - np.einsum("...ij,...jk,...kl->...il", L2, L1, L0))
    dens = np.einsum("...ii->...", comm).real / 6.0   # eps contraction /3!
    B = float(dens.sum() * dx**3 / (24 * np.pi**2) * 6.0)

    report("WINDING TEXTURE: the 3D knot and its charge")
    report("=" * 64)
    report(f"  grid {n}^3, box +/-{L}, profile f(0) = pi -> 0")
    report(f"  topological charge B = {B:+.4f}   (target: +1, "
           f"pi_3(SU(3)) = Z)")
    ok = abs(abs(B) - 1.0) < 0.15   # sign = hedgehog orientation convention
    report(f"  unit winding: {'CONFIRMED' if ok else 'grid too coarse'}")
    report("-" * 64)
    report("  With the level-3 WZ term (presence cited: CS/WZW. Level")
    report("  derived by counting: bookkeeping.py), exchange phase =")
    report("  e^{i pi * 3 * B} = e^{3 pi i} = -1: the B = 1 knot is a")
    report("  FERMION.  Computed here: existence + unit charge.")
    report("  Protection: holds for fixed-norm maps. The substrate's")
    report("  linear spinor is not one and does NOT protect it")
    report("  (knot_charge.py).  Protection is claimed of the confined")
    report("  phase (DERIVATION_PROGRAMS #4).")
    return B


if __name__ == "__main__":
    run()
