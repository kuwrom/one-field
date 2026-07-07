"""
winding_texture.py (statistics): the explicit 3D knot and
its topological charge, computed.

THE CONSTRUCTION: the confined child is SU(3); 3D textures are
classified by pi_3(SU(3)) = Z.  Build the explicit unit hedgehog
(SU(2) hedgehog embedded in the upper-left block of SU(3)):

    U(x) = exp(i f(r) x_hat . sigma),  f(0) = pi, f(inf) = 0,

and compute its winding number by the standard integral

    B = (1/24 pi^2) INT eps_{ijk} Tr[(U+ dU_i)(U+ dU_j)(U+ dU_k)].

RESULT: B = 1 numerically.  This is the knot as a genuine
topological object, the thing 1D could not protect (we watched
windings mix there), 3D protects absolutely.

THE STATISTICS STEP (Witten 1983, conditional as labeled in
notes/statistics.md): with the level-3 WZ term supplied by the
CS/WZW correspondence of the confined sector, a B = 1 texture
acquires exchange phase e^{i pi k B} = e^{3 pi i} = -1: FERMION.
What is computed here: the texture EXISTS and carries unit charge.
What remains cited: the WZ term's presence at level 3 in the
effective action (the CS/WZW correspondence supplies it; a
first-principles substrate derivation is the labeled remainder).

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
    report("  With the level-3 WZ term (CS/WZW of the confined sector),")
    report("  exchange phase = e^{i pi * 3 * B} = e^{3 pi i} = -1:")
    report("  the B = 1 knot is a FERMION.  Computed here: existence +")
    report("  unit charge.  Cited (labeled): the WZ term at level 3.")
    report("  Contrast with 1D (probes/stationary.py): there the")
    report("  winding MIXED, no pi_3 protection; in 3D it cannot.")
    return B


if __name__ == "__main__":
    run()
