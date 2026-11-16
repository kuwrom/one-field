"""The lump in all three components: what binds it, what makes it a clock.

lump_in_empty_field.py solves the lump on the clock orientation. This
probe examines the other directions: the second variation of E - mu N
on the lump in all six real components, the lump at g1 = 0, the
potential on every internal orientation, and the ground state at fixed
N reached from random three component data.

PREMISES. The nonrelativistic limit of the carrier (m = hbar = 1), the
self-bound sign of self_bound_vacuum.py (g0 = -1, g1 = 1/sqrt 2,
kappa = 1), the empty field outside the lump.

STATEMENTS (the numbers below check each one).

  S1 (second variation). With x = (Re psi, Im psi) and W(r) the real
     Hessian of w - mu |psi|^2 on the profile, the operator in the
     partial wave l is L_l = -d^2/dr^2 + l(l+1)/r^2 + W(r) on r x.
     On the lump: one negative eigenvalue (l = 0, the amplitude
     direction), one zero at l = 0 (the phase), one zero at l = 1
     (translation), and nothing else at or below zero. With dN/dmu < 0
     on the branch this is the stable count (Grillakis, Shatah,
     Strauss). Every relative direction is positive, so the lump is a
     minimum of the energy at fixed N against every change of internal
     orientation.
  S2 (what binds). At g1 = 0 the lump exists with a = g0/2 and its
     l = 0 kernel is five dimensional: the phase and the four relative
     rotations of U(3). The common term binds. The g1 term locks the
     orientation to the clock and opens the relative gaps.
  S3 (bounded). Over internal orientations the quartic coefficient is
     lowest on the clock, where it is a = g0/2 - g1/6. With the sextic
     the carrier potential U = (m c^2/2) R + w is nonnegative on every
     orientation iff m c^2 >= 2 |mu0|. Then the empty field is the
     global minimum.
  S4 (the clock is reached). In one dimension, gradient flow at fixed N
     from smooth random three component data ends on the clock lump.
     One Z3 Fourier mode carries the norm (q = 1, and q = 2 from the
     conjugate data, which is degenerate with it), and the energy equals
     that of the closed form lump.

Usage: python3 tests/probes/lump_stability.py     (about ten seconds)
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import interaction_matrix  # noqa: E402
from gravity_tensor import interaction_hessian, real_hessian  # noqa: E402
from arena_rulers import to_real  # noqa: E402
from self_bound_vacuum import energy as local_energy, hessian as fd_hessian  # noqa: E402
import lump_in_empty_field as lf  # noqa: E402

G0, KAPPA = lf.G0, lf.KAPPA
W3 = lf.W3


def local_hessian(f, mu, g1):
    """Real 6x6 Hessian of w - mu |psi|^2 at psi = f w^a / sqrt 3."""
    phi = f * W3 / math.sqrt(3)
    A, B = interaction_hessian(phi, G0, g1, mu)
    R = f * f
    A = A + KAPPA * (R * R * np.eye(3) + 2 * R * np.outer(phi, phi.conj()))
    B = B + 2 * KAPPA * R * np.outer(phi, phi)
    return real_hessian(A, B)


def partial_wave(profile, mu, g1, ell, rmax=30.0, n=900, k=14):
    """Lowest eigenvalues of L_l, each with its weight on the clock plane."""
    h = rmax / (n + 1)
    r = h * np.arange(1, n + 1)
    f = profile(r)
    # W depends on r through f alone: tabulate on amplitude, interpolate
    grid = np.linspace(0.0, f.max() * 1.001, 240)
    table = np.array([local_hessian(x, mu, g1) for x in grid])
    Wr = np.empty((n, 6, 6))
    for i in range(6):
        for j in range(6):
            Wr[:, i, j] = np.interp(f, grid, table[:, i, j])
    N = 6 * n
    main = np.zeros(N)
    blocks = {}
    for i in range(6):
        main[i::6] = 2 / h**2 + ell * (ell + 1) / r**2 + Wr[:, i, i]
        for j in range(i + 1, 6):
            blocks[j - i] = blocks.get(j - i, np.zeros(N))
            blocks[j - i][i::6] = Wr[:, i, j]
    off = -np.ones(N - 6) / h**2
    data, offsets = [main, off, off], [0, 6, -6]
    for d, v in blocks.items():
        data += [v[:N - d], v[:N - d]]
        offsets += [d, -d]
    L = diags(data, offsets, format="csc")
    vals, vecs = eigsh(L, k=k, sigma=-2.0, which="LM")
    clock = to_real(W3 / math.sqrt(3))
    plane = np.stack([clock, to_real(1j * W3 / math.sqrt(3))])
    share = []
    for v in vecs.T:
        u = v.reshape(n, 6)
        share.append(np.sum((u @ plane.T)**2) / np.sum(u * u))
    order = np.argsort(vals)
    return vals[order], np.array(share)[order]


def orientation_floor(g1, trials=4000, seed=0):
    """Lowest quartic coefficient w/R^2 over random internal orientations."""
    rng = np.random.default_rng(seed)
    best = np.inf
    for _ in range(trials):
        z = rng.normal(size=3) + 1j * rng.normal(size=3)
        z /= np.linalg.norm(z)
        best = min(best, local_energy(to_real(z), G0, g1, 0.0))
    return best


def flow_to_ground(N_target, seed, conjugate=False, L=60.0, n=384, dtau=0.05, steps=30000):
    """S4: 1D gradient flow at fixed N from smooth random data."""
    dx = L / n
    x = (np.arange(n) - n // 2) * dx
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    rng = np.random.default_rng(seed)
    white = rng.normal(size=(3, n)) + 1j * rng.normal(size=(3, n))
    psi = np.exp(-(x / 5)**2) * np.fft.ifft(np.fft.fft(white, axis=1) * (abs(k) < 2.0), axis=1)
    if conjugate:
        psi = psi.conj()

    def norm(p):
        return p * math.sqrt(N_target / (np.sum(abs(p)**2) * dx))
    psi = norm(psi)
    denom = 1 + 0.5 * dtau * k * k
    for _ in range(steps):
        R = np.sum(abs(psi)**2, axis=0)
        nl = np.einsum("xab,bx->ax", interaction_matrix(psi, G0, lf.G1), psi) + KAPPA * R**2 * psi
        psi = norm(np.fft.ifft((np.fft.fft(psi, axis=1) - dtau * np.fft.fft(nl, axis=1)) / denom,
                               axis=1))
    rho = abs(psi)**2
    R = rho.sum(axis=0)
    px = np.fft.ifft(1j * k * np.fft.fft(psi, axis=1), axis=1)
    E = np.sum(0.5 * abs(px)**2) * dx + np.sum(
        0.5 * G0 * R**2 + KAPPA * R**3 / 3 + lf.G1 * sum(
            rho[j] * np.real(np.conj(psi[(j + 1) % 3]) * psi[(j + 2) % 3])
            for j in range(3))) * dx
    modes = np.array([np.sum(np.conj(W3)[:, None]**q * psi, axis=0) / 3 for q in range(3)])
    weight = np.sum(abs(modes)**2, axis=1)
    return dict(E=E, mode=int(np.argmax(weight)), share=weight.max() / weight.sum())


def lump_1d(N_target):
    """Closed form 1D clock lump at given N: chemical potential and energy."""
    def integrals(mu):
        f, fp = lf.profile_1d(mu)

        def q(g):
            return quad(g, -60, 60, limit=400, epsabs=0, epsrel=1e-12)[0]
        return (q(lambda x: f(x)**2), q(lambda x: 0.5 * fp(x)**2),
                q(lambda x: lf.A * f(x)**4 + KAPPA * f(x)**6 / 3))
    mu = brentq(lambda m: integrals(m)[0] - N_target, 0.999 * lf.MU0, 1e-4 * lf.MU0, xtol=1e-13)
    N, K, P = integrals(mu)
    return mu, K + P


def run(report=print):
    rng = np.random.default_rng(4)
    # The analytic local Hessian against finite differences of the potential.
    f_test, mu_test = 0.83, -0.11
    fd = fd_hessian(lambda y: local_energy(y, G0, lf.G1, KAPPA) - mu_test * float(y @ y),
                    to_real(f_test * W3 / math.sqrt(3)))
    hess_err = np.max(abs(fd - local_hessian(f_test, mu_test, lf.G1)))
    assert hess_err < 1e-5

    # S1: the lump at mu = mu0/2.
    mu = 0.5 * lf.MU0
    d = lf.lump(mu)
    spectra = {ell: partial_wave(d["profile"], mu, lf.G1, ell) for ell in (0, 1, 2)}
    tol = 4e-3
    counts, rel_min = {}, {}
    for ell, (vals, share) in spectra.items():
        counts[ell] = (int(np.sum(vals < -tol)), int(np.sum(abs(vals) <= tol)))
        rel_min[ell] = float(vals[share < 0.5].min())
    assert counts == {0: (1, 1), 1: (0, 1), 2: (0, 0)}
    assert min(rel_min.values()) > 0.02
    threshold = -2 * mu                                # the continuum starts at -2 mu
    bound_rel = {ell: int(np.sum((s < 0.5) & (v < threshold)))
                 for ell, (v, s) in spectra.items()}
    dN = lf.lump(0.52 * lf.MU0)["N"] - lf.lump(0.48 * lf.MU0)["N"]
    assert dN > 0                                      # N grows as mu falls: dN/dmu < 0

    # S2: the same lump with the phase lock switched off.
    a0 = G0 / 2
    mu0_free = a0 * (-3 * a0 / (2 * KAPPA)) / 2
    d0 = lf.lump(0.5 * mu0_free, a=a0)
    v0, s0 = partial_wave(d0["profile"], 0.5 * mu0_free, 0.0, 0)
    kernel_free = int(np.sum(abs(v0) <= tol))
    assert int(np.sum(v0 < -tol)) == 1 and kernel_free == 5

    # S3: the clock is the lowest orientation. The empty field is the minimum.
    floor = orientation_floor(lf.G1)
    assert floor >= lf.A - 1e-9 and floor < lf.A + 0.02
    floor_free = orientation_floor(0.0)
    assert abs(floor_free - a0) < 1e-9

    def carrier_min(mc2):
        worst = np.inf
        for _ in range(2000):
            z = rng.normal(size=3) + 1j * rng.normal(size=3)
            z *= math.sqrt(rng.uniform(0, 3) * lf.R0) / np.linalg.norm(z)
            worst = min(worst, 0.5 * mc2 * float(np.vdot(z, z).real)
                        + local_energy(to_real(z), G0, lf.G1, KAPPA))
        clock = 0.5 * mc2 * lf.R0 + lf.A * lf.R0**2 + KAPPA * lf.R0**3 / 3
        return min(worst, clock)
    edge = 2 * abs(lf.MU0)
    assert carrier_min(1.02 * edge) > 0 > carrier_min(0.98 * edge)

    # S4: fixed-N ground state from random data.
    N1 = 1.26359
    mu1, E1 = lump_1d(N1)
    flows = [flow_to_ground(N1, 1), flow_to_ground(N1, 1, conjugate=True)]
    for fl in flows:
        assert fl["share"] > 1 - 1e-8 and abs(fl["E"] - E1) < 2e-4 * abs(E1)
    assert {flows[0]["mode"], flows[1]["mode"]} == {1, 2}

    report(f"Local Hessian against finite differences: {hess_err:.1e}.")
    for ell, (vals, share) in spectra.items():
        report(f"S1 l = {ell}: lowest eigenvalues {np.round(vals[:5], 4)}; (negative, zero) = "
               f"{counts[ell]}. Lowest relative eigenvalue {rel_min[ell]:.4f}. Relative "
               f"states below the continuum edge {threshold:.4f}: {bound_rel[ell]}.")
    report(f"S1 branch: N(0.52 mu0) - N(0.48 mu0) = {dN:+.2f}, dN/dmu < 0.")
    report(f"S2 g1 = 0: lump at a = {a0:.3f}, N = {d0['N']:.2f}; l = 0 lowest eigenvalues "
           f"{np.round(v0[:7], 4)}. Kernel dimension {kernel_free} (phase + four rotations).")
    report(f"S3 quartic floor over orientations: {floor:.5f} (clock a = {lf.A:.5f}). At g1 = 0: "
           f"{floor_free:.5f}. U >= 0 on every orientation iff m c^2 >= {edge:.5f}.")
    for seed, fl in zip(("random data", "conjugate data"), flows):
        report(f"S4 {seed}: Z3 mode q = {fl['mode']}, share {fl['share']:.10f}, "
               f"E = {fl['E']:.7f} (closed form {E1:.7f}).")
    report("RESULT: the lump has the stable count (one negative direction and the two "
           "symmetry zeros) and is a minimum against every change of internal "
           "orientation. The common term binds it. The g1 term makes it a clock. The "
           "empty field is the global minimum for m c^2 >= 2|mu0|, and random data at "
           "fixed N flows to the clock lump.")
    return dict(hessian_error=hess_err, counts=counts, relative_min=rel_min,
                bound_relative=bound_rel, kernel_free=kernel_free, floor=floor,
                floor_free=floor_free, edge=edge, flows=flows, E_1d=E1, mu_1d=mu1)


if __name__ == "__main__":
    run()
