"""Clock and ruler fields on the Z3-NLS substrate: the arena count.

A metric independent arena in the sense of Carone, Erlich and Vaman
(arXiv:1610.08521, Eqs. 2.1-2.11) needs D scalars X^I with a constant
target metric eta_IJ, entering only through derivatives (the global
shift X^I -> X^I + Delta^I stated below their Eq. 2.1), so that the
static gauge X^I ~ x^I (Eq. 2.7) solves the field equations with the
other fields at the minimum of V (text after Eq. 2.11). X^0 is the
clock, X^1..X^{D-1} the rulers. In 3+1 dimensions D = 4.

PREMISES. The committed functional of tests/substrate/nls_soliton.py,
alpha = 0, g1 = g0/sqrt(2); V = (g0/2) R^2 + g1 sum_k rho_k
Re(psi*_{k+1} psi_{k+2}), R = rho_tot. Flat kinetic term. Clock
phi_a = sqrt(rho) w^a (R = 3 rho), mu = (3 g0 - g1) rho. Hessians are
taken in (Re psi, Im psi) of V - mu R.

LEMMA 1 (linearity). A continuous internal point symmetry of
(1/2)|grad psi|^2 + V is an isometry of the flat target R^6. Since
sum_k rho_k sqrt(rho_{k+1} rho_{k+2}) <= R^2/3, V >= (g0/2 - g1/3)
R^2, a positive definite quartic for g1 < 3 g0/2. Every symmetry
then fixes psi = 0 and lies in O(6) (in U(3) once the first order
Schrodinger term is kept).

LEMMA 2 (rank). A ruler shift must commute with the other shifts
(abelian Delta^I, flat eta_IJ). Commuting antisymmetric 6x6 matrices
span at most 3 dimensions (maximal torus of SO(6), and of U(3)).
Hence any potential in this class admits at most three independent
shift symmetric fields: D <= 3, one clock and at most two rulers.

THEOREM (committed potential). The phases enter V only through
cos(phi_{k+2} - phi_{k+1}) with coefficients g1 rho_k sqrt(rho_{k+1}
rho_{k+2}) that are nonzero on the clock, so the symmetry algebra of V
inside so(6) is the single common phase u(1): one clock, zero rulers.
Proof of "no ruler direction": X^i = x^i sweeps the field along a
whole orbit, so V must be constant along that orbit. Such an orbit
consists of critical points of V - mu R, so its tangent is a zero
mode of the Hessian. The relative sector Hessian eigenvalues are
(5 -/+ sqrt 13) g1 rho / 2, each twice, all positive. There is no flat
direction that could become a ruler.

LEMMA 3 (pairing, any Z3 x U(1)-invariant V). g' = w^-1 (cyclic
shift) fixes the clock and acts on the four relative directions with
no fixed vector. J = (2 g' + 1)/sqrt 3 is a complex structure that
commutes with the Hessian, so relative zero modes come in pairs. The
radial stiffness is positive, so rulers are relative zero modes. Since
they come in pairs, three rulers force all four relative directions to
be zero modes: the relative block vanishes and the phases unlock.
Locked (rank 4, committed): zero rulers. One pair unlocked: at most
two. Locking and three rulers are strictly incompatible.

SCOPE. This closes the clock and ruler (static gauge) form only.
Carone (arXiv:1911.13166, Sec. III.A, IV, Eqs. 3.2-3.4, 4.1-4.2): the
massless pole sits where the X^I vev vanishes, and a determinant
action with no X gives the same pole. That form is the carrier written
without a background metric (tests/probes/composite_metric.py), and
zero rulers is its massless point.

RESULT. The continuous symmetry algebra of the committed V inside
so(6) is one dimensional, the common phase. The relative Hessian
block, (5 -/+ sqrt 13) g1 rho / 2 each twice, has no zero, so the
computed arena is one clock and zero rulers. The U(1)^3 extension
gives two relative zero modes, and three rulers are excluded.

Usage: python3 tests/probes/arena_rulers.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import G1_OVER_G0, interaction_matrix  # noqa: E402

W = np.exp(2j * np.pi / 3)
SHIFT = np.roll(np.eye(3), 1, axis=1)          # (S psi)_k = psi_{k+1}


def to_complex(x):
    return x[:3] + 1j * x[3:]


def to_real(z):
    return np.r_[z.real, z.imag]


def real_rep(U):
    return np.block([[U.real, -U.imag], [U.imag, U.real]])


def potential(x, g0=1.0, g1=G1_OVER_G0):
    psi = to_complex(x)
    rho = abs(psi)**2
    ex = sum(rho[k] * np.real(psi[(k+1) % 3].conj() * psi[(k+2) % 3])
             for k in range(3))
    return 0.5 * g0 * rho.sum()**2 + g1 * ex


def gradient_live(x, g0=1.0, g1=G1_OVER_G0):
    """Real gradient from the committed force M[psi] psi = dV/dpsi*."""
    psi = to_complex(x)
    return 2 * to_real(interaction_matrix(psi, g0, g1) @ psi)


def fd_gradient(f, x, h=1e-6):
    return np.array([(f(x + h*e) - f(x - h*e)) / (2*h) for e in np.eye(len(x))])


def fd_hessian(f, x, h=1e-4):
    n, e = len(x), np.eye(len(x))
    H = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            H[i, j] = (f(x + h*e[i] + h*e[j]) - f(x + h*e[i] - h*e[j])
                       - f(x - h*e[i] + h*e[j]) + f(x - h*e[i] - h*e[j])) / (4*h*h)
    return 0.5 * (H + H.T)


def clock_frame(phi):
    """Radial, phase, and orthonormal relative (4-dim) directions."""
    er, et = to_real(phi), to_real(1j * phi)
    er, et = er / np.linalg.norm(er), et / np.linalg.norm(et)
    P = np.eye(6) - np.outer(er, er) - np.outer(et, et)
    w, v = np.linalg.eigh(P)
    return er, et, v[:, w > 0.5]


def stationary_hessian(f, phi):
    """Hessian of f - mu |psi|^2 at phi, mu from the radial gradient."""
    x = to_real(phi)
    g = fd_gradient(f, x)
    er, et, Q = clock_frame(phi)
    mu = (g @ er) / (2 * np.linalg.norm(x))
    residual = np.linalg.norm(g - 2 * mu * x)
    return fd_hessian(f, x) - 2 * mu * np.eye(6), mu, residual


def symmetry_algebra_dim(f, samples=40, seed=1):
    """dim of {T in so(6): grad f(x) . T x = 0 for all x}."""
    rng = np.random.default_rng(seed)
    E = np.eye(6)
    gens = [np.outer(E[i], E[j]) - np.outer(E[j], E[i])
            for i in range(6) for j in range(i + 1, 6)]
    rows = [[fd_gradient(f, x) @ (T @ x) for T in gens]
            for x in rng.normal(size=(samples, 6))]
    s = np.linalg.svd(np.array(rows), compute_uv=False)
    null = int(np.sum(s < 1e-7 * s[0]))
    return null, s


def random_invariant_potential(rng, g0=1.0):
    """A generic Z3 x U(1)-invariant quartic (cyclic average, charge 0)."""
    C = rng.normal(size=(3,) * 4) + 1j * rng.normal(size=(3,) * 4)

    def f(x):
        out = 0.5 * g0 * np.sum(abs(to_complex(x))**2)**2
        for k in range(3):
            p = np.linalg.matrix_power(SHIFT, k) @ to_complex(x)
            out += 0.3 * np.real(np.einsum("abcd,a,b,c,d", C, p.conj(),
                                           p.conj(), p, p))
        return out
    return f


def run(report=print):
    g0, g1, rho = 1.0, G1_OVER_G0, 1.0
    r, mu_clock = g1 * rho, (3 * g0 - g1) * rho
    phi = math.sqrt(rho) * W**np.arange(3)
    rng = np.random.default_rng(7)

    # V is the committed functional: its gradient is the live force.
    x = rng.normal(size=6)
    force_error = np.max(abs(fd_gradient(potential, x) - gradient_live(x)))
    assert force_error < 1e-6

    # Continuous internal symmetries of V inside so(6): the u(1) alone.
    sym_dim, _ = symmetry_algebra_dim(potential)
    assert sym_dim == 1
    assert abs(fd_gradient(potential, x) @ (real_rep(1j*np.eye(3)) @ x)) < 1e-6

    # Clock Hessian: full six, fixed rho_tot five, relative four.
    H, mu, residual = stationary_hessian(potential, phi)
    assert abs(mu - mu_clock) < 1e-6 and residual < 1e-6
    er, et, Q = clock_frame(phi)
    full = np.linalg.eigvalsh(H)
    sphere_basis = np.c_[et, Q]
    sphere = np.linalg.eigvalsh(sphere_basis.T @ H @ sphere_basis)
    Hrel = Q.T @ H @ Q
    rel = np.linalg.eigvalsh(Hrel)
    expected = np.repeat([(5 - math.sqrt(13)) * r / 2,
                          (5 + math.sqrt(13)) * r / 2], 2)
    assert np.max(abs(rel - expected)) < 1e-5
    assert abs(er @ H @ er - 4 * mu_clock) < 1e-5
    assert np.sum(abs(sphere) < 1e-5) == 1          # the clock direction
    assert np.min(rel) > 0.4                        # no ruler direction

    # Z3' complex structure on the relative block (Lemma 3).
    G = real_rep(np.conj(W) * SHIFT)
    assert np.linalg.norm(G @ to_real(phi) - to_real(phi)) < 1e-12
    Grel = Q.T @ G @ Q
    assert np.linalg.norm(Q @ Grel - G @ Q) < 1e-12   # preserves the block
    assert np.min(abs(np.linalg.eigvals(Grel) - 1)) > 0.5   # no fixed vector
    J = (2 * Grel + np.eye(4)) / math.sqrt(3)
    j_error = np.linalg.norm(J @ J + np.eye(4))
    commutator = np.linalg.norm(Hrel @ J - J @ Hrel)
    assert j_error < 1e-12 and commutator < 1e-5

    # Pairing holds for every Z3 x U(1)-invariant extension sampled.
    pair_gaps = []
    for _ in range(12):
        f = random_invariant_potential(rng)
        Hf, _, res_f = stationary_hessian(f, phi)
        assert res_f < 1e-5
        ev = np.linalg.eigvalsh(Q.T @ Hf @ Q)
        pair_gaps.append(max(ev[1] - ev[0], ev[3] - ev[2]))
    assert max(pair_gaps) < 1e-4

    # The most symmetric Z3 x U(1) extension: V(rho_k) only, U(1)^3.
    def diagonal(x):
        rk = abs(to_complex(x))**2
        return 0.5 * g0 * rk.sum()**2 + 0.4 * np.sum(rk**2)
    diag_dim, _ = symmetry_algebra_dim(diagonal)
    Hd, _, _ = stationary_hessian(diagonal, phi)
    diag_null = int(np.sum(abs(np.linalg.eigvalsh(Q.T @ Hd @ Q)) < 1e-5))
    assert diag_dim == 3 and diag_null == 2

    locked_rank = int(np.sum(abs(rel) > 1e-5))
    rulers_committed, three_rulers_possible = 4 - locked_rank, locked_rank <= 1
    assert locked_rank == 4 and rulers_committed == 0 and not three_rulers_possible

    report(f"Force check |dV - M psi|={force_error:.1e}. Mu={mu:.9f}.")
    report(f"Symmetry algebra of V in so(6): dim {sym_dim} (common phase).")
    report(f"Clock Hessian (6): {np.round(full, 6)}")
    report(f"Fixed rho_tot (5): {np.round(sphere, 6)}. One zero, the clock.")
    report(f"Relative (4): {np.round(rel, 6)} = (5 -/+ sqrt13) g1 rho/2, "
           "paired, none zero.")
    report(f"Z3' structure: |J^2+1|={j_error:.1e}, |[H,J]|={commutator:.1e}. "
           f"Random extension pair gaps <= {max(pair_gaps):.1e}.")
    report(f"U(1)^3 extension: symmetry dim {diag_dim}, relative zero modes "
           f"{diag_null}: at most one clock and two rulers.")
    report("RESULT: one clock, zero rulers. Three rulers need the whole "
           "relative block to vanish (unlocked phases) and exceed the rank "
           "bound 3 in any case: strictly incompatible with the locked "
           "circulant.")
    return dict(mu=mu, force_error=force_error, symmetry_dim=sym_dim,
                clock_hessian=full, sphere_hessian=sphere,
                relative_hessian=rel, relative_expected=expected,
                j_error=j_error, commutator=commutator,
                pair_gaps=pair_gaps, diagonal_symmetry_dim=diag_dim,
                diagonal_relative_null=diag_null, locked_rank=locked_rank,
                rulers=rulers_committed,
                three_rulers_possible=three_rulers_possible)


if __name__ == "__main__":
    run()
