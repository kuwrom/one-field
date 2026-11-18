"""Lambda route A: a volume element built from the carrier.

THE CONSTRUCTION. The carrier's common phase current J = Im sum_a
conj(phi_a) d phi_a gives the exact four form Phi = dJ ^ dJ, a density
of weight one built from the complex structure alone. A Weyl invariant
metric g^ = h (Phi / sqrt det h)^(1/2) has sqrt det g^ = |Phi|, and any
constant times Phi is a boundary term: the abelian case of the
Pontryagin density unimodular gravity (Jirousek, arXiv:2301.01662
Eqs. 50-58, "in the Abelian case P represents the Pfaffian of the
matrix F"), and the measure of two measure theories (Guendelman and
Kaganovich, gr-qc/9905029 Eq. 4).

ONE BILINEAR. With H_mu nu = sum_a conj(d_mu phi_a) d_nu phi_a, a
Hermitian 4x4 matrix, the composite metric is h = 2 Re H and the
current's curl is dJ = 2 Im H (antisymmetric). Phi = Pf(dJ) d^4x. The
metric and the would be volume element are the real and imaginary
parts of one object.

STATEMENTS (the numbers below check each one).

  A1 (the form). Phi transforms as a density of weight one under a
     linear change of coordinates, and vanishes on every configuration
     whose components are proportional to one complex function, the
     clock lump included.
  A2 (the two measure condition). The measure fields must not appear in
     the Lagrangian coupled to Phi (gr-qc/9905029: "We assume in what
     follows that L does not contain explicitly the measure fields").
     Here the same fields build h and U, so the shift symmetry that
     makes L = M a first integral is absent.  Stated from the source,
     not computed; A3 and A4 close the route before it matters.
  A3 (classical: no Lorentzian vacuum). The image of Im H lies in the
     span of {Re d phi_a, Im d phi_a}, and Re H is positive with that
     span as its image.  Phi != 0 needs Im H of rank four, so the span
     is all of R^4 and h is positive definite.  Wherever the classical
     volume element is nonzero the metric is Euclidean, and g^, being
     conformal to h, is too.  Checked on random configurations, and the
     converse: rank h < 4 gives Phi = 0.
  A4 (quantum: <Phi> = 0). In a Gaussian state with
     <conj(d_mu phi_a) d_nu phi_b> = delta_ab C_mu nu, Wick's theorem
     gives <Phi> = sum of eps contractions of C x C.  For C real
     symmetric, the Lorentz invariant vacuum's C proportional to eta
     included, every term vanishes: <Phi> = 0 exactly, with no parity
     argument needed.  A nonzero <Phi> needs Im C of rank four, an
     antisymmetric correlator that fixes a frame.

RESULT. Route A is excluded. The carrier's own volume element is nonzero
only where its metric is Euclidean, and its expectation vanishes in
every vacuum with a symmetric gradient correlator. The gauge field
version (Pontryagin density of an independent Yang-Mills field) needs
<F F~> != 0, a CP-odd condensate the vacuum of theta bar = 0 does not
carry (interference/couplings.py, strong CP). Stated, not computed.

Usage: python3 tests/probes/composite_volume.py     (a second)
"""

import itertools

import numpy as np

EPS = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS[_p] = np.linalg.det(np.eye(4)[list(_p)])


def bilinear(grad):
    """grad[a, mu] = d_mu phi_a (complex): H = sum_a conj(d phi_a) (x) d phi_a."""
    return np.einsum("am,an->mn", grad.conj(), grad)


def phi_density(grad):
    """Phi = Pf(dJ) with dJ = 2 Im H: (1/8) eps F F."""
    F = 2 * np.imag(bilinear(grad))
    return 0.125 * np.einsum("abcd,ab,cd->", EPS, F, F)


def a1(rng):
    grad = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
    L = rng.normal(size=(4, 4))                       # x -> L x: d_mu -> (L^-1)^T d
    grad2 = grad @ np.linalg.inv(L)
    weight = phi_density(grad2) / phi_density(grad) * np.linalg.det(L)
    f = rng.normal(size=4) + 1j * rng.normal(size=4)  # all components prop. to one function
    clock = np.outer(np.exp(2j * np.pi * np.arange(3) / 3), f)
    return dict(weight_one=float(weight), clock_phi=float(abs(phi_density(clock))))


def a3(rng, n=2000):
    worst_min_eig, rank_def = np.inf, 0.0
    for _ in range(n):
        grad = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
        if abs(phi_density(grad)) > 1e-6:
            h = 2 * np.real(bilinear(grad))
            worst_min_eig = min(worst_min_eig, np.linalg.eigvalsh(h).min())
    for _ in range(200):                              # gradients spanning only 3 directions
        basis = rng.normal(size=(3, 4))
        grad = (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) @ basis
        rank_def = max(rank_def, abs(phi_density(grad)))
    return dict(min_eig_when_phi=float(worst_min_eig), phi_when_rank3=float(rank_def))


def wick_phi(C, n_comp=3):
    """<Phi> for <conj(d_mu phi_a) d_nu phi_b> = delta_ab C_mu nu (complex Gaussian)."""
    # Phi = (1/8) eps F F, F = -i (H - H^T); <H_mn H_rs> = sum_ab [d_aa d_bb C_mn C_rs + d_ab d_ab C_ms C_rn]
    I = np.eye(n_comp)
    n_same = np.einsum("aa,bb->", I, I)               # n^2
    n_cross = np.einsum("ab,ab->", I, I)              # n
    HH = n_same * np.einsum("mn,rs->mnrs", C, C) + n_cross * np.einsum("ms,rn->mnrs", C, C)
    # <F_mn F_rs> = -<(H - H^T)_mn (H - H^T)_rs>
    FF = -(HH - HH.transpose(1, 0, 2, 3) - HH.transpose(0, 1, 3, 2) + HH.transpose(1, 0, 3, 2))
    return complex(0.125 * np.einsum("mnrs,mnrs->", EPS, FF))


def a4(rng):
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    S = rng.normal(size=(4, 4))
    real_sym = S @ S.T
    A = rng.normal(size=(4, 4))
    herm = real_sym + 1j * (A - A.T)                  # Im part antisymmetric, rank four
    return dict(lorentz=abs(wick_phi(0.3 * eta)), real_symmetric=abs(wick_phi(real_sym)),
                frame_fixing=abs(wick_phi(herm)))


def run(report=print, seed=11):
    rng = np.random.default_rng(seed)
    r1, r3, r4 = a1(rng), a3(rng), a4(rng)
    out = dict(a1=r1, a3=r3, a4=r4)
    assert abs(r1["weight_one"] - 1) < 1e-10 and r1["clock_phi"] < 1e-12
    assert r3["min_eig_when_phi"] > 0 and r3["phi_when_rank3"] < 1e-10
    assert r4["lorentz"] < 1e-12 and r4["real_symmetric"] < 1e-10 and r4["frame_fixing"] > 1e-3
    report(f"A1 Phi is a density of weight one (Phi' det L / Phi = {r1['weight_one']:.12f}); "
           f"on a clock configuration |Phi| = {r1['clock_phi']:.1e}.")
    report("A2 the two measure condition (measure fields absent from L) fails by construction: "
           "the same fields build h and U.  Stated from the source.")
    report(f"A3 wherever Phi != 0 the smallest eigenvalue of h is {r3['min_eig_when_phi']:.3e} > 0: "
           f"h is Euclidean there. With gradients in three directions |Phi| <= "
           f"{r3['phi_when_rank3']:.1e}.")
    report(f"A4 Wick <Phi>: Lorentz invariant vacuum {r4['lorentz']:.1e}, any real symmetric "
           f"correlator {r4['real_symmetric']:.1e}, a frame fixing Hermitian correlator "
           f"{r4['frame_fixing']:.3f}.")
    report("RESULT: route A is excluded. The carrier's volume element is nonzero only where its "
           "metric is Euclidean, and vanishes in expectation in every vacuum with a symmetric "
           "gradient correlator.")
    return out


if __name__ == "__main__":
    run()
