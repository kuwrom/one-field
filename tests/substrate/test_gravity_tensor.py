"""Independent checks of the clock operator and the tensor projection."""

import math

import numpy as np
import pytest

from probes import gravity_tensor as gt
from nls_soliton import _bdg_analysis, interaction_matrix


@pytest.mark.parametrize("g0,g1,rho", [(1, 1/math.sqrt(2), 1), (2, .3, .7)])
def test_clock_hessian_matches_live_force_and_symplectic_dynamics(g0, g1, rho):
    phi, mu, r, A, B = gt.clock_background(g0, g1, rho)
    assert np.linalg.norm(interaction_matrix(phi, g0, g1) @ phi - mu*phi) < 1e-12
    assert A == pytest.approx(A.conj().T, abs=1e-12)
    assert B == pytest.approx(B.T, abs=1e-12)
    for step in (2e-4, 1e-5):
        An, Bn = gt.live_force_jacobian(phi, g0, g1, mu, step)
        assert A == pytest.approx(An, abs=2e-7)
        assert B == pytest.approx(Bn, abs=2e-7)
    # Canonical Hessian has one global phase flat direction, not a gauge field.
    H = gt.real_hessian(A, B)
    J = np.block([[np.zeros((3, 3)), np.eye(3)],
                  [-np.eye(3), np.zeros((3, 3))]]) / 2
    assert np.linalg.matrix_rank(J) == 6
    phase = np.r_[-phi.imag, phi.real]
    assert np.linalg.norm(H @ phase) < 1e-12
    p = .37
    generator = J @ (H + p*p*np.eye(6))
    expected = np.sort(np.r_[-gt.clock_frequencies(p, mu, r),
                              gt.clock_frequencies(p, mu, r)])
    assert np.sort(np.linalg.eigvals(generator).imag) == pytest.approx(expected, abs=1e-11)


def test_hessian_recovers_symmetric_background_and_generic_derivatives():
    g0, g1, rho = 1.0, .41, 1.2
    mu = (3*g0 + 2*g1)*rho
    A, B = gt.interaction_hessian(np.full(3, math.sqrt(rho)), g0, g1, mu)
    legacy = _bdg_analysis(g0, g1, rho)
    assert A == pytest.approx(legacy["A"], abs=1e-12)
    # B can be independently reconstructed from the published symmetric formula.
    assert B == pytest.approx(rho*(np.full((3, 3), g0+g1)-g1*np.eye(3)), abs=1e-12)
    phi = np.array([.3+.2j, -.1+.7j, 1.1-.8j])
    A, B = gt.interaction_hessian(phi, g0, g1, mu)
    An, Bn = gt.live_force_jacobian(phi, g0, g1, mu)
    assert A == pytest.approx(An, abs=1e-8)
    assert B == pytest.approx(Bn, abs=1e-8)


def test_clock_aligned_fourier_operator_and_raw_phase_label():
    phi, mu, r, A, B = gt.clock_background()
    D = np.diag(phi)
    U = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/math.sqrt(3)
    transform = D @ U
    assert transform.conj().T @ A @ transform == pytest.approx(np.diag([mu, 2*r, r/2]), abs=1e-12)
    assert transform.conj().T @ B @ transform.conj() == pytest.approx(
        np.array([[mu, 0, 0], [0, 0, r/2], [0, r/2, 0]]), abs=1e-12)
    assert abs(U.conj().T @ (1j*phi)) == pytest.approx([0, math.sqrt(3), 0], abs=1e-12)


def test_tt_projection_is_covariant_and_tensor_starts_at_second_order():
    rng = np.random.default_rng(20260928)
    gradient = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    Q = gt.neutral_tensor(gradient)
    k = np.array([.2, -.4, .7])
    T = gt.tt_project(Q, k)
    assert abs(np.trace(T)) < 1e-12
    assert np.linalg.norm(T @ k) < 1e-12
    assert gt.tt_project(T, k) == pytest.approx(T, abs=1e-12)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    assert gt.tt_project(rotation @ Q @ rotation.T, rotation @ k) == pytest.approx(
        rotation @ T @ rotation.T, abs=1e-12)
    assert gt.neutral_tensor(np.roll(gradient, 1, axis=0)) == pytest.approx(Q, abs=1e-12)
    assert gt.neutral_tensor(gradient*np.exp(.4j)) == pytest.approx(Q, abs=1e-12)
    for amplitude in (.1, .01):
        assert gt.neutral_tensor(amplitude*gradient) == pytest.approx(amplitude**2*Q, abs=1e-12)
    # Collinear pairs carry no TT vertex even at the continuum threshold.
    assert gt.tt_project(np.outer(k, k), k) == pytest.approx(np.zeros((3, 3)), abs=1e-12)
