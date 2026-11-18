"""Clock vacuum fluctuation operator and neutral tensor response.

This is a Gaussian calculation for the *implemented* three complex scalar
NLS, not a tensor completion of it. The Hessian is differentiated from
the interaction in nls_soliton.py and checked against its live force.

Write psi = exp(-i mu t) (phi + eta), phi_a = sqrt(rho) w^a,
w = exp(2 pi i/3), r = g1 rho, mu = (3 g0 - g1) rho. Then

    i d_t (eta, eta*) = [[A+epsilon, B], [-B*, -A*-epsilon]] (eta, eta*)

with epsilon = |p|^2/2. In the condensate aligned frame eta_a/w^a,
the unitary Z3 transform gives A = diag(mu, 2r, r/2), B_00 = mu,
B_12 = B_21 = r/2. The three positive frequencies are

    omega_0 = sqrt(epsilon (epsilon + 2 mu)),
    omega_+/- = sqrt((epsilon + 5r/4)^2 - r^2/4) +/- 3r/4.

The q=0 label here is in the aligned frame. The physical phase tangent
i phi lies in raw component Fourier q=1 (q=2 for the conjugate clock).
The first order NLS has three canonical pairs, a nondegenerate
symplectic form, and a *global* U(1), not a local metric gauge freedom.

Q_ij = Re sum_a d_i psi_a* d_j psi_a - delta_ij |grad psi|^2/3
is internally neutral and spatially rank two. On a uniform background
it starts at order eta^2. There is no one quasiparticle Q vertex.
For standard canonical bosonic quantization of the quadratic NLS, its
zero temperature positive frequency spectral measure is therefore a
two quasiparticle continuum. This quantization is an explicit premise.
The classical substrate simulations alone do not establish it.

The spectral density below sums the two orthonormal TT polarizations,
per unit volume, with convention integral d^3p/(2pi)^3 |vertex|^2
delta(Omega E_p E_{k-p}) (no extra 2pi in the frequency delta).
Neutral pairs are 0,0 (identical: factor 1/2) and +,- (distinct: 1).
Their coherence vertex is u_p v_q + v_p u_q. Angular integration of
the frequency delta reduces the answer to one radial quadrature.
An externally imposed TT source is only a diagnostic, not a new field.

This excludes a separate tensor pole at Gaussian order about this
uniform vacuum. It does not exclude an interacting bound state or a
different phase. No bridge mass or gravitational coupling is inferred.

RESULT. Computed about the uniform clock vacuum at g0 > 0, the gas sign:
the analytic Hessian matches the live force, with one global U(1) zero
direction and no local metric gauge constraint. The neutral tensor Q
starts at order eta^2, and at k = 0 its infrared TT density scales as
Omega^4, a two quasiparticle continuum. Under standard canonical bosonic
quantization a separate massless spin 2 pole at Gaussian order about
this vacuum is excluded.

Usage: python3 tests/probes/gravity_tensor.py
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import G1_OVER_G0, interaction_matrix  # noqa: E402


def interaction_hessian(phi, g0, g1, mu):
    """Analytic Wirtinger A=dF/dpsi, B=dF/dpsi* for alpha=0.

    F = interaction_matrix(psi) @ psi - mu psi. This function also
    accepts nonstationary backgrounds for independent derivative checks.
    """
    phi = np.asarray(phi, dtype=complex)
    if phi.shape != (3,):
        raise ValueError("expected three complex background components")
    A = g0 * (np.vdot(phi, phi).real * np.eye(3)
              + np.outer(phi, phi.conj())) - mu * np.eye(3)
    B = g0 * np.outer(phi, phi)
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        A[a, a] += g1 * np.real(phi[b].conjugate() * phi[c])
        A[a, b] += g1 / 2 * (phi[a] * phi[c].conjugate()
                            + abs(phi[c])**2 + phi[b].conjugate() * phi[c])
        A[a, c] += g1 / 2 * (phi[a] * phi[b].conjugate()
                            + abs(phi[b])**2 + phi[c].conjugate() * phi[b])
        B[a, b] += g1 / 2 * (phi[a] * phi[c] + phi[b] * phi[c])
        B[a, c] += g1 / 2 * (phi[a] * phi[b] + phi[c] * phi[b])
    return A, B


def live_force_jacobian(phi, g0, g1, mu, step=1e-5):
    """Independent central differences of the actual nonlinear force."""
    def force(v):
        return interaction_matrix(v, g0, g1) @ v - mu * v

    dx = np.column_stack([(force(phi + step * e) - force(phi - step * e))
                          / (2 * step) for e in np.eye(3)])
    dy = np.column_stack([(force(phi + 1j * step * e)
                           - force(phi - 1j * step * e))
                          / (2 * step) for e in np.eye(3)])
    return (dx - 1j * dy) / 2, (dx + 1j * dy) / 2


def clock_background(g0=1.0, g1=G1_OVER_G0, rho=1.0):
    r, mu = g1 * rho, (3 * g0 - g1) * rho
    if rho <= 0 or r <= 0 or mu <= 0:
        raise ValueError("clock calculation requires rho, g1*rho, mu > 0")
    phi = math.sqrt(rho) * np.exp(2j * np.pi * np.arange(3) / 3)
    A, B = interaction_hessian(phi, g0, g1, mu)
    return phi, mu, r, A, B


def bdg_matrix(A, B, momentum):
    """Six component Nambu operator. Momentum is a spatial 3-vector."""
    momentum = np.asarray(momentum, dtype=float)
    if momentum.shape != (3,):
        raise ValueError("momentum must be a spatial 3-vector")
    Ap = A + np.dot(momentum, momentum) / 2 * np.eye(3)
    return np.block([[Ap, B], [-B.conj(), -Ap.conj()]])


def real_hessian(A, B):
    """Hessian of E-mu N in (Re eta, Im eta). Generator is J H / 2."""
    return 2 * np.block([[(A + B).real, -(A - B).imag],
                         [(A + B).imag, (A - B).real]])


def clock_frequencies(p, mu, r):
    """Ordered by channel: phase, lower relative, upper relative."""
    epsilon = np.asarray(p, dtype=float)**2 / 2
    relative = np.sqrt((epsilon + 5 * r / 4)**2 - r**2 / 4)
    return np.array([np.sqrt(epsilon * (epsilon + 2 * mu)),
                     relative - 3 * r / 4, relative + 3 * r / 4])


def neutral_tensor(gradient):
    """Q for gradient[a,i]=d_i psi_a, with three internal components."""
    gradient = np.asarray(gradient, dtype=complex)
    if gradient.shape != (3, 3):
        raise ValueError("gradient must have internal, spatial shape (3, 3)")
    Q = (gradient.conj().T @ gradient).real
    return Q - np.eye(3) * np.trace(Q) / 3


def tt_project(tensor, momentum):
    """Transverse, traceless part of a symmetric spatial tensor."""
    momentum = np.asarray(momentum, dtype=float)
    norm = np.linalg.norm(momentum)
    if momentum.shape != (3,) or norm == 0:
        raise ValueError("TT projection requires a nonzero spatial 3-vector")
    P = np.eye(3) - np.outer(momentum, momentum) / norm**2
    tensor = np.asarray(tensor)
    symmetric = (tensor + tensor.T) / 2
    projected = P @ symmetric @ P
    return projected - P * np.trace(projected) / 2


def _pair_parameters(channel, mu, r):
    if channel == "phase":
        return mu, mu, 0.5
    if channel == "relative":
        return 5 * r / 4, r / 2, 1.0
    raise ValueError("channel must be phase or relative")


def pair_threshold(k, mu, r, channel="phase"):
    """Infimum of continuum support; TT vertex vanishes at the edge."""
    a, b, _ = _pair_parameters(channel, mu, r)
    return 2 * math.sqrt((k**2 / 8 + a)**2 - b**2)


def tensor_spectral_density(k, frequency, mu, r, channel="phase"):
    """Gaussian Q^TT spectral density for one neutral pair channel.

    Frequencies of the relative pair are D(p)+s and D(q)-s, so their
    splitting cancels. Both channels reduce to E(p)=sqrt((p^2/2+a)^2-b^2).
    No Lorentzian peak broadening, fitted pole, or adjustable mass is used.
    k=0 denotes the limit with a fixed choice of TT projection axis.
    """
    if k < 0 or frequency < 0 or mu <= 0 or r <= 0:
        raise ValueError("require k, frequency >= 0 and mu, r > 0")
    a, b, symmetry = _pair_parameters(channel, mu, r)

    def energy(p):
        epsilon = p * p / 2
        return math.sqrt(epsilon * (epsilon + 2 * a) + (a - b) * (a + b))

    def inverse_energy(e):
        # Rationalization avoids cancellation near the phase mode origin.
        epsilon = (e * e - (a - b) * (a + b)) / (math.hypot(e, b) + a)
        return math.sqrt(max(0.0, 2 * epsilon))

    if frequency <= 2 * energy(k / 2):
        return 0.0
    if k == 0:
        p = inverse_energy(frequency / 2)
        e = frequency / 2
        return symmetry * p**5 * b**2 / (15 * math.pi**2 * e * (p*p/2 + a))

    # Exact collinear boundaries of the allowed radial interval.
    if frequency < energy(0) + energy(k):
        lo = brentq(lambda p: energy(p) + energy(k - p) - frequency,
                    0, k / 2)
        hi = k - lo
    else:
        lo = brentq(lambda p: energy(p) + energy(k + p) - frequency,
                    0, inverse_energy(frequency - energy(0)))
        hi = lo + k

    def coherence(p, e):
        u = math.sqrt((p*p/2 + a + e) / (2 * e))
        v = -b / math.sqrt(2 * e * (p*p/2 + a + e))
        return u, v

    def integrand(p):
        ep = energy(p)
        eq = frequency - ep
        q = inverse_energy(eq)
        if p == 0 or q == 0:
            return 0.0
        cosine = np.clip((k*k + p*p - q*q) / (2 * k * p), -1, 1)
        vertex_sq = p**4 * (1 - cosine*cosine)**2 / 2
        up, vp = coherence(p, ep)
        uq, vq = coherence(q, eq)
        # Angular delta Jacobian: q/(k p E'(q)) = eq/(k p (q^2/2+a)).
        measure = p * eq / (4 * math.pi**2 * k * (q*q/2 + a))
        return symmetry * measure * vertex_sq * (up*vq + vp*uq)**2

    return quad(integrand, lo, hi, epsabs=1e-12, epsrel=2e-8)[0]


def run(report=print):
    phi, mu, r, A, B = clock_background()
    An, Bn = live_force_jacobian(phi, 1.0, G1_OVER_G0, mu)
    jacobian_error = max(np.max(abs(A - An)), np.max(abs(B - Bn)))
    assert jacobian_error < 1e-8
    H = real_hessian(A, B)
    hessian_eigenvalues = np.linalg.eigvalsh(H)
    phase_ward_error = np.linalg.norm(A @ phi - B @ phi.conj())
    assert phase_ward_error < 1e-12
    assert np.count_nonzero(abs(hessian_eigenvalues) < 1e-10) == 1
    assert min(hessian_eigenvalues) > -1e-10

    # Matrix spectrum is checked independently against the analytic branches.
    spectrum_error = 0.0
    for p in (0.01, 0.1, 0.6, 1.7, 4.0):
        expected = np.sort(np.r_[-clock_frequencies(p, mu, r),
                                  clock_frequencies(p, mu, r)])
        for direction in ((0, 0, 1), (1, 2, -2)):
            direction = np.array(direction) / np.linalg.norm(direction)
            eig = np.linalg.eigvals(bdg_matrix(A, B, p * direction))
            assert max(abs(eig.imag)) < 1e-10
            spectrum_error = max(spectrum_error, max(abs(np.sort(eig.real) - expected)))
    assert spectrum_error < 1e-10
    gaps = clock_frequencies(0, mu, r)
    report(f"Clock: mu={mu:.9f}. Gaps={gaps}; 3 canonical scalar branches.")
    report(f"Live force Jacobian error={jacobian_error:.2e}; "
           f"matrix/analytic spectrum error={spectrum_error:.2e}.")
    report("One global U(1) zero direction. No local metric gauge constraint.")

    gradient = np.array([[1+2j, 2-1j, 0.3], [0.2j, 1, -0.7j], [1, 2j, 3]])
    Q = neutral_tensor(gradient)
    assert np.allclose(neutral_tensor(np.roll(gradient, 1, axis=0)), Q)
    assert np.allclose(neutral_tensor(0.01 * gradient), 0.01**2 * Q)
    qt = tt_project(Q, (0, 0, 1))
    assert np.linalg.norm(qt) > 0
    assert abs(np.trace(qt)) < 1e-12 and np.linalg.norm(qt[:, 2]) < 1e-12
    report("Neutral Q has a nonzero TT projection, but its uniform vacuum "
           "vertex starts quadratically: no one quasiparticle residue.")

    k = 1.2
    center = float(clock_frequencies(k, mu, r)[0])
    thresholds = {c: pair_threshold(k, mu, r, c) for c in ("phase", "relative")}
    windows = (0.02, 0.01, 0.005)
    weights = [quad(lambda f: tensor_spectral_density(k, f, mu, r),
                    center - width/2, center + width/2, epsabs=1e-12)[0]
               for width in windows]
    densities = np.array(weights) / windows
    assert all(w > 0 for w in weights)
    assert max(densities) / min(densities) < 1.02
    ir_frequencies = np.array([0.02, 0.01, 0.005]) * mu
    ir_density = np.array([tensor_spectral_density(0, f, mu, r)
                           for f in ir_frequencies])
    ir_power = math.log(ir_density[-2] / ir_density[-1]) / math.log(2)
    assert abs(ir_power - 4) < 0.001
    report(f"At |k|={k}: continuum edges={thresholds}; "
           f"phase frequency={center:.9f}.")
    report(f"Shrinking window weights about that frequency: {weights}; "
           "weight/width stays finite (continuum, not a delta pole).")
    report(f"At k=0: infrared TT spectral density scales as Omega^{ir_power:.6f}.")
    report("RESULT: no separate massless spin-2 pole at Gaussian order. "
           "Interacting composites and other backgrounds remain untested.")
    return dict(mu=mu, r=r, gaps=gaps, jacobian_error=jacobian_error,
                spectrum_error=spectrum_error, phase_ward_error=phase_ward_error,
                hessian_eigenvalues=hessian_eigenvalues, thresholds=thresholds,
                window_widths=windows, window_weights=weights,
                infrared_power=ir_power)


if __name__ == "__main__":
    run()
