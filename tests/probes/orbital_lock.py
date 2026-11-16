"""Orbital lock: can a component direction be tied to a spatial one?

PREMISES. The committed substrate (tests/substrate/nls_soliton.py):
three complex components psi_a, spatial scalars, energy

    E = int { 1/2 sum_a |grad psi_a|^2 + (g0/2) rho_tot^2
              + g1 sum_k rho_k Re(conj(psi_{k+1}) psi_{k+2}) },

g1 = g0/sqrt(2) > 0, alpha = 0. Exact symmetry: G_int x SO(3) x
translations x Galilean boosts, with G_int = S3 x U(1) (cyclic Z3 on
the index, the reflection a -> -a at alpha = 0) and time reversal.
The only continuous internal factor is the global U(1).

STATEMENTS (proof sketches; the numbers below check them):

  L1 (Schur). SO(3) acts on no internal index, so every invariant
     coefficient tensor C^A_ij equals c^A delta_ij: no local term, and
     no correlation function of a rotation symmetric state, ties a
     component index to a spatial one.
  T2 (no isotropic lock). SO(3) is simple and not virtually abelian,
     G_int is virtually abelian. By Goursat a stabilizer projecting
     onto SO(3) contains {1} x SO(3). A locked state has a proper
     closed spatial residual H: axial (inside O(2)) or T, O.
  T3 (locks break translations). Invariance under translations modulo
     G_int forces a common plane wave (R^3 connected, U(1) the only
     continuous factor): uniform in its rest frame, hence unlocked.
  T4 (global minimum). Wolstenholme: g1-term >= -(g1/2) e2(rho) >=
     -(g1/6) rho_tot^2 pointwise, so E - E_clock >= 1/2 int |grad psi|^2
     at fixed N. Every locked state costs at least its gradient energy.
  B  (boost). In the rest frame sum_a rho_a u_a = 0 with rho_a > 0, so
     T_ij = sum_a rho_a u_a^i u_a^j has rank <= 2 for ANY densities.

RESULT (numbers from run()):

  1. Bilinear coefficient tensors t_{ab,ij}: Z3 x SO(3) has 3
     invariants, all spatially delta. Twisted C3 (index shift with a
     120 degree turn about (1,1,1)): 27, neutral contractions
     anisotropic at O(1). Twisted T (add the axis half turns): 9,
     locked per component but every neutral contraction is delta.
  2. Rest frame tensor s_min/s_max <= 2e-16. The Z3 counterflow gives
     T = (3/2) rho Q^2 (1 - n n), uniaxial.  3. T4 bound holds.
  4. Counterflow q = (Q,0), (0,Q), (-Q,-Q), Q = 2 pi/L: ansatz excess
     1/2 rho sum q^2 + (3/2) g1 rho^2 to 1e-15. Residual after
     per component mu 0.707, 0.354 at g1/2 (not stationary). Gradient
     flow changes the windings through nodes, relaxed excess 1.098
     (L = 8), 0.340 (L = 16): a wall remnant, no bulk locked phase.
  5. Explicit triad locks about the clock, kappa = 0.3, p = 0.1:
     kinetic, generic C3 triad: phonon 17%, relative branch 37%
     anisotropy, p-independent. Density gradient (boost invariant),
     generic: phonon anisotropy ~ p^2 (rho_s = rho by boosts), relative
     branch 43%. Orthonormal triad (residual T): phonon ~ p^2
     (kinetic) or p^4 (density): anisotropy only at rank four.
  6. p-wave pair operator Phi_c = eps_abc psi_a grad psi_b: boost
     invariant to 1e-15, zero on a common twist, and sum |Phi_c|^2 =
     rho^2 sum_{a<b} |q_a - q_b|^2 on a counterflow. Lambda sum |Phi|^2
     makes the relative stiffness rho/2 + 3 lambda rho^2, soft at
     lambda = -1/(6 rho), and its minimum is again a coplanar lock.

No parameter is fitted. kappa = 0.3 and lambda are probes of a
hypothetical term, not substrate values.

Usage: python3 tests/probes/orbital_lock.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import G1_OVER_G0, interaction_matrix  # noqa: E402
from gravity_tensor import clock_background  # noqa: E402

G0, G1, RHO = 1.0, G1_OVER_G0, 1.0
W = np.exp(2j * np.pi / 3)
SHIFT = np.roll(np.eye(3), 1, axis=0)          # index shift a -> a+1
AXIS_CYCLE = np.roll(np.eye(3), 1, axis=0)     # x -> y -> z, 120 deg


def _rotation(axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * K @ K


def invariant_subspace(generators):
    """Common fixed vectors of 81 x 81 matrices acting on t_{ab,ij}."""
    stack = np.vstack([g - np.eye(81) for g in generators])
    _, s, vt = np.linalg.svd(stack)
    return vt[np.sum(s > 1e-9):]


def _act(internal, spatial):
    return np.kron(np.kron(internal, internal), np.kron(spatial, spatial))


def schur_check():
    """Invariants of bilinear coefficient tensors under three groups."""
    I3 = np.eye(3)
    R1, R2 = _rotation((1, 2, 3), 0.7), _rotation((-2, 1, 0.5), 1.3)
    groups = {
        "Z3 x SO(3)": [_act(SHIFT, I3), _act(I3, R1), _act(I3, R2)],
        "twisted C3": [_act(SHIFT, AXIS_CYCLE)],
        "twisted T": [_act(SHIFT, AXIS_CYCLE),
                      _act(I3, np.diag([1., -1, -1])),
                      _act(I3, np.diag([-1., 1, -1]))],
    }
    circulants = [I3, SHIFT, SHIFT.T]
    out = {}
    for name, gens in groups.items():
        basis = invariant_subspace(gens)
        neutral, component = 0.0, 0.0
        for v in basis:
            t = v.reshape(3, 3, 3, 3)
            for c in circulants:
                s = np.einsum("ab,abij->ij", c, t)
                s = (s + s.T) / 2
                neutral = max(neutral, np.linalg.norm(
                    s - np.eye(3) * np.trace(s) / 3))
            for a in range(3):
                s = (t[a, a] + t[a, a].T) / 2
                component = max(component, np.linalg.norm(
                    s - np.eye(3) * np.trace(s) / 3))
        out[name] = dict(dim=len(basis), neutral_aniso=neutral,
                         component_aniso=component)
    return out


def rest_frame_rank(trials=2000, seed=0):
    rng = np.random.default_rng(seed)
    worst_rest, lab_ratios = 0.0, []
    for _ in range(trials):
        rho = rng.uniform(0.1, 2.0, 3)
        q = rng.standard_normal((3, 3))
        lab = np.einsum("a,ai,aj->ij", rho, q, q)
        u = q - rho @ q / rho.sum()
        rest = np.einsum("a,ai,aj->ij", rho, u, u)
        s_rest = np.linalg.svd(rest, compute_uv=False)
        s_lab = np.linalg.svd(lab, compute_uv=False)
        worst_rest = max(worst_rest, s_rest[-1] / s_rest[0])
        lab_ratios.append(s_lab[-1] / s_lab[0])
    n = np.ones(3) / math.sqrt(3)
    u0 = np.array([1., -1, 0]) / math.sqrt(2)
    Q = 0.8
    qz3 = np.array([Q * _rotation(n, 2 * math.pi * a / 3) @ u0
                    for a in range(3)])
    Tz3 = RHO * qz3.T @ qz3
    z3_error = np.linalg.norm(Tz3 - 1.5 * RHO * Q**2 * (np.eye(3)
                                                       - np.outer(n, n)))
    return dict(rest_min_over_max=worst_rest, lab_median_ratio=float(np.median(lab_ratios)),
                z3_uniaxial_error=z3_error)


def _g1_density(psi):
    rho = np.abs(psi)**2
    return sum(rho[k] * np.real(np.conj(psi[(k + 1) % 3]) * psi[(k + 2) % 3])
               for k in range(3))


def pointwise_bound(samples=200000, seed=1):
    rng = np.random.default_rng(seed)
    psi = (rng.standard_normal((3, samples))
           + 1j * rng.standard_normal((3, samples)))
    rho = np.abs(psi)**2
    s = _g1_density(psi)
    e2 = rho[0] * rho[1] + rho[1] * rho[2] + rho[2] * rho[0]
    clock = np.sqrt(1 / 3) * W**np.arange(3)
    return dict(min_ratio=float(np.min(s / rho.sum(0)**2)),
                wolstenholme_slack=float(np.min(s + e2 / 2)),
                clock_ratio=float(_g1_density(clock[:, None])[0]))


class Box:
    def __init__(self, L, N):
        self.L, self.N, self.dx = L, N, L / N
        x = np.arange(N) * self.dx
        self.X, self.Y = np.meshgrid(x, x, indexing="ij")
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.KX, self.KY = np.meshgrid(k, k, indexing="ij")
        self.K2 = self.KX**2 + self.KY**2

    def grad(self, f):
        F = np.fft.fft2(f, axes=(-2, -1))
        return (np.fft.ifft2(1j * self.KX * F, axes=(-2, -1)),
                np.fft.ifft2(1j * self.KY * F, axes=(-2, -1)))

    def energy(self, psi, g1=G1):
        F = np.fft.fft2(psi, axes=(1, 2))
        kin = 0.5 * np.sum(self.K2 * np.abs(F)**2) / self.N**2
        rt = np.sum(np.abs(psi)**2, axis=0)
        return float(kin + np.sum(0.5 * G0 * rt**2 + g1 * _g1_density(psi))
                     ) * self.dx**2

    def h_psi(self, psi, g1=G1):
        lap = np.fft.ifft2(-self.K2 * np.fft.fft2(psi, axes=(1, 2)),
                           axes=(1, 2))
        M = interaction_matrix(psi, G0, g1)
        return -0.5 * lap + np.einsum("...ij,j...->i...", M, psi)

    def winding(self, psi):
        def turn(z):
            return round(float(np.angle(np.roll(z, -1) / z).sum())
                         / (2 * math.pi))
        return [(turn(p[:, 0]), turn(p[0, :])) for p in psi]


def counterflow(L, N, steps, seed=0):
    box = Box(L, N)
    Q = 2 * math.pi / L
    qs = np.array([(Q, 0), (0, Q), (-Q, -Q)])
    psi = np.array([math.sqrt(RHO) * W**a
                    * np.exp(1j * (qs[a, 0] * box.X + qs[a, 1] * box.Y))
                    for a in range(3)])
    V = L * L
    e_clock = (4.5 * G0 - 1.5 * G1) * RHO**2 * V
    excess = (box.energy(psi) - e_clock) / V
    predicted = 0.5 * RHO * np.sum(qs**2) + 1.5 * G1 * RHO**2

    def residual(g1):
        hp = box.h_psi(psi, g1)
        mu = np.sum(np.conj(psi) * hp, axis=(1, 2)).real / np.sum(
            np.abs(psi)**2, axis=(1, 2))
        return float(np.linalg.norm(hp - mu[:, None, None] * psi)
                     / np.linalg.norm(psi))

    res_full, res_half = residual(G1), residual(G1 / 2)
    n_tot = np.sum(np.abs(psi)**2)
    rng = np.random.default_rng(seed)
    p = psi * (1 + 1e-3 * rng.standard_normal(psi.shape))
    dt = 0.02
    for _ in range(steps):
        p = p - dt * np.einsum("...ij,j...->i...",
                               interaction_matrix(p, G0, G1), p)
        p = np.fft.ifft2(np.fft.fft2(p, axes=(1, 2)) / (1 + 0.5 * dt * box.K2),
                         axes=(1, 2))
        p *= math.sqrt(n_tot / np.sum(np.abs(p)**2))
    return dict(L=L, excess=excess, predicted=predicted,
                residual=res_full, residual_half_g1=res_half,
                windings_start=box.winding(psi), windings_end=box.winding(p),
                relaxed_excess=(box.energy(p) - e_clock) / V,
                min_rho_end=float(np.min(np.abs(p)**2)))


def triad(beta):
    n = np.ones(3) / math.sqrt(3)
    u = np.array([1., -1, 0]) / math.sqrt(2)
    v = np.cross(n, u)
    return np.array([math.cos(beta) * n + math.sin(beta)
                     * (math.cos(2 * math.pi * a / 3) * u
                        + math.sin(2 * math.pi * a / 3) * v) for a in range(3)])


def triad_spectrum(p, E, kappa, kind, phi, A, B):
    """Positive BdG frequencies with an explicit lock term.

    kind 'kinetic': kappa sum_a (e_a.grad conj psi_a)(e_a.grad psi_a).
    kind 'density': kappa sum_a (e_a.grad rho_a)^2 (boost invariant).
    """
    Aa = A + np.eye(3) * (p @ p) / 2
    Bb = B.copy()
    for a in range(3):
        s = (E[a] @ p)**2
        if kind == "kinetic":
            Aa[a, a] += kappa * s
        else:
            Aa[a, a] += 2 * kappa * s * abs(phi[a])**2
            Bb[a, a] += 2 * kappa * s * phi[a]**2
    H = np.block([[Aa, Bb], [-Bb.conj(), -Aa.conj()]])
    e = np.linalg.eigvals(H).real
    return np.sort(e[e > 1e-12])


def triad_locks(kappa=0.3, seed=2):
    phi, mu, r, A, B = clock_background()
    rng = np.random.default_rng(seed)
    dirs = rng.standard_normal((120, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    out = {}
    for kind in ("kinetic", "density"):
        for label, beta in (("orthonormal (T)", math.acos(1 / math.sqrt(3))),
                            ("generic C3", 0.5)):
            E = triad(beta)
            gap = triad_spectrum(np.zeros(3), E, kappa, kind, phi, A, B)[1]
            row = {}
            for p in (0.01, 0.03, 0.1):
                c = [triad_spectrum(p * d, E, kappa, kind, phi, A, B)[0] / p
                     for d in dirs]
                m = [(triad_spectrum(p * d, E, kappa, kind, phi, A, B)[1]
                      - gap) / p**2 for d in dirs]
                row[p] = (max(c) / min(c) - 1, max(m) / min(m) - 1)
            slope = math.log(row[0.03][0] / row[0.01][0]) / math.log(3.0)
            out[(kind, label)] = dict(phonon=row[0.1][0], relative=row[0.1][1],
                                      phonon_power=slope)
    return out


def pwave(L=8.0, N=32):
    box = Box(L, N)
    Q = 2 * math.pi / L
    qs = np.array([(Q, 0), (0, Q), (-Q, -Q)])

    def field(extra):
        return np.array([math.sqrt(RHO) * W**a * np.exp(
            1j * ((qs[a, 0] + extra[0]) * box.X + (qs[a, 1] + extra[1]) * box.Y))
            for a in range(3)])

    def pair_energy(psi):
        gx, gy = box.grad(psi)
        total = 0.0
        for c in range(3):
            a, b = (c + 1) % 3, (c + 2) % 3
            for g in (gx, gy):
                total += np.sum(np.abs(psi[a] * g[b] - psi[b] * g[a])**2)
        return float(total) * box.dx**2 / L**2

    rest = pair_energy(field((0, 0)))
    boosted = pair_energy(field((Q, 2 * Q)))
    predicted = RHO**2 * sum(np.sum((qs[a] - qs[b])**2)
                             for a, b in ((0, 1), (1, 2), (2, 0)))
    common = pair_energy(np.array([math.sqrt(RHO) * W**a * np.exp(
        1j * Q * box.X) for a in range(3)]))
    return dict(pair_density=rest, predicted=predicted,
                boost_change=abs(boosted - rest), common_twist=common,
                lambda_critical=-1 / (6 * RHO))


def run(report=print):
    s = schur_check()
    for name, v in s.items():
        report(f"{name:11s}: {v['dim']:2d} invariants, neutral anisotropy "
               f"{v['neutral_aniso']:.1e}, per component {v['component_aniso']:.1e}")
    assert s["Z3 x SO(3)"]["dim"] == 3
    assert s["Z3 x SO(3)"]["component_aniso"] < 1e-10
    assert s["twisted C3"]["neutral_aniso"] > 0.1
    assert s["twisted T"]["neutral_aniso"] < 1e-10
    assert s["twisted T"]["component_aniso"] > 0.1

    b = rest_frame_rank()
    report(f"Rest frame counterflow tensor: max s_min/s_max = "
           f"{b['rest_min_over_max']:.1e} (rank <= 2). Lab frame generic "
           f"rank 3 (median ratio {b['lab_median_ratio']:.2f}); Z3 counterflow "
           f"T = (3/2) rho Q^2 (1 - n n) to {b['z3_uniaxial_error']:.1e}.")
    assert b["rest_min_over_max"] < 1e-12 and b["z3_uniaxial_error"] < 1e-12

    pb = pointwise_bound()
    report(f"g1 density / rho_tot^2 >= {pb['min_ratio']:.6f} over samples "
           f"(bound -1/6 = {-1/6:.6f}, clock {pb['clock_ratio']:.6f}); "
           f"Wolstenholme slack >= {pb['wolstenholme_slack']:.1e}.")
    assert pb["min_ratio"] >= -1 / 6 - 1e-12
    assert pb["wolstenholme_slack"] >= -1e-12
    assert abs(pb["clock_ratio"] + 1 / 6) < 1e-12

    flows = [counterflow(8.0, 32, 12000), counterflow(16.0, 64, 12000)]
    for f in flows:
        report(f"Counterflow L={f['L']:.0f}: excess {f['excess']:.6f} "
               f"(predicted {f['predicted']:.6f}). Residual {f['residual']:.3f}, "
               f"at g1/2 {f['residual_half_g1']:.3f}. Windings "
               f"{f['windings_start']} -> {f['windings_end']}. Relaxed excess "
               f"{f['relaxed_excess']:.4f} (x L = {f['relaxed_excess'] * f['L']:.2f})")
        assert abs(f["excess"] - f["predicted"]) < 1e-9
        assert 1.8 < f["residual"] / f["residual_half_g1"] < 2.2
        assert f["windings_end"] != f["windings_start"]
        assert f["relaxed_excess"] > 0
    assert flows[1]["relaxed_excess"] < flows[0]["relaxed_excess"]

    t = triad_locks()
    for (kind, label), v in t.items():
        report(f"{kind:7s} lock, {label:15s}: phonon aniso {v['phonon']:.2e} "
               f"(p^{v['phonon_power']:.2f} for p <= 0.03), relative branch "
               f"curvature aniso {v['relative']:.2e} at p = 0.1")
    assert t[("kinetic", "generic C3")]["phonon"] > 0.05
    assert t[("density", "generic C3")]["relative"] > 0.05
    assert abs(t[("kinetic", "orthonormal (T)")]["phonon_power"] - 2) < 0.2
    assert abs(t[("density", "generic C3")]["phonon_power"] - 2) < 0.2

    pw = pwave()
    report(f"p-wave pair density {pw['pair_density']:.6f} (predicted "
           f"{pw['predicted']:.6f}). Boost change {pw['boost_change']:.1e}; "
           f"common twist {pw['common_twist']:.1e}. Relative stiffness "
           f"vanishes at lambda = {pw['lambda_critical']:.4f}.")
    assert abs(pw["pair_density"] - pw["predicted"]) < 1e-9
    assert pw["boost_change"] < 1e-9 and pw["common_twist"] < 1e-9

    report("RESULT: no invariant term locks index to space. No lock state is "
           "isotropic, homogeneous, or below the clock. Axial locks leave "
           "O(kappa) rank-2 anisotropy. Tetrahedral locks push it to rank "
           "four and carry no rank-2 tensor.")
    return dict(schur=s, boost=b, bound=pb, counterflow=flows,
                triad=t, pwave=pw)


if __name__ == "__main__":
    run()
