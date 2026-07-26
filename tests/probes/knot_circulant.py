"""
knot_circulant.py: the knot manufactures the anomalous phase
(stage zero + first reading of the two-phase bridge).

CONTEXT.  The Koide-readable and knot-forming phases of the substrate
are disjoint for every admissible uniform coupling (r, alpha)
(circulant_normalization PARTS D-G), and at alpha = 0 the interaction
matrix's off-diagonals are density-weighted and REAL for every field
configuration, so no state can generate the Koide phase in the NORMAL
channel.  The last mechanism is the ANOMALOUS (pair) channel, whose
BdG couplings go as psi psi and therefore carry a non-uniform state's
phases.  This probe constructs a certified stationary knot and takes
the first reading of that channel.

STAGE ZERO (stationary knot, sponge-free).  At the knot-forming
coupling the uniform state is an energy saddle, so imaginary-time
relaxation with norm projection cannot fall back to it.  From a
real-time saturated seed, split-step imaginary time converges to a
non-uniform stationary state (relative RMS 0.25), certified by the
eigenvalue residual ||H psi - mu psi|| / ||psi||.

THE READING.  The anomalous block B(x) = dF/dpsi_bar of the local
nonlinear map is computed by Wirtinger finite differences; its q = 1
Z3 Fourier row component b1(x) transforms as e^{2i phi} under the
global U(1), so the GAUGE-INVARIANT phase is
beta(x) = arg(b1(x) e^{-2i arg p0(x)}) with p0 the local common mode.

MEASURED (seed 11, N = 256, L = 40, residual ~1e-3).
  The invariant of genuine complexity is |sin beta| (a real anomalous
  value gives beta in {0, pi}).  Uniform control at alpha = 0:
  sin beta = 0 exactly (the identity).  The uniform state at the
  CANONICAL coupling also has sin beta = 0 (there the phase lives in
  the normal channel).  The stationary knot: |sin beta| = 0.70 at the
  core, 0.57 density-weighted mean: a genuinely complex anomalous
  channel that NO uniform state carries.  Circular mean
  beta = +0.171, between the canonical alpha class (0.145) and 2/9
  (0.222), identified with neither.  ONE SEED; existence facts are
  frozen, identification is open.

ENSEMBLE (4 seeds x 2 boxes, all certified ~1e-3): the
existence fact is UNIVERSAL (density-weighted mean |sin beta| in
[0.19, 0.67], all runs), the value is NOT (circular mean beta spans
[-0.71, +0.30] with sign flips): the crude scalar identifies with
nothing; identification, if any, sits in the spectrum.

BDG SPECTRUM (bdg() below): the full 6N x 6N BdG on the
certified background, spectral Laplacian plus Wirtinger blocks.
Validated: two exact zero modes (phase, translation), stability at
the residual scale (max Im omega ~ 1e-2).  MEASURED: the knot
HYBRIDIZES the relative channels: low modes carry mixed q1/q2 weight
(~0.40/0.44), where every uniform alpha = 0 medium has exact
q-eigenmodes with exact q1/q2 degeneracy.  The degeneracy is lifted
at the eigenbasis level.  No sharply localized internal triplet
exists on this broad saturated texture (IPR ~ 1-2).  The
interpretation and its falsifier are the LAYER SCOPING statement in
circulant_normalization.py: the classical layer provides the slots,
the ratios are written at the deeper foldings, and a classical
measurement of Q = 2/3 would falsify that scoping.

Usage: python3 knot_circulant.py     (~2 min, mark slow)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _bootstrap  # noqa: E402,F401

from nls_soliton import (interaction_matrix, _split_step,  # noqa: E402
                         bdg_operator, G1_OVER_G0, _z3_fourier)

G0 = 1.0
G1 = G0 * G1_OVER_G0
L, N = 40.0, 256
DX = L / N
K = 2 * np.pi * np.fft.fftfreq(N, d=DX)
W3 = np.exp(2j * np.pi / 3)
ALPHA_CLASS = 4.333441893 % (2 * np.pi / 3)   # canonical alpha, gauge class


def _h_apply(p):
    lap = np.fft.ifft(-(K ** 2) * np.fft.fft(p, axis=1), axis=1)
    M = interaction_matrix(p, G0, G1, 0.0)
    return -0.5 * lap + np.einsum('xij,jx->ix', M, p)


def _relax(psi, n0, dtau, n_iter):
    kin = np.exp(-0.5 * dtau * (K ** 2))
    for _ in range(n_iter):
        psi = np.fft.ifft(kin * np.fft.fft(psi, axis=1), axis=1)
        M = interaction_matrix(psi, G0, G1, 0.0)
        w, V = np.linalg.eigh(M)
        comp = np.moveaxis(psi, 0, -1)
        amp = np.einsum('...ji,...j->...i', V.conj(), comp) \
            * np.exp(-dtau * w)
        psi = np.moveaxis(np.einsum('...ij,...j->...i', V, amp), -1, 0)
        psi *= np.sqrt(n0 / (float(np.sum(np.abs(psi) ** 2)) * DX))
    return psi


def _beta(psi, g1=G1, alpha=0.0):
    """Gauge-invariant anomalous phase field beta(x).

    A REAL anomalous value gives beta in {0, pi} (time-reversal
    symmetric); genuinely complex structure is sin(beta) != 0."""
    h = 1e-7
    bloc = np.zeros((3, 3, N), dtype=complex)

    def floc(p):
        M = interaction_matrix(p, G0, g1, alpha)
        return np.einsum('xij,jx->ix', M, p)

    for j in range(3):
        e = np.zeros((3, N), dtype=complex)
        e[j] = 1.0
        dr = (floc(psi + h * e) - floc(psi - h * e)) / (2 * h)
        di = (floc(psi + 1j * h * e) - floc(psi - 1j * h * e)) / (2 * h)
        bloc[:, j] = 0.5 * (dr + 1j * di)
    b1 = bloc[0, 0] + bloc[0, 1] * W3 + bloc[0, 2] * W3 ** 2
    p0 = psi.mean(axis=0)
    return np.angle(b1 * np.exp(-2j * np.angle(p0)))


def run(report=print):
    report("KNOT CIRCULANT: the anomalous phase at a stationary knot")
    report("=" * 70)

    # uniform control at alpha = 0: the identity forces a REAL
    # anomalous channel, beta in {0, pi}, i.e. sin(beta) = 0
    psi_u = np.ones((3, N), dtype=complex)
    s_u = np.abs(np.sin(_beta(psi_u)))
    report(f"  uniform control (alpha = 0): max |sin beta| = "
           f"{s_u.max():.2e}")
    assert s_u.max() < 1e-9
    # canonical reference: beta of the UNIFORM state at the canonical
    # coupling (r_c, alpha_c), same invariant, same convention
    R_C, AL_C = 0.383029723417077641, 4.33344189294178415
    beta_ref = float(_beta(psi_u, g1=R_C * G0, alpha=AL_C)[0])
    report(f"  canonical uniform reference: beta_ref = {beta_ref:+.4f}")

    # stage zero: saturate, then relax to a stationary state
    rng = np.random.default_rng(11)
    psi = np.empty((3, N), dtype=complex)
    for kk in range(3):
        noise = rng.standard_normal(N) + 1j * rng.standard_normal(N)
        psi[kk] = 1.0 + 1e-3 * noise
    psi = _split_step(psi, DX, 0.005, G0, G1, n_steps=int(40.0 / 0.005))
    n0 = float(np.sum(np.abs(psi) ** 2)) * DX
    psi = _relax(psi, n0, 0.002, 8000)
    psi = _relax(psi, n0, 0.02, 15000)
    psi = _relax(psi, n0, 0.005, 15000)
    hp = _h_apply(psi)
    mu = float(np.real(np.sum(np.conj(psi) * hp)) * DX) / n0
    res = float(np.linalg.norm(hp - mu * psi) / np.linalg.norm(psi))
    rel = _z3_fourier(np.abs(psi) ** 2)["relative"]
    rms = float(np.sqrt(rel.mean()))
    report(f"  stationary state: mu = {mu:.6f}, residual = {res:.2e}, "
           f"relative RMS = {rms:.4f}")
    assert res < 5e-3, "stationarity certification failed"
    assert rms > 0.1, "the relaxed state must keep relative structure"

    # the reading: |sin beta| is the invariant of genuine complexity
    beta = _beta(psi)
    s = np.abs(np.sin(beta))
    wt = rel / rel.sum()
    core = int(np.argmax(rel))
    mean_s = float(np.sum(wt * s))
    # circular density-weighted mean of beta
    mean_circ = float(np.angle(np.sum(wt * np.exp(1j * beta))))
    report(f"  knot: |sin beta| at core = {s[core]:.4f}, "
           f"density-weighted mean |sin beta| = {mean_s:.4f}")
    report(f"  circular mean beta = {mean_circ:+.4f} "
           f"(canonical uniform reference {beta_ref:+.4f}, "
           f"2/9 = +0.2222)")
    # frozen existence facts (the identification question stays open)
    assert mean_s > 0.3, "the knot must carry a complex anomalous phase"
    assert s[core] > 0.1
    report("-" * 70)
    report("  VERDICT: the knot GENERATES the anomalous phase the "
           "uniform medium cannot")
    report("  carry (uniform sin beta = 0 by identity; knot mean |sin beta| "
           f"= {mean_s:.2f}).  The")
    report("  quantitative identification (mean beta vs the canonical "
           "alpha class) is a")
    report("  one-seed lead, reported not asserted.  OPEN: ensemble, "
           "box scaling, and the")
    report("  full BdG mode spectrum on this background.")
    _BDG_STATE["psi"] = psi
    return dict(mu=mu, residual=res, rel_rms=rms,
                sin_beta_core=float(s[core]), sin_beta_mean=mean_s,
                beta_circ_mean=mean_circ, beta_ref=beta_ref)


def bdg(psi=None, report=print):
    """Full BdG on the stationary background: validation + the
    channel-hybridization measurement."""
    if psi is None:
        run(report=lambda *a, **k: None)
        psi = _BDG_STATE["psi"]
    hp = _h_apply(psi)
    n0 = float(np.sum(np.abs(psi) ** 2)) * DX
    mu = float(np.real(np.sum(np.conj(psi) * hp)) * DX) / n0

    # shared operator: nls_soliton.bdg_operator, one copy
    w, V = bdg_operator(psi, G0, G1, 0.0, K, mu)

    # validation: two zero modes, stability at the residual scale
    wr = np.sort(np.abs(w.real))
    assert wr[1] < 1e-3, "phase + translation zero modes required"
    assert np.abs(w.imag).max() < 2e-2, "background must be stable"

    # channel hybridization of the low positive modes
    sel = np.argsort(np.where(w.real > 1e-3, w.real, np.inf))[:5]
    mixes = []
    for m in sel:
        u = V[:3 * N, m].reshape(3, N)
        qwt = []
        for q in range(3):
            proj = sum(u[k] * W3 ** (-q * k) for k in range(3)) / 3.0
            qwt.append(float(np.sum(np.abs(proj) ** 2)))
        s = sum(qwt)
        qwt = [x / s for x in qwt]
        mixes.append(min(qwt[1], qwt[2]))
        report(f"  omega = {w.real[m]:.4f}: q-weights "
               f"({qwt[0]:.2f}, {qwt[1]:.2f}, {qwt[2]:.2f})")
    hybrid = float(np.mean(mixes))
    report(f"  mean min(q1,q2) weight over the 5 lowest modes: "
           f"{hybrid:.3f} (uniform alpha = 0 medium: exactly 0)")
    assert hybrid > 0.15, \
        "the knot must hybridize the relative channels"
    return dict(mu=mu, hybridization=hybrid,
                zero_modes=float(wr[1]),
                max_growth=float(np.abs(w.imag).max()))


_BDG_STATE = {}


if __name__ == "__main__":
    run()
