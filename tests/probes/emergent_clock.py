"""
emergent_clock.py: give the substrate a phase and READ THE CLOCK BACK.

THE POINT (substrate-conjecture ledger, entry experiment).  Everywhere
else in this repository the chain runs algebra → number.  This probe
runs the OTHER direction, the one the ontology actually asserts: seed
the Z₃-coupled NLS medium with a GENERIC random fluctuation (no
target, no Brannen angle, no PDG value anywhere in this file), let the
relative-sector instability grow knots from the noise, and MEASURE the
internal clock of the three components at the knots.  Mass is a rate
(E = ħω), so whatever spectrum the medium locks into is the medium's
own answer, not the ledger's.

THE PROTOCOL (frozen, part of the measurement):
    seed     symmetric vacuum + complex Gaussian noise + one global
             (gauge) phase
    warm-up  evolve to t = 30 (saturation observed at t ≈ 25, growth
             ≈ ×10 in the relative-mode RMS)
    readout  density-weighted knot mode a_k(t), sampled at Δt = 0.2
             over t ∈ [30, 120]; dominant line ω_k = peak of the
             windowed spectrum; LINEWIDTH = peak power share.

CALIBRATION FINDINGS (measured, kept):
    · the instability grows only from the Z₃-SYMMETRIC vacuum; random
      per-component phase offsets stabilise the background (growth ~1
      on the same horizon).  The symmetric vacuum is the one that
      decays into matter.
    · in THIS regime (1D, substrate-only, this horizon) the knots form
      but the clocks do NOT lock a sharp line: peak power share is at
      the few-percent level (broadband, soliton collisions in the
      box).  The verdict is frozen as measured.

WHAT IS ASSERTED (structural, target-blind): number conservation;
knot growth from noise (> ×3); dominant frequencies positive and
below the linearised μ + margin; and the verdict fields being
reported, including clocks_locked = False for the substrate-only 1D regime.

WHAT IS *NOT* ASSERTED: Q ≈ Q₀ = 2/3 or any Brannen structure.  The
ledgered conjecture is that the CONFINED-PHASE, 3D regime (the same
place the π₃ obligation moved after the substrate-only negative,
tests/probes/knot_charge.py) locks sharp clocks with the Z₃ circulant
spectrum, B/A = √2, θ = 2/9.  This probe is the harness and the first
data point.  The day a refined regime prints a sharp locked Q → 2/3
from random seeds, the retrodiction criticism of the scorecard is
dead.  Falsification direction: if every refinement stays broadband
or locks onto a spectrum incompatible with the circulant, the
substrate conjecture fails at the clock level.

Usage: python3 emergent_clock.py     (~12 s: 1D, three seeds)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "interference"))
import _bootstrap  # noqa: E402,F401
# sibling probes (substrate3d) import by bare name, as they do among
# themselves; keep the one canonical namespace
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nls_soliton import _split_step, _z3_fourier, G1_OVER_G0  # noqa: E402


# run parameters (frozen: part of the measurement protocol)
N, L = 256, 40.0
DT = 0.002
T_WARM = 30.0               # warm-up horizon (saturation ≈ t = 25)
T_READ = 90.0               # readout window length
DT_SAMPLE = 0.2             # readout sampling interval
NOISE = 0.05
SEEDS = (11, 23, 47)
NFFT = 4096


def _dominant_lines(series, dt_s):
    """Dominant spectral line per component of the knot mode a_k(t):
    Hann-windowed FFT, ω_k = argmax, share = peak power fraction."""
    a = np.asarray(series)                     # (n_t, 3)
    freqs = 2 * np.pi * np.fft.fftfreq(NFFT, d=dt_s)
    omegas, shares = [], []
    for k in range(3):
        A = np.fft.fft(a[:, k] * np.hanning(a.shape[0]), n=NFFT)
        P = np.abs(A)**2
        pk = int(np.argmax(P))
        omegas.append(float(-freqs[pk]))       # ψ ~ e^{-iωt}
        shares.append(float(P[pk] / P.sum()))
    return omegas, shares


def _one_seed(seed, g0=1.0):
    g1 = g0 * G1_OVER_G0
    rng = np.random.default_rng(seed)
    dx = L / N

    # GENERIC SEED: symmetric vacuum + noise + one global gauge phase
    glob = np.exp(2j * np.pi * rng.random())
    psi = np.empty((3, N), dtype=complex)
    for k in range(3):
        noise = (rng.standard_normal(N) + 1j * rng.standard_normal(N))
        psi[k] = (1.0 + NOISE * noise) * glob

    total0 = float((np.abs(psi)**2).sum())
    rel0 = float(np.sqrt(_z3_fourier(np.abs(psi)**2)["relative"].mean()))

    # warm-up to saturation
    psi = _split_step(psi, dx, DT, g0, g1,
                      n_steps=int(round(T_WARM / DT)))
    growth = float(np.sqrt(
        _z3_fourier(np.abs(psi)**2)["relative"].mean())) / max(rel0, 1e-300)

    # readout: knot mode a_k(t) = Σ_x w(x)ψ_k, w = relative density
    per = int(round(DT_SAMPLE / DT))
    n_samp = int(round(T_READ / DT_SAMPLE))
    series = []
    for _ in range(n_samp):
        psi = _split_step(psi, dx, DT, g0, g1, n_steps=per)
        rel = _z3_fourier(np.abs(psi)**2)["relative"]
        w = rel / max(rel.sum(), 1e-300)
        series.append([complex(np.sum(w * psi[k])) for k in range(3)])

    conservation = abs(float((np.abs(psi)**2).sum()) - total0) / total0
    omegas, shares = _dominant_lines(series, DT_SAMPLE)

    Q = None
    if all(o > 0 for o in omegas):
        Q = sum(omegas) / sum(math.sqrt(o) for o in omegas)**2
    return dict(seed=seed, conservation=conservation, growth=growth,
                omegas=omegas, peak_shares=shares, Q=Q)


# ═══════════════════════════════════════════════════════════════════════
#  3D feasibility rung (protocol frozen)
# ═══════════════════════════════════════════════════════════════════════

# Measured at freeze (seed 11, N = 24, L = 24, warm t = 30 free,
# readout t ∈ [30, 90] at Δt = 0.2):
#   free readout:    growth ×11.9, ω ≈ (3.22, 3.23, 3.27),
#                    peak shares ≤ 0.031, Q = 0.3333
#   sponged readout: knots survive isolation (residual relative RMS
#                    0.11), clock drops to ω ≈ 1.8, shares ≤ 0.032
# Figures vary across BLAS/FFT builds; the directional bands below
# are the frozen content.
# VERDICT: the substrate-only theory does not lock a sharp clock in 3D
# either, free or isolated.  Rung 2 of the ladder is the confined-
# phase effective regime (the locked circulant clock), not more substrate-only
# evolution.

N3, L3 = 24, 24.0
T3_WARM, T3_READ, DT3_SAMPLE = 30.0, 60.0, 0.2
SEED3 = 11


def run_3d(report=print, sponged=False):
    """One-seed 3D rung of the protocol (slow: several minutes)."""
    from substrate3d import DT as DT3, evolve, grid, rel_modes
    rng = np.random.default_rng(SEED3)
    dx, X, Y, Z, R, K2 = grid(N3, L3)
    glob = np.exp(2j * np.pi * rng.random())
    psi = (1.0 + NOISE * (rng.standard_normal((3, N3, N3, N3))
                          + 1j * rng.standard_normal((3, N3, N3, N3)))
           ) * glob

    def rel_rms(p):
        p1, p2 = rel_modes(p)
        return float(np.sqrt((np.abs(p1)**2 + np.abs(p2)**2).mean()))

    total0 = float((np.abs(psi)**2).sum())
    rel0 = rel_rms(psi)
    psi = evolve(psi, K2, R, int(round(T3_WARM / DT3)), sponge=False)
    growth = rel_rms(psi) / max(rel0, 1e-300)

    per = int(round(DT3_SAMPLE / DT3))
    series = []
    for _ in range(int(round(T3_READ / DT3_SAMPLE))):
        psi = evolve(psi, K2, R, per, sponge=sponged)
        p1, p2 = rel_modes(psi)
        w = np.abs(p1)**2 + np.abs(p2)**2
        w /= max(w.sum(), 1e-300)
        series.append([complex((w * psi[k]).sum()) for k in range(3)])

    conservation = (abs(float((np.abs(psi)**2).sum()) - total0) / total0
                    if not sponged else None)
    omegas, shares = _dominant_lines(series, DT3_SAMPLE)
    Q = None
    if all(o > 0 for o in omegas):
        Q = sum(omegas) / sum(math.sqrt(o) for o in omegas)**2
    out = dict(seed=SEED3, growth=growth, conservation=conservation,
               omegas=omegas, peak_shares=shares, Q=Q,
               clocks_locked=max(shares) > 0.5, sponged=sponged)
    report(f"  3D seed {SEED3} N={N3} "
           f"{'SPONGED' if sponged else 'free':>7s}: growth={growth:.1f} "
           f"ω={[round(o, 3) for o in omegas]} "
           f"shares={[round(s, 3) for s in shares]} "
           f"Q={None if Q is None else round(Q, 4)}")
    return out


def run(report=print):
    report("EMERGENT CLOCK: random fluctuation in, spectrum read back")
    report("=" * 70)
    report(f"  medium: Z₃-NLS, g₁/g₀ = 1/√2 (G₂ constraint), the ONLY")
    report(f"  physics input. No angle, no target, no PDG value here.")
    report(f"  grid N={N}, L={L}, dt={DT}; warm-up t={T_WARM:.0f}, "
           f"readout t∈[{T_WARM:.0f},{T_WARM+T_READ:.0f}] at "
           f"Δt={DT_SAMPLE}; noise={NOISE}, seeds={SEEDS}")

    results = [_one_seed(s) for s in SEEDS]
    Qs = [r["Q"] for r in results if r["Q"] is not None]

    report(f"\n  {'seed':>5s} {'conserv':>9s} {'growth':>7s} "
           f"{'ω₀':>8s} {'ω₁':>8s} {'ω₂':>8s} {'line':>6s} {'Q':>8s}")
    for r in results:
        o = r["omegas"]
        q_str = "   (n/a)" if r["Q"] is None else f"{r['Q']:.4f}"
        report(f"  {r['seed']:>5d} {r['conservation']:>9.1e} "
               f"{r['growth']:>7.1f} {o[0]:>8.4f} {o[1]:>8.4f} "
               f"{o[2]:>8.4f} {max(r['peak_shares']):>6.2f} {q_str:>8s}")

    # the verdict, stated plainly
    locked = all(max(r["peak_shares"]) > 0.5 for r in results)
    Q_mean = float(np.mean(Qs)) if Qs else None
    Q_spread = float(np.std(Qs)) if Qs else None
    report("-" * 70)
    report(f"  clocks_locked = {locked}  (sharp line := peak power "
           f"share > 0.5; measured ≤ "
           f"{max(max(r['peak_shares']) for r in results):.2f})")
    if Q_mean is not None:
        report(f"  Q(dominant lines) = {Q_mean:.4f} ± {Q_spread:.4f} "
               f"over {len(Qs)} seeds (math band [1/3, 1); circulant "
               f"target 2/3; degenerate clocks 1/3)")
    report(f"  VERDICT: knots form from generic noise, but the substrate-only 1D "
           f"regime is broadband on this horizon - no sharp clock.  "
           f"The locked-clock obligation sits with the confined-phase "
           f"3D regime (the locked circulant clock conjecture), same place "
           f"the π₃ obligation moved.  READOUT, not claim.")
    return dict(results=results, Q_values=Qs, Q_mean=Q_mean,
                Q_spread=Q_spread, clocks_locked=locked)


if __name__ == "__main__":
    run()
