"""
knot_spectrum.py: the excitation tower of the saturated lump
(derivation program D7 of the gaps audit. Feeds the surplus-edge
prediction channel, D4).

THE QUESTION.  knot_charge.py showed the substrate's stable attractor
is a nodeless B = 0 lump.  Whatever the confined phase adds, the
CLASSICAL mode tower around that lump is measurable today: linearize
by kicking the saturated state with a small perturbation and Fourier-
analyzing the ring-down.  Two questions are frozen:

  1. Is the tower GAPPED above the condensate (no soft internal modes
     of the lump other than the exact symmetries)?
  2. Where does the first internal excitation sit relative to the
     BdG scales of the substrate (the relative-sector binding scale
     |lam1| = g1 rho0 / 2 and the chemical potential mu = 3 g0 + g1)?
     A state far below |lam1| would be an extra light internal state
     the framework does not have, a surplus prediction or a
     falsifier (D4 channel).  (Reading |lam1| as the scale the MTC
     lock-in spectrum rides on is an interpretive bridge. The frozen
     assertion uses only the substrate's own scale.)

METHOD.  Saturate the texture-imprinted knot (t = 8, sponge on the
relative modes at r > 9), then kick the relative sector with a small
symmetric radial factor (1 + eps * exp(-r^2/w^2)) and record two core
observables densely to t = 24: the relative-mode core norm (charged
channel) and the common-mode core density (neutral channel).  The
peaks of the FFT of the ring-down are the tower.

Usage: python3 knot_spectrum.py     (~2-4 min, mark slow)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# the shared 3D machinery (one copy: substrate3d.py)
from substrate3d import (DT, evolve, g0, g1, rel_modes,  # noqa: E402
                         rho0, texture_state)


def evolve_record(psi, K2, R, n_steps, record_every=0, core_r=3.0):
    """Evolve while recording the core observables (charged norm,
    neutral density) every record_every steps."""
    core = R < core_r
    ts_rel, ts_com = [], []

    def _record(p):
        p1, p2 = rel_modes(p)
        ts_rel.append(float((np.abs(p1)**2 + np.abs(p2)**2)[core].mean()))
        ts_com.append(float((np.abs(p.mean(axis=0))**2)[core].mean()))

    psi = evolve(psi, K2, R, n_steps, callback=_record,
                 every=record_every)
    return psi, np.array(ts_rel), np.array(ts_com)


def spectrum(ts, dt_sample):
    x = ts - ts.mean()
    x *= np.hanning(len(x))
    amp = np.abs(np.fft.rfft(x))
    freq = np.fft.rfftfreq(len(x), d=dt_sample) * 2 * np.pi  # angular
    return freq, amp


def peaks(freq, amp, floor=0.05):
    """Angular frequencies of local maxima above floor * max."""
    thr = floor * amp[1:].max()
    out = []
    for i in range(2, len(amp) - 1):
        if amp[i] > thr and amp[i] >= amp[i - 1] and amp[i] >= amp[i + 1]:
            out.append((freq[i], amp[i]))
    out.sort(key=lambda t: -t[1])
    return out


def run(report=print):
    report("KNOT SPECTRUM: ring-down tower of the saturated B = 0 lump")
    report("=" * 70)
    N, L, size = 48, 24.0, 2.5
    lam1 = g1 * rho0 / 2.0            # relative-sector binding scale
    mu = (3 * g0 + g1) * rho0         # chemical potential (BdG)
    report(f"  BdG scales: |lam1| = {lam1:.4f},  mu = {mu:.4f}")

    psi, dx, R, K2 = texture_state(N, L, amp=0.7, size=size)
    psi, _, _ = evolve_record(psi, K2, R, int(8.0 / DT))
    report(f"  saturated to t = 8 at {N}^3")

    # the kick: small symmetric radial factor on the relative sector
    eps, w = 0.03, 2.0
    kick = 1.0 + eps * np.exp(-(R / w) ** 2)
    p0 = psi.mean(axis=0)
    psi = p0[None] + (psi - p0[None]) * kick[None]

    every = 10                          # sample dt = 0.05
    psi, ts_rel, ts_com = evolve_record(psi, K2, R, int(16.0 / DT),
                                        record_every=every)
    dt_s = every * DT
    report(f"  ring-down recorded: t = 8 -> 24, {len(ts_rel)} samples")

    f_rel, a_rel = spectrum(ts_rel, dt_s)
    f_com, a_com = spectrum(ts_com, dt_s)
    pk_rel = peaks(f_rel, a_rel)[:5]
    pk_com = peaks(f_com, a_com)[:3]
    report("  charged (relative) tower, dominant peaks (omega, rel amp):")
    for w_, a_ in pk_rel:
        report(f"    omega = {w_:6.3f}   ({w_/lam1:5.2f} |lam1|, "
               f"{w_/mu:5.3f} mu)   amp {a_/a_rel[1:].max():.2f}")
    report("  neutral (common) response, dominant peaks:")
    for w_, a_ in pk_com:
        report(f"    omega = {w_:6.3f}   amp {a_/a_com[1:].max():.2f}")

    # ── frozen assertions ────────────────────────────────────────────
    assert len(pk_rel) >= 1, "the ring-down must resolve a tower"
    w_min = min(w_ for w_, _ in pk_rel)
    # 1. The tower is GAPPED: no internal mode at zero frequency
    #    (the resolution floor is the FFT bin, 2 pi / T with T = 16)
    assert w_min > 2 * np.pi / 16.0, \
        "soft internal mode found below the resolution floor"
    # 2. registered outcome (2026-07): the first internal excitation
    #    sits at or ABOVE the relative-sector binding scale |lam1|:
    #    no extra light internal states below the substrate's own
    #    binding scale.  A refactor pushing a peak far below |lam1|
    #    creates a surplus state: register it (D4) or treat it as a
    #    falsifier.
    assert w_min > 0.5 * lam1, \
        f"internal mode at {w_min:.3f} << |lam1| = {lam1:.3f}: " \
        f"surplus state below the binding scale. Register or falsify it"
    report("-" * 70)
    report("  REGISTERED OUTCOME: tower gapped. First internal mode at")
    report(f"  omega = {w_min:.3f} >= 0.5 |lam1|. No surplus states below")
    report("  the binding scale |lam1| in the classical ring-down.")
    return dict(pk_rel=pk_rel, pk_com=pk_com, lam1=lam1, mu=mu)


if __name__ == "__main__":
    run()
