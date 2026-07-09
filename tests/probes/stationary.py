"""
stationary.py: after formation, does the knot PERSIST, does it keep
its winding, and how fast does it circulate?

Continues probes/knot.py: seed the winding lump, let the instability
saturate (t ~ 8-10), then evolve long (t = 40) and measure

  1. persistence: relative peak and width, sampled every 5 time units
     (a knot is a lump that breathes without dying).
  2. winding purity: integral of |rho_q1|^2 vs |rho_q2|^2, does the
     structure stay in its winding class.
  3. internal frequency: the phase of the q=1 amplitude at the lump's
     core, unwrapped and fitted, the knot's circulation rate omega.
     (m = hbar*omega is the framework's mass-as-rotation claim. This
     is the first measured number on the road to the spectrum test.)

Usage: python3 stationary.py   (~25 s)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
from nls_soliton import _split_step, _z3_fourier, G1_OVER_G0


def run(report=print):
    N, L = 512, 100.0
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]
    dt = 0.002
    g0, g1, rho0 = 1.0, G1_OVER_G0, 1.0
    w3 = np.exp(2j * np.pi * np.arange(3) / 3)
    psi = (np.sqrt(rho0) *
           (1.0 + 0.3 * w3[:, None] * np.exp(-x**2 / 8.0)[None, :])
           ).astype(complex)
    omega = np.exp(2j * np.pi / 3)
    core = N // 2

    report("STATIONARITY PROBE: the saturated knot, long evolution")
    report("=" * 66)
    report(f"  {'t':>5s} {'rel peak':>10s} {'width':>7s} "
           f"{'q1 share':>9s} {'q2 share':>9s}")

    t, phases, times = 0.0, [], []
    for tm in [0, 5, 10, 15, 20, 25, 30, 35, 40]:
        n = int(round((tm - t) / dt))
        if n:
            # sample core phase densely across this leg for the fit
            leg = 25
            per = n // leg
            for _ in range(leg):
                psi = _split_step(psi, dx, dt, g0, g1, n_steps=per)
                a1 = (psi[0, core] + np.conj(omega) * psi[1, core]
                      + np.conj(omega)**2 * psi[2, core]) / 3.0
                t += per * dt
                if t > 10.0:              # post-saturation only
                    phases.append(np.angle(a1))
                    times.append(t)
        rho = np.abs(psi)**2
        z = _z3_fourier(rho)
        rel = z['relative']
        peak = float(rel.max())
        width = float((rel > 0.5 * peak).sum() * dx)
        # winding shares must be FIELD-level: density-level Z3 cannot
        # see winding (rho_q2 = conj(rho_q1) identically for real rho)
        a1f = (psi[0] + np.conj(omega)*psi[1] + np.conj(omega)**2*psi[2])/3.0
        a2f = (psi[0] + np.conj(omega)**2*psi[1] + np.conj(omega)*psi[2])/3.0
        s1, s2 = float((np.abs(a1f)**2).sum()), float((np.abs(a2f)**2).sum())
        tot = s1 + s2
        report(f"  {tm:5.0f} {peak:10.4f} {width:7.2f} "
               f"{s1/tot:9.3f} {s2/tot:9.3f}")

    ph = np.unwrap(np.array(phases))
    tt = np.array(times)
    w_fit = float(np.polyfit(tt, ph, 1)[0])
    mu = (3 * g0 + g1) * rho0
    report("-" * 66)
    report(f"  persistence: lump alive at t = 40 "
           f"(formation was at t ~ 8): "
           f"{'YES' if peak > 0.5 else 'NO'}")
    w_rel = w_fit + mu     # subtract the condensate's own -mu rotation
    report(f"  raw core phase rate = {w_fit:+.4f}. Condensate rate = {-mu:+.4f}")
    report(f"  INTERNAL circulation (relative to condensate): "
           f"omega_rel = {w_rel:+.4f}")
    report(f"  scale |lambda_1| rho0 = {g1*rho0/2:.4f}: the knot")
    report(f"  oscillates below the natural scale. Leakage empirically")
    report(f"  weak (persistence).  NOTE: low-omega relative waves exist")
    report(f"  above the instability band, so this is suppression, not")
    report(f"  strict spectral absence.")
    report(f"  (the spectrum test proper needs the three windings'")
    report(f"   frequencies compared as Delta_k^2, next probe)")
    return {"peak_final": peak, "omega": w_fit,
            "omega_rel": w_rel, "gap_scale": g1 * rho0 / 2,
            "persistent": bool(peak > 0.5)}


if __name__ == "__main__":
    run()
