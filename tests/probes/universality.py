"""
universality.py: is the standing wave THE survivor, regardless of seed?

Completes the selectivity question at the level the earlier probes
located it (the SATURATED state): evolve the circulating winding seed
and the static v-mark seed to t = 40 and compare their endpoints:
field-level winding shares, internal frequency, peak.  If both land
on the same class of standing structure, the survivor is UNIQUE and
seed-independent: every loop forgets its label and becomes the same
standing wave.  That is the strongest dynamical form of the
forgetting/selection claim available in 1D.

Usage: python3 universality.py  (~35 s)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
from nls_soliton import _split_step, _z3_fourier, G1_OVER_G0


def endpoint(psi0, g0, g1, x, dx, dt, T=40.0):
    psi = psi0.copy()
    omega = np.exp(2j*np.pi/3)
    core = psi.shape[1] // 2
    mu = (3*g0 + g1) * 1.0
    # evolve to T, sampling core phase over the last 20 time units
    n_pre = int(20.0/dt)
    psi = _split_step(psi, dx, dt, g0, g1, n_steps=n_pre)
    phases, times, t = [], [], 20.0
    for _ in range(200):
        psi = _split_step(psi, dx, dt, g0, g1, n_steps=50)
        t += 50*dt
        a1 = (psi[0, core] + np.conj(omega)*psi[1, core]
              + np.conj(omega)**2*psi[2, core]) / 3.0
        phases.append(np.angle(a1))
        times.append(t)
    w_rel = float(np.polyfit(times, np.unwrap(phases), 1)[0]) + mu
    rho = np.abs(psi)**2
    rel = _z3_fourier(rho)['relative']
    peak = float(rel.max())
    a1f = (psi[0] + np.conj(omega)*psi[1] + np.conj(omega)**2*psi[2])/3.0
    a2f = (psi[0] + np.conj(omega)**2*psi[1] + np.conj(omega)*psi[2])/3.0
    s1 = float((np.abs(a1f)**2).sum())
    s2 = float((np.abs(a2f)**2).sum())
    return {"peak": peak, "share1": s1/(s1+s2), "w_rel": w_rel}


def run(report=print):
    N, L = 384, 90.0
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]
    dt = 0.002
    g0, g1 = 1.0, G1_OVER_G0
    G = np.exp(-x**2/8.0)
    w3 = np.exp(2j*np.pi*np.arange(3)/3)
    v = np.array([1.0, -0.5, -0.5])
    seeds = {
        "winding (circulating)": (1.0 + 0.3*w3[:, None]*G[None, :]),
        "v-mark (static bias)":  (1.0 + 0.3*v[:, None]*G[None, :]),
    }
    report("UNIVERSALITY: do all seeds become the same standing wave?")
    report("=" * 64)
    report(f"  {'seed':<24s} {'peak(40)':>9s} {'q1 share':>9s} "
           f"{'omega_rel':>10s}")
    out = {}
    for name, s in seeds.items():
        e = endpoint(s.astype(complex), g0, g1, x, dx, dt)
        out[name] = e
        report(f"  {name:<24s} {e['peak']:9.3f} {e['share1']:9.3f} "
               f"{e['w_rel']:+10.4f}")
    a, b = out.values()
    same = (abs(a['share1'] - b['share1']) < 0.05
            and abs(a['w_rel'] - b['w_rel']) < 0.1
            and 0.3 < a['peak']/max(b['peak'], 1e-9) < 3.0)
    report("-" * 64)
    report(f"  same standing attractor: {'YES' if same else 'NO'} "
           f"(shares equal, frequencies match, peaks same order)")
    return same


if __name__ == "__main__":
    run()
