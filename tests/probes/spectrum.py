"""
spectrum.py: what does the knot ring at?

Spectroscopy of the saturated lump: evolve past saturation, record
the winding amplitudes at the core over a long window, remove the
condensate's carrier rotation (e^{-i mu t}), FFT, and read the
discrete internal lines.

PRE-STATED MEANINGS (written before looking):
  (a) three lines with ratios ~ Delta_k of the circulant fringe
      (0.040 : 0.580 : 2.379, i.e. 1 : 14.4 : 59) -> the mass
      spectrum is dynamical: the framework's strongest possible
      result;
  (b) one line + its harmonics -> the 1D lump is a single-mode
      oscillator; the circulant spectrum lives at the MTC layer,
      and the literal-PDE reading is constrained (NOT a kill: the
      papers derive the spectrum at the topological layer);
  (c) broadband -> turbulent state, inconclusive; a true stationary
      profile is needed first.

Usage: python3 spectrum.py   (~35 s)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
from nls_soliton import _split_step, G1_OVER_G0

FRINGE = np.array([0.04034990821920736, 0.5802119201475366,
                   2.3794381716332555])          # Delta_k of the circulant


def run(report=print):
    N, L = 384, 90.0
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]
    dt = 0.002
    g0, g1, rho0 = 1.0, G1_OVER_G0, 1.0
    mu = (3*g0 + g1) * rho0
    w3 = np.exp(2j*np.pi*np.arange(3)/3)
    psi = (np.sqrt(rho0) *
           (1.0 + 0.3*w3[:, None]*np.exp(-x**2/8.0)[None, :])
           ).astype(complex)
    core = N // 2
    omega = np.exp(2j*np.pi/3)

    # reach saturation
    psi = _split_step(psi, dx, dt, g0, g1, n_steps=5000)   # t = 10

    # record window t = 10 .. 66
    every = 8
    n_rec = 3500
    ts, A1, A2 = [], [], []
    t = 10.0
    for i in range(n_rec):
        psi = _split_step(psi, dx, dt, g0, g1, n_steps=every)
        t += every*dt
        a1 = (psi[0, core] + np.conj(omega)*psi[1, core]
              + np.conj(omega)**2*psi[2, core]) / 3.0
        a2 = (psi[0, core] + np.conj(omega)**2*psi[1, core]
              + np.conj(omega)*psi[2, core]) / 3.0
        carrier = np.exp(1j*mu*t)          # remove condensate rotation
        ts.append(t)
        A1.append(a1*carrier)
        A2.append(a2*carrier)

    ts = np.array(ts)
    dts = ts[1] - ts[0]
    win = np.hanning(len(ts))

    def lines(series, n_top=6):
        S = np.fft.fftshift(np.fft.fft(np.array(series)*win))
        f = np.fft.fftshift(np.fft.fftfreq(len(ts), d=dts))*2*np.pi
        P = np.abs(S)**2
        P /= P.max()
        # top local maxima above threshold
        idx = [i for i in range(2, len(P)-2)
               if P[i] > P[i-1] and P[i] > P[i+1] and P[i] > 0.01]
        idx.sort(key=lambda i: -P[i])
        return [(float(f[i]), float(P[i])) for i in idx[:n_top]]

    report("KNOT SPECTROSCOPY: internal lines of the saturated lump")
    report("=" * 66)
    report(f"  window t = 10..{ts[-1]:.0f}, resolution "
           f"d_omega = {2*np.pi/(ts[-1]-ts[0]):.3f}")
    for name, series in [("a1 (winding +1)", A1), ("a2 (winding -1)", A2)]:
        report(f"  {name} core amplitude, dominant lines "
               f"(omega, rel. power):")
        for fr, p in lines(series):
            report(f"      omega = {fr:+7.3f}   power {p:8.3f}")
    report("-" * 66)
    report("  reference ratios: circulant fringe Delta_k = "
           f"{FRINGE[0]:.4f} : {FRINGE[1]:.4f} : {FRINGE[2]:.4f}")
    report("  (pre-stated meanings in the docstring; interpret against")
    report("   (a) three-line fringe / (b) single-mode / (c) broadband)")


if __name__ == "__main__":
    run()
