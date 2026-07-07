"""
knot.py: does matter FORM in the true substrate?

Uses the upstream equation exactly (interference/nls_soliton.py):

    i dpsi_k/dt = -1/2 psi_k'' + [g0*rho_total
                   + g1*Re(conj(psi_{k+1}) psi_{k+2})] psi_k

with the G2 constraint g1 = g0/sqrt(2).  The BdG analysis upstream
shows the relative modes carry lambda_1 = -g1*rho0/2 < 0: a
MODULATIONAL INSTABILITY, the framework's claimed matter-formation
mechanism.  This probe tests it nonlinearly, with two experiments:

  EXP 1 (formation + control): seed a localized WINDING lump
    (psi_k ~ 1 + eps*exp(2pi i k/3)*G(x), a circulating loop in
    component space).  Under the true coupling, does it grow and
    localize (knot forms)?  CONTROL: same seed with g1 = 0, the
    instability is absent, so dispersal there proves the G2 exchange
    term is what binds.

  EXP 2 (selectivity): winding lump (circulating) vs v-biased lump
    (static mark, v = (1,-1/2,-1/2)).  The narrative requires the
    circulating loop to persist and the static mark to radiate.

Usage: python3 knot.py   (~30 s)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
from nls_soliton import _split_step, _z3_fourier, G1_OVER_G0


def _measures(psi, x, c):
    rho = np.abs(psi)**2
    z = _z3_fourier(rho)
    rel = z['relative']
    peak = float(rel.max())
    half = rel > 0.5 * peak
    width = float(half.sum() * (x[1] - x[0])) if peak > 0 else 0.0
    common_dip = float((z['common'][c] - z['common'].mean()).mean())
    return peak, width, common_dip


def evolve_report(psi0, g0, g1, x, dx, dt, t_marks, label, report=print):
    psi = psi0.copy()
    c = np.abs(x) < 8.0
    report(f"  {label}")
    report(f"    {'t':>5s} {'rel peak':>10s} {'width':>8s} "
           f"{'common dip':>12s} {'norm':>10s}")
    t = 0.0
    p0 = None
    for tm in t_marks:
        n = int(round((tm - t) / dt))
        if n:
            psi = _split_step(psi, dx, dt, g0, g1, n_steps=n)
        t = tm
        peak, width, dip = _measures(psi, x, c)
        norm = float((np.abs(psi)**2).sum() * dx)
        if p0 is None:
            p0 = peak
        report(f"    {t:5.1f} {peak:10.5f} {width:8.2f} "
               f"{dip:+12.6f} {norm:10.3f}")
    return peak / p0, psi


def run(report=print):
    N, L = 512, 100.0
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]
    dt = 0.002
    g0 = 1.0
    g1 = g0 * G1_OVER_G0
    rho0 = 1.0
    G = np.exp(-x**2 / 8.0)
    eps = 0.3
    t_marks = [0, 2, 4, 8, 12]

    # winding lump: circulating loop in component space (q = 1)
    w3 = np.exp(2j * np.pi * np.arange(3) / 3)
    psi_wind = np.sqrt(rho0) * (1.0 + eps * w3[:, None] * G[None, :]) \
        .astype(complex)
    # v-biased lump: static mark
    v = np.array([1.0, -0.5, -0.5])
    psi_bias = np.sqrt(rho0) * (1.0 + eps * v[:, None] * G[None, :]) \
        .astype(complex)

    report("KNOT PROBE: the true substrate (exchange coupling)")
    report("=" * 64)
    report("EXP 1: winding lump, true coupling vs control (g1 = 0)")
    growth_true, _ = evolve_report(psi_wind, g0, g1, x, dx, dt,
                                   t_marks, "true (g1 = g0/sqrt(2)):",
                                   report)
    growth_ctrl, _ = evolve_report(psi_wind, g0, 0.0, x, dx, dt,
                                   t_marks, "control (g1 = 0):", report)
    report(f"  relative-peak growth: true x{growth_true:.2f} vs "
           f"control x{growth_ctrl:.2f}")
    report(f"  => knot formation "
           f"{'OBSERVED (G2 coupling binds)' if growth_true > 1.5 > growth_ctrl else 'NOT OBSERVED'}")

    report("")
    report("EXP 2: circulating winding vs static v-mark (selectivity)")
    g_wind, _ = evolve_report(psi_wind, g0, g1, x, dx, dt, t_marks,
                              "winding (circulating):", report)
    g_bias, _ = evolve_report(psi_bias, g0, g1, x, dx, dt, t_marks,
                              "v-mark (static bias):", report)
    report(f"  winding growth x{g_wind:.2f} vs mark growth x{g_bias:.2f}")
    sel = g_wind / max(g_bias, 1e-9)
    verdict = ('SELECTIVE: circulation survives, mark does not'
               if sel > 1.5 else 'not selective at this amplitude/time')
    report(f"  selectivity ratio (winding/mark): {sel:.2f}  {verdict}")
    return {"growth_true": growth_true, "growth_ctrl": growth_ctrl,
            "selectivity": sel}


if __name__ == "__main__":
    run()
