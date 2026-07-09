"""
skyrmion_3d.py: 3D knot FORMATION in the one-field substrate.

Marries the repo's two probes, in full 3D dynamics:

  knot.py (1D)            : Z3-NLS modulational-instability formation.
                            a winding seed GROWS under the true G2
                            coupling and disperses under g1 = 0.
  winding_texture.py (3D) : the explicit pi_3 skyrmion texture, static.

Here the winding texture is seeded at SMALL amplitude (eps = 0.3, a
knot.py-style Gaussian envelope) into the relative (q = 1, 2) sector
of the Z3-coupled NLS and evolved with the repo's exact equation
(g1 = g0/sqrt(2), no sponge: the run is norm-conserving and
conservation is tracked as a diagnostic).  This is the review
session's sim3d experiment, made a house probe.

FOUR MEASUREMENTS, each frozen in CI (values from the 48^3 run):

  1. FORMATION vs CONTROL.  The relative-mode rms grows x15.6 by
     t = 12 under the true coupling and the peak saturates near the
     condensate scale. Under g1 = 0 the rms stays flat (x1.00, no
     growth channel exists) and the peak collapses x0.009 (pure
     dispersal).  Matter forms BECAUSE of the G2 coupling, the 1D
     result (knot.py), now in 3D.  (Without a sponge the run
     conserves norm, so rms is the growth metric and peak the
     dispersal metric.)

  2. BdG GROWTH MATCH.  The e-folding rate of the relative rms in
     the growth window t in [3, 8] is 0.28 = 0.80 sigma, against the
     Bogoliubov instability rate sigma = g1 rho0 / 2 = 0.354
     (nls_soliton._bdg_analysis: lam1 = -g1 rho0 / 2).  Below sigma
     and within a factor ~2, as a band-limited seed must be: the
     formation mechanism is the computed BdG instability, not an
     artifact.

  3. THE SHADOW.  The common (q = 0) mode dips at the core while the
     knot grows (common_min 1.000 -> 0.956): the gravitational
     shadow, sourced by formation itself.

  4. THE CHARGE RIDES ALONG, THEN DOES ITS OWN THING.  The seed
     carries the pi_3 winding (B = -0.77 at 48^3. The estimator
     under-reads on coarse grids, knot_charge.py part 1). B is
     tracked through formation and simply REPORTED: it wanders to
     ~0 while the density knot grows, the charge non-protection that
     knot_charge.py registers.  Not re-asserted here.

Usage: python3 skyrmion_3d.py     (~70 s at 48^3, mark slow)
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# the shared 3D machinery (one copy: substrate3d.py)
from substrate3d import (DT, evolve, g0, g1, grid, omega,  # noqa: E402
                         rel_modes, rho0)

SIGMA_BDG = g1 * rho0 / 2.0        # max MI growth rate = |lam1|


def build_seed(N, L, eps=0.3, tex_scale=1.5, env_scale=8.0):
    """The winding texture at small amplitude in the relative sector
    (exactly the review's sim3d seed)."""
    dx, X, Y, Z, R, K2 = grid(N, L)
    f = np.pi * np.exp(-R / tex_scale)
    nx, ny, nz = X / R, Y / R, Z / R
    c, s = np.cos(f), np.sin(f)
    chi1 = c + 1j * s * nz
    chi2 = s * (1j * nx - ny)
    amp = np.exp(-R**2 / env_scale)
    phi1, phi2 = eps * amp * chi1, eps * amp * chi2
    psi = np.empty((3, N, N, N), dtype=complex)
    for k in range(3):
        psi[k] = np.sqrt(rho0) * (1.0 + omega**k * phi1
                                  + omega**(2 * k) * phi2)
    return psi, dx, R, K2


def z3_modes(psi):
    """Field relative modes (shared rel_modes) + density modes
    (upstream _z3_fourier convention): common density and the
    relative density power."""
    phi1, phi2 = rel_modes(psi)
    rho = np.abs(psi)**2
    common = (rho[0] + rho[1] + rho[2]) / 3.0
    q1 = (rho[0] + omega * rho[1] + omega**2 * rho[2]) / 3.0
    q2 = (rho[0] + omega**2 * rho[1] + omega * rho[2]) / 3.0
    relative = np.abs(q1)**2 + np.abs(q2)**2
    return phi1, phi2, common, relative


def topo_charge(phi1, phi2, dx):
    """pi_3 winding of the normalized relative spinor."""
    norm = np.sqrt(np.abs(phi1)**2 + np.abs(phi2)**2)
    m = norm > 0.05 * norm.max()
    nrm = np.maximum(norm, 1e-12)
    n4 = np.stack([phi1.real / nrm, phi1.imag / nrm,
                   phi2.real / nrm, phi2.imag / nrm], axis=-1)
    g = [np.gradient(n4[..., a], dx, axis=ax)
         for ax in range(3) for a in range(4)]
    M = np.stack([n4, np.stack(g[0:4], -1), np.stack(g[4:8], -1),
                  np.stack(g[8:12], -1)], axis=-1)
    dens = np.linalg.det(M) / (2 * np.pi**2)
    return float(np.where(m, dens, 0.0).sum() * dx**3)


def checkpoint(psi, dx, total0, track_B=True):
    phi1, phi2, common, relative = z3_modes(psi)
    return dict(
        B=topo_charge(phi1, phi2, dx) if track_B else None,
        peak=float(relative.max()),
        rms=float(np.sqrt(relative.mean())),
        cmin=float(common.min()),
        cons=abs(float((np.abs(psi)**2).sum()) - total0) / total0)


def formation_run(N, L, g1_run, t_marks, report, tag, track_B=True):
    """track_B=False for the control: the normalized-spinor estimator
    is meaningless on a dispersed field (no localized map to S^3)."""
    psi, dx, R, K2 = build_seed(N, L)
    total0 = float((np.abs(psi)**2).sum())
    t, recs = 0.0, []
    for tm in t_marks:
        n = int(round((tm - t) / DT))
        if n:
            psi = evolve(psi, K2, R, n, g1_eff=g1_run, sponge=False)
        t = tm
        r = checkpoint(psi, dx, total0, track_B)
        r["t"] = t
        recs.append(r)
        b_txt = f"B={r['B']:+.3f}  " if track_B else ""
        report(f"  [{tag}] t={t:5.1f}  {b_txt}"
               f"rel_peak={r['peak']:.4f}  rel_rms={r['rms']:.4f}  "
               f"common_min={r['cmin']:.4f}  cons={r['cons']:.1e}")
    return recs


def growth_rate(recs, t_lo, t_hi):
    """e-folding rate of rel_rms between two checkpoint times."""
    a = next(r for r in recs if r["t"] == t_lo)
    b = next(r for r in recs if r["t"] == t_hi)
    return float(np.log(b["rms"] / a["rms"]) / (t_hi - t_lo))


def run(report=print, N=48, L=24.0):
    report("3D SKYRMION FORMATION: seeded texture vs g1 = 0 control")
    report("=" * 70)
    report(f"  grid {N}^3, box +/-{L/2}, dt = {DT}. g0 = {g0}, "
           f"g1 = g0/sqrt(2) = {g1:.4f}")
    report(f"  BdG: max MI growth rate sigma = g1 rho0/2 = "
           f"{SIGMA_BDG:.4f}")

    marks_true = [0, 1.5, 3, 4.5, 6, 8, 10, 12]
    marks_ctrl = [0, 3, 6, 9]
    tr = formation_run(N, L, g1, marks_true, report, "true")
    ct = formation_run(N, L, 0.0, marks_ctrl, report, "ctrl",
                       track_B=False)

    growth = tr[-1]["rms"] / tr[0]["rms"]          # growth metric
    ctrl_rms = ct[-1]["rms"] / ct[0]["rms"]
    decay = ct[-1]["peak"] / ct[0]["peak"]         # dispersal metric
    sig_meas = growth_rate(tr, 3.0, 8.0)
    report("-" * 70)
    report(f"  formation: rel rms x{growth:.1f} (true, t=12) vs "
           f"x{ctrl_rms:.2f} (ctrl, t=9). Ctrl peak x{decay:.3f} "
           f"(dispersed)")
    report(f"  BdG match: measured e-folding {sig_meas:.3f} vs "
           f"sigma = {SIGMA_BDG:.3f} "
           f"({sig_meas/SIGMA_BDG:.2f} sigma_BdG)")
    report(f"  shadow: common_min {tr[0]['cmin']:.4f} -> "
           f"{tr[-1]['cmin']:.4f}")
    report(f"  charge: B {tr[0]['B']:+.3f} (seed) -> "
           f"{tr[-1]['B']:+.3f} at t=12. The charge is not protected")
    report("  (registered mechanism: knot_charge.py)")

    # ── frozen assertions (measured 2026-07 at 48^3) ─────────────────
    # 1. formation vs control: growth needs the coupling
    assert growth > 5.0, "true coupling must grow the knot (rms)"
    assert ctrl_rms < 1.5, "g1 = 0 has no growth channel (rms flat)"
    assert decay < 0.05, "g1 = 0 seed must disperse (peak collapse)"
    # 2. BdG growth match: below sigma (band-limited seed) and within
    #    a factor ~2 of it (the mechanism IS the computed instability.
    #    measured 0.28 = 0.80 sigma)
    assert 0.5 * SIGMA_BDG < sig_meas <= 1.1 * SIGMA_BDG, \
        f"growth {sig_meas:.3f} inconsistent with BdG sigma {SIGMA_BDG:.3f}"
    # 3. The shadow deepens as the knot forms (measured 0.956)
    assert tr[-1]["cmin"] < 0.99 < tr[0]["cmin"] + 1e-9, \
        "common mode must dip at the core"
    # 4. The seed carries the charge (0.77 at 48^3: coarse-grid
    #    under-read, see knot_charge part 1). Norm is conserved
    assert abs(tr[0]["B"]) > 0.7, "seed must carry the pi_3 charge"
    assert all(r["cons"] < 1e-6 for r in tr), "norm must be conserved"
    return dict(true=tr, ctrl=ct, growth=growth, decay=decay,
                ctrl_rms=ctrl_rms, sigma_measured=sig_meas,
                sigma_bdg=SIGMA_BDG)


if __name__ == "__main__":
    run()
