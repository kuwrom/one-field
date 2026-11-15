"""The Higgs pole mass the quartic boundary implies, used as data on the vent.

THE QUESTION. The engine fixes the Higgs quartic at the Planck scale by
the bridge echo, lambda(M_Pl) = -N_bridge alpha_G2^2 E[v^2] (1 - h10)
= -0.01245 (root.py WEB["lambda_MPl"]), and transports it to the top
scale with the Standard Model equations (higgs_ew.py), reaching
lambda(m_t) = 0.126508. What pole mass does that quartic mean?

PREMISES. The Higgs read as an elementary doublet from M_Pl down, so
the Standard Model transport applies over the whole interval (the
alternative to the engine's transmuted reading, regulator_free_ratios.py
R6). The engine's
couplings at m_t (higgs_ew.py). The NNLO MS-bar to pole matching of
Buttazzo et al., arXiv:1307.3536v2, eq. 55:
    lambda(M_t) = 0.12711 + 0.00206 (M_h/GeV - 125.66)
                  - 0.00004 (M_t/GeV - 173.10) +- 0.00030 (theory).

STATEMENTS (the numbers below check each one).
  H1  At tree level sqrt(2 lambda(m_t)) v = 123.85 GeV.
  H2  With the NNLO matching, M_h = 125.36 +- 0.15 GeV at the engine's
      M_t = 172.508 GeV, against 125.20 +- 0.11 (PDG): +0.85 sigma.
  H3  The boundary does the work: lambda(M_Pl) = 0 gives 129.39 GeV,
      and M_h moves by 0.35 GeV per 0.001 in lambda(M_Pl).  Matching
      within 0.2 GeV therefore needs the boundary within 5 percent.
  H4  The bridge coupling enters the boundary squared. Taken as the
      MS-bar coupling at M_Pl (two loop running up from pi/32 at mu*,
      1.5 percent below the action value 1/(24 pi)), the boundary gives
      M_h = 125.48 GeV (+1.6 sigma).

  H5  The vent factor (1 - h10) in the boundary is one entry of a
      21-term menu, fixed by the measured Higgs mass. In this transport
      the menu has two readings within 1 sigma: 1/6 gives 125.04 GeV
      (-0.85 sigma) and h10 = 2/9 gives 125.36 GeV (+0.86 sigma). No
      vent gives 124.09 GeV.

RESULT. The boundary matches the Higgs mass. The measured mass
constrains the vent and leaves two readings, 1/6 and h10.  The pole
matching is the extra stage higgs_ew.py excludes ("identified at the
depth-4 Web endpoint rather than by an extra matching stage"). The
Higgs mass enters the chain as data on the vent, and the engine's
outputs contain no M_h.

Usage: python3 tests/probes/higgs_mass.py     (about a minute: the
engine runs once and the transport is solved again five times)
"""

import contextlib
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import couplings, gravity, higgs_ew, masses, mixing, root  # noqa: E402

M_H_PDG = (root.PDG_EW["m_H"], 0.11)
TH = 0.00030 / 0.00206                     # GeV, the matching's stated theory error


def pole_mass(lam_mt, M_t):
    return 125.66 + (lam_mt - 0.12711 + 0.00004 * (M_t - 173.10)) / 0.00206


def transport():
    with contextlib.redirect_stdout(io.StringIO()):
        R = root.derive()
        m = masses.derive(R)
        c = couplings.derive(R, m)
        g = gravity.derive(R, m, mixing.derive(R, m), c)
    mt = m["m_t"] / 1e3
    lmz = math.log(mt / g["M_Z_derived"])
    g1 = higgs_ew._run_gauge_one_loop(
        math.sqrt(4 * math.pi * g["alpha_em_MZ"] / (1 - g["sin2W"])) * math.sqrt(5 / 3), 41 / 10, lmz)
    g2 = higgs_ew._run_gauge_one_loop(math.sqrt(4 * math.pi * g["alpha_em_MZ"] / g["sin2W"]), -19 / 6, lmz)
    g3 = math.sqrt(4 * math.pi * c["alpha_s_mt"])
    yt, t_pl = g["yt_mt_transport"], math.log(root.M_Pl_GeV / mt)

    def lam_mt(lam_pl, n=8000):
        li = 0.0
        for _ in range(3):
            up = higgs_ew._run_rge(np.array([g1, g2, g3, yt, li]), 0.0, t_pl, higgs_ew._beta_sm, n)
            li = higgs_ew._run_rge(np.array([up[0], up[1], up[2], up[3], lam_pl]),
                                   t_pl, 0.0, higgs_ew._beta_sm, n)[4]
        return li
    return g, m, mt, lam_mt


def run(report=print):
    g, m, mt, lam_mt = transport()
    lam_pl = g["lam_Pl_boundary"]
    lam0 = lam_mt(lam_pl)
    tree = math.sqrt(2 * lam0) * m["v_EW_pred_GeV"]
    Mh = pole_mass(lam0, mt)
    pull = (Mh - M_H_PDG[0]) / math.hypot(TH, M_H_PDG[1])
    scan = {lp: pole_mass(lam_mt(lp), mt) for lp in (0.0, -0.006, -0.02)}
    slope = (pole_mass(lam_mt(lam_pl + 1e-3), mt) - Mh) / 1e-3
    b0, b1 = 32 / 3, 232 / 3
    mu_star = root.M_Pl_GeV * math.exp(-(root.S_lepton - float(root.C2_26)))
    L, a, n = math.log(root.M_Pl_GeV / mu_star), math.pi / 32, 20000
    f = lambda x: -(b0 / (2 * math.pi)) * x * x - (b1 / (8 * math.pi**2)) * x**3
    h = L / n
    for _ in range(n):
        k1 = f(a); k2 = f(a + h * k1 / 2); k3 = f(a + h * k2 / 2); k4 = f(a + h * k3)
        a += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    Mh_ms = pole_mass(lam_mt(lam_pl * (a * 24 * math.pi)**2), mt)
    base = lam_pl / (1 - float(root.h10))
    vents = {"none": 0.0, "1/6": 1 / 6, "h10": float(root.h10), "1/4": 0.25}
    vent_mh = {k: pole_mass(lam_mt(base * (1 - x)), mt) for k, x in vents.items()}
    sig = math.hypot(TH, M_H_PDG[1])
    landings = sorted(k for k, v in vent_mh.items() if abs(v - M_H_PDG[0]) < sig)
    assert landings == ["1/6", "h10"]
    assert abs(lam0 - 0.126508) < 2e-6 and abs(tree - 123.85) < 0.02
    assert abs(Mh - 125.36) < 0.02 and abs(pull - 0.85) < 0.05
    assert abs(scan[0.0] - 129.39) < 0.05 and 300 < slope < 360
    report(f"H1 lambda(M_Pl) = {lam_pl:+.5f} -> lambda(m_t) = {lam0:.6f}. Tree sqrt(2 lambda) v = {tree:.2f} GeV")
    report(f"H2 NNLO pole matching at M_t = {mt:.3f} GeV: M_h = {Mh:.2f} +- {TH:.2f} GeV, "
           f"PDG {M_H_PDG[0]} +- {M_H_PDG[1]}: {pull:+.2f} sigma")
    report("H3 boundary scan: " + ", ".join(f"lambda(M_Pl) = {k:+.3f} -> {v:.2f} GeV" for k, v in scan.items())
           + f"; dM_h/dlambda(M_Pl) = {slope:.0f} GeV")
    report(f"H4 MS-bar bridge coupling at M_Pl ({100 * (a * 24 * math.pi - 1):+.1f}%): M_h = {Mh_ms:.2f} GeV")
    report("H5 vent menu: " + ", ".join(f"{k} {v:.2f} GeV" for k, v in vent_mh.items())
           + f". Within 1 sigma: {landings}. The measured mass constrains the vent.")
    return dict(vent_mh=vent_mh, landings=landings, lam_pl=lam_pl, lam_mt=lam0, tree=tree, M_h=Mh, pull=pull, scan=scan,
                slope=slope, M_h_msbar=Mh_ms, M_t=mt)


if __name__ == "__main__":
    run()
