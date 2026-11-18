"""Ratios that do not depend on the regulator, compared with data.

THE QUESTION. The engine's outputs are ratios: the unit cancels and the
reading holds at every web depth. The gravity lane carries one premise
that is a scale, the covariant regulator t0. It is d11 l_Pl^2 in the
diagnostic, and Sigma/(2 pi) l_Pl^2 with G anchored, the boundary
identification of root.py (sakharov_normalisation.py,
composite_pole.py). Which ratios of the induced action do not depend on
that premise, and how do they compare with data?

PREMISES. One loop of the ledger's rows under one regulator weight
rho(t) on the heat trace proper time, the premise of composite_pole.py
and ledger_kernel.py. Rows light against the regulator scale. The
engine's own couplings, transported by higgs_ew.py.

STATEMENTS (the numbers below check each one).

  R1 (moments). An operator of mass dimension 2n in the induced action
     carries the moment M_{3-n} = int dt t^{n-3} rho(t) e^{-t m^2}: the
     vacuum term M_3, the Einstein term and every scalar mass term M_2,
     the logarithmic terms M_1. Ratios inside one class are free of the
     regulator's scale and of its shape. Ratios across classes keep the
     shape. Checked with four weights (Pauli-Villars at L = 10 m, 30 m,
     proper time at t0 = 1e-3, 1e-4 m^-2): the induced mass term over
     the Einstein stiffness read from the full tensor kernel is 6 per
     real component for every weight to 1e-6, while the vacuum term
     over the stiffness squared, in units of -576 pi^2, is 1.013 and
     1.002 under the proper time cutoff and 2.33 and 2.14 under
     Pauli-Villars, tending to 1 and to 2.107 as the regulator scale
     grows.
  R2 (dimension zero). Lambda G = -(pi/2) str(1)/Sigma^2 times a shape
     factor. Its vanishing needs str(1) = 0. The counted ledger has
     str(1) = 135 (composite_metric.py C6), the one list 128. The
     traceless form removes the term (C2): Lambda is an integration
     constant there.
  R3 (dimension two, the Higgs). The induced coefficient of H^dag H
     over the induced 1/(16 pi G) is C_H/Sigma for every weight, with
     C_H = str d(M^2)/d(H^dag H) over the rows that couple to the
     Higgs: 12 lambda + (9/2) g^2 + (3/2) g'^2 - 12 y_t^2 = 2 C_V, C_V
     the Veltman combination. The data, m_H = 125.2 GeV against the
     engine's M_Pl, need C_H/Sigma = -8 pi m_H^2/M_Pl^2 = -2.6e-33.
     With the engine's couplings the Standard Model rows give
     C_V = -0.218 at the regulator scale M_Pl/sqrt(d11) and -0.206 at
     M_Pl. C_V vanishes only at 2.0e23 GeV. At the regulator scale it
     would need y_t = 0.334 against the engine's 0.385. The induced term
     is tachyonic with |m| = 0.018 M_Pl. A Higgs mass term could cancel
     it only to 31 digits. The engine reads v_EW = M_Pl exp(-S) with no
     such term. The combination is Veltman's, d str M^2 / d(phi^2/2)
     with weights +1 per real scalar, +3 per massive vector, -4 per
     Dirac fermion (Masina and Quiros arXiv:1308.1242). The gauge
     weights hold because every row shares one weight rho(t), the
     common cutoff that Holthausen, Lim and Lindner arXiv:1112.2415
     name as the condition for a scheme free statement. Visser
     gr-qc/0204062 Eqs. 20, 21 put str m^2(H) and 1/G on the same
     moment. Independent transports agree: the central couplings of
     Buttazzo et al. arXiv:1307.3536 give C_V = -0.218 at M_Pl, and
     Hamada, Kawai and Oda arXiv:1210.2538 find the zero near 1e23 GeV
     for m_t = 173 GeV, and at M_Pl only for m_t = 169.8 GeV. Among the
     ledger's other rows only the bridge has an edge to the Higgs
     (root.py WEB["lambda_MPl"]), and that edge is the quartic's. Closing the ratio needs the bridge
     rows to supply C_H = +0.44 at the regulator scale, +0.0024 per
     channel.
  R4 (dimension four). The engine's boundary lambda(M_Pl) = -0.01245 is
     forced by the bridge echo, and the transported beta_lambda at M_Pl
     is small. The logarithmic class of the Higgs boundary is near
     zero. The dimension two class is not.
  R5 (the bridge's quadratic term). The engine carries the bridge's
     Higgs edge for the quartic only. Read as the one loop threshold of
     N real channels with portal kappa, (kappa/2) S^2 H^dag H, under the
     same proper time weight, the quartic is
     delta_lambda = -N kappa^2 M_1/(64 pi^2) and the quadratic term is
     C_H = N kappa, with M_1 = E_1(t0 m_S^2) the logarithmic moment at
     the bridge mass. One kappa gives both. The quartic's sign is the
     engine's, and C_H(bridge) = 8 pi sqrt(N |delta_lambda| / M_1). For
     a bridge mass anywhere from m_e to the regulator scale this lies
     between +3.8 and +81, against the +0.436 the ratio needs. The sign
     is right, and the size is more than eight times too large. Closing needs M_1 = 7.5e3, a
     mass e^(-7.5e3) below the regulator scale.
     For an elementary Higgs no constraint removes the term either.
     Read as an elementary field, the engine's transport coordinate is a
     linear doublet: four real rows in the ledger, and a finite quartic
     boundary carried by the Standard Model running. H^dag H would then
     be a dynamical operator at the regulator scale. A fixed norm field, as
     in the O(N) sigma models where the mass operator is a constant,
     would carry three rows and no quartic.
  R6 (the engine's reading). R3 to R5 exclude a light elementary Higgs
     at the regulator scale, and the engine reads none there. Every light
     scale it derives is the renormalization group invariant scale of
     G2, M_Pl exp(-S) with S = 2 pi/(b0 alpha_G2(M_Pl)) and its vents,
     one loop dimensional transmutation (root.py). The electroweak scale
     is v_EW = M_Pl exp(-S_quark) (masses.py), checked below to the
     engine's own digits. The engine produces no m_H. The lambda(M_Pl)
     transport is a coordinate, not a second physical mass (higgs_ew.py).
     Read as an elementary doublet's boundary and matched to the pole
     mass, it comes out at 125.36 GeV only through a vent that takes m_H
     as data, so m_H is an input on that route
     (tests/probes/higgs_mass.py). The ledger of interference/gravity.py
     counts no elementary Higgs rows at the regulator scale, since the
     condensing channel's constituents are already among its fermion
     rows. If the knots were read as lumps of the carrier, they would
     leave the ledger too. But no carrier mass makes a lump a point
     particle (tests/probes/closed_event_scale.py S5), so the knots stay
     rows (tests/probes/boundary_content.py P4). So the light scales sit
     in the logarithmic class. At the regulator scale the dimension two
     class carries Planckian terms only, and there is nothing to cancel.
     The measured couplings then fix which state of the condensing
     channel the 125 GeV boson is. Its couplings to W, Z and the
     fermions follow the masses that need electroweak breaking, by gauge
     invariance, so kappa_V = kappa_F = 1 for the radial mode of the
     electroweak condensate, a linear composite of the condensing
     channel (kappa_V = 1.035 +/- 0.031, ATLAS, Nature 607, 52 (2022)).
     The other readings fail. An elementary doublet at the regulator
     scale needs C_H to vanish there, and nothing in the ledger makes it
     vanish. A pseudo Goldstone Higgs needs f > 1 TeV (xi < 0.06, arXiv:2407.10950
     Eq. 4.8). No light scalar at all gives kappa_V = 0. And a dilaton of
     the transmuted scale, whose couplings scale as v/f (Goldberger,
     Grinstein and Skiba, PRL 100, 111802 (2008); Bellazzini et al.
     arXiv:1209.3299), is the fluctuation of the scale itself, so it
     moves mu*, M_V and m_t together and couples to gluons with
     kappa_g = +23/2 (R7, tests/probes/scale_ties.py). The radial mode's
     cubic coupling is set by the condensate's potential, which the
     engine does not compute.
  R7 (the gluon coupling). The engine matches G2 to SU(3) at mu* with
     the coset threshold M_V = g v_break sqrt(C2(7)/6), v_break^2 =
     E[v^2] v_EW^2 (couplings.py). By the low energy theorem the 125 GeV
     state couples to gluons through d(1/alpha_s)/d ln v. At beta order 1,
     the order of the threshold matching, the top alone gives
     -1/(3 pi) = -0.106 and the coset threshold 21/(12 pi) = +0.557 per
     unit ln v: if the state moved M_V, kappa_g = -17/4 and gluon fusion
     would run at 18 times the Standard Model rate, against
     kappa_g = 0.92 +/- 0.08 (CMS, Nature 607, 60 (2022)). The engine's
     two loop chain, which carries no two loop decoupling term, gives
     -3.52, with the same verdict. This holds for any state that moves M_V with mu*
     fixed, elementary or composite. A dilaton moves mu* as well and
     gets kappa_g = +23/2 (tests/probes/scale_ties.py). A threshold that
     is not a particle still couples this way when its mass follows a
     field. The radial mode of the electroweak condensate, with G2
     broken by its own field, gives kappa_g = 1.00. The measured
     coupling therefore requires the connector that carries
     v_break^2 = E[v^2] v_EW^2 to transport the vacuum value and not the
     field (couplings.py, THE GLUON CONDITION). The condensate on the G2
     side then has a radial mode of its own, a colour singlet coupling to
     gluons through the threshold at 21/4 times the top's per unit
     ln v_break, with mass 2 alpha_G2 v_break = 34 GeV
     (tests/probes/g2_radial_mode.py).

RESULT. Every ratio of the gravity lane whose data value is a
structural zero (F0/F2 = -2, gamma = 1, c_T = c, the Eotvos parameter)
passes at every order. Among the regulator free ratios that are
compared with a measured number, dimension zero is removed by the
traceless form, and dimension two is the Higgs. The engine derives
v_EW = M_Pl exp(-S). Under one regulator the induced Higgs mass term
must vanish at the regulator scale. The Standard Model rows do not make
it vanish, and the bridge read through its quartic edge overshoots the
other way. The engine places no light elementary scalar there: its
light scales are transmuted and the 125 GeV boson is the radial mode of
the electroweak condensate (R6). So the ratio excludes a reading the
engine does not make.

Usage: python3 tests/probes/regulator_free_ratios.py     (about half a
minute standalone: the engine runs once)
"""

import contextlib
import io
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import couplings, gravity, higgs_ew, masses, mixing, root  # noqa: E402
import composite_metric as cm  # noqa: E402
import composite_pole as cp  # noqa: E402

PV4 = cp.PV4
STEPS = 4000                                     # RK4 steps of the transport


def weights(m2=1.0):
    """Four regulator weights: (J, induced potential V(M^2)) per real component."""
    out = {}
    for lam in (10.0, 30.0):
        J = cp.pauli_villars(m2, lam)[0]

        def V(M2, lam=lam):
            return sum(c * (M2 + s * lam * lam)**2 * math.log(M2 + s * lam * lam)
                       for c, s in PV4) / (64 * math.pi**2)
        out[f"Pauli-Villars L = {lam:g} m"] = (J, V)
    for t0 in (1e-3, 1e-4):
        J = cp.proper_time(m2, t0)[0]

        def V(M2, t0=t0):
            return -expn(3, t0 * M2) / (32 * math.pi**2 * t0**2)
        out[f"proper time t0 = {t0:g}"] = (J, V)
    return out


def moments(m2=1.0, q=1e-2):
    """R1: induced mass term and vacuum term against the kernel's stiffness."""
    rows = {}
    e0 = np.array([q, 0.0, 0.0, 0.0])
    for name, (J, V) in weights(m2).items():
        F2 = cp.form_factors(e0, J)[0]
        stiffness = 4 * F2 / q**2                 # Einstein stiffness per real component
        h = 1e-4 * m2
        mass_term = (V(m2 + h) - V(m2 - h)) / (2 * h)
        rows[name] = dict(mass_over_stiffness=mass_term / stiffness,
                          vacuum_over_stiffness2=V(m2) / stiffness**2 / (-576 * math.pi**2))
    return rows


def engine():
    """The solved engine's dictionaries, run once when standalone."""
    with contextlib.redirect_stdout(io.StringIO()):
        R = root.derive()
        m = masses.derive(R)
        c = couplings.derive(R, m)
        g = gravity.derive(R, m, mixing.derive(R, m), c)
    return m, c, g


def transport_start(m, c, g):
    """(g1, g2, g3, y_t, lambda) at m_t in the engine's transport coordinates."""
    m_t = m["m_t"] / 1e3
    g1 = math.sqrt(4 * math.pi * g["alpha_em_MZ"] / (1 - g["sin2W"])) * math.sqrt(5 / 3)
    g2 = math.sqrt(4 * math.pi * g["alpha_em_MZ"] / g["sin2W"])
    L = math.log(m_t / g["M_Z_derived"])
    y = np.array([higgs_ew._run_gauge_one_loop(g1, 41 / 10, L),
                  higgs_ew._run_gauge_one_loop(g2, -19 / 6, L),
                  math.sqrt(4 * math.pi * c["alpha_s_mt"]),
                  g["yt_mt_transport"], g["lambda_mt_transport"]])
    return m_t, y


def at(mu, start, steps=STEPS):
    m_t, y0 = start
    return higgs_ew._run_rge(y0, 0.0, math.log(mu / m_t), higgs_ew._beta_sm, steps)


def veltman(y):
    """C_V = 6 lambda + (9/4) g^2 + (3/4) g'^2 - 6 y_t^2, g' = sqrt(3/5) g1."""
    g1, g2, _, yt, lam = y
    return 6 * lam + 2.25 * g2**2 + 0.75 * (0.6 * g1**2) - 6 * yt**2


def run(report=print, res=None):
    out = {}
    # R1: the moment classes.
    rows = moments()
    out["mass_over_stiffness"] = {k: r["mass_over_stiffness"] for k, r in rows.items()}
    out["vacuum_over_stiffness2"] = {k: r["vacuum_over_stiffness2"] for k, r in rows.items()}
    for name, r in rows.items():
        report(f"R1 [{name}]: mass term / Einstein stiffness = {r['mass_over_stiffness']:.6f}; "
               f"vacuum term / stiffness^2 = {r['vacuum_over_stiffness2']:.4f} x (-576 pi^2)")

    # R2: dimension zero.
    m, c, g = engine() if res is None else (res["m"], res["c"], res["g"])
    led = cm.ledger()
    sigma = led["sigma"]
    out["str1"], out["sigma"], out["lambda_G"] = led["str1"], sigma, led["lambda_G"]
    report(f"R2 dimension zero: str(1) = {led['str1']}, Lambda G = {led['lambda_G']:+.4f} "
           "(proper time shape). The traceless form removes the term.")

    # R3: dimension two, the Higgs.
    start = transport_start(m, c, g)
    M_Pl = root.M_Pl_GeV
    mu_reg = M_Pl / math.sqrt(root.d11)
    y_reg, y_Pl = at(mu_reg, start), at(M_Pl, start)
    y_check = at(M_Pl, start, steps=5 * STEPS)
    out["transport_check"] = float(np.max(abs(y_Pl - y_check)))
    CV_reg, CV_Pl = veltman(y_reg), veltman(y_Pl)
    mu_zero = math.exp(brentq(lambda L: veltman(at(math.exp(L), start)),
                              math.log(1e21), math.log(1e25), xtol=1e-3))
    g1, g2, _, _, lam = y_reg
    yt_needed = math.sqrt((6 * lam + 2.25 * g2**2 + 0.45 * g1**2) / 6)
    m_H = root.PDG_EW["m_H"]
    data = -8 * math.pi * m_H**2 / M_Pl**2
    induced = 2 * CV_reg / sigma                      # C_H / Sigma
    m_induced = math.sqrt(abs(induced) / (16 * math.pi))
    bridge_need = -2 * CV_reg
    out.update(CV_reg=CV_reg, CV_Pl=CV_Pl, mu_zero=mu_zero, yt_reg=float(y_reg[3]),
               yt_needed=yt_needed, data_ratio=data, induced_ratio=induced,
               m_induced_over_MPl=m_induced, bridge_need=bridge_need,
               bridge_per_channel=bridge_need / root.N_bridge)
    report(f"R3 dimension two: C_V = {CV_reg:+.4f} at mu = M_Pl/sqrt(d11) = {mu_reg:.3e} GeV, "
           f"{CV_Pl:+.4f} at M_Pl. Zero at {mu_zero:.2e} GeV (transport check "
           f"{out['transport_check']:.1e}).")
    report(f"   C_H/Sigma = {induced:+.3e} induced against {data:+.3e} from m_H = {m_H} GeV: "
           f"the induced term is tachyonic, |m| = {m_induced:.4f} M_Pl. At the regulator scale "
           f"C_V = 0 needs y_t = {yt_needed:.4f} against {y_reg[3]:.4f}. The bridge rows would "
           f"need C_H = {bridge_need:+.4f}, {out['bridge_per_channel']:+.5f} per channel.")

    # R4: dimension four.
    beta = higgs_ew._beta_sm(y_Pl)
    out["lambda_Pl"], out["beta_lambda_Pl"] = float(y_Pl[4]), float(beta[4])
    report(f"R4 dimension four: lambda(M_Pl) = {y_Pl[4]:+.5f} (forced boundary "
           f"{g['lam_Pl_boundary']:+.5f}), d lambda/d ln mu = {beta[4]:+.2e} there.")
    # R5: the bridge's quadratic term from its quartic edge.
    N, delta = root.N_bridge, -g["lam_Pl_boundary"]
    scales = {"m_e": (root.M_E_ANCHOR_MEV / 1e3 / mu_reg)**2, "100 GeV": (100.0 / mu_reg)**2,
              "1e-2 of the regulator": 1e-2, "the regulator scale": 1.0}
    bridge = {k: 8 * math.pi * math.sqrt(N * delta / expn(1, x)) for k, x in scales.items()}
    M1_needed = (8 * math.pi)**2 * N * delta / bridge_need**2
    out.update(bridge_portal=bridge, M1_needed=M1_needed)
    report("R5 bridge read as a one loop portal through its quartic edge: C_H(bridge) = "
           + ", ".join(f"{v:+.2f} ({k})" for k, v in bridge.items())
           + f" against {bridge_need:+.3f}. Closing needs M_1 = {M1_needed:.0f}.")
    # R6: the engine's light scales are transmuted.
    S_quark = root.WEB.state["S_quark"]
    out["v_over_MPl"] = m["v_EW_pred_GeV"] / M_Pl
    out["transmuted"] = math.exp(-S_quark)
    report(f"R6 v_EW/M_Pl = {out['v_over_MPl']:.6e} = exp(-S_quark) = {out['transmuted']:.6e}: "
           "the light scale is the transmuted one, a logarithmic reading.")
    # R7: the gluon coupling through the G2 -> SU(3) threshold.
    from interference import couplings as cpl
    mu_star = M_Pl * math.exp(-(root.S_lepton - float(root.C2_26)))
    a_star = math.pi / 32
    g_G2 = math.sqrt(4 * math.pi * a_star)
    v0, mt0 = m["v_EW_pred_GeV"], m["m_t"] / 1e3

    def run1(a0, mu0, mu1, nf):                     # beta order 1
        return 1.0 / (1.0 / a0 + cpl._b0(nf) * math.log(mu1 / mu0))

    def inv_alpha(v, track_MV, track_mt, mu=125.0, run=run1):
        MV = g_G2 * (v if track_MV else v0) * math.sqrt(float(root.E_v2) / 3)
        mt = mt0 * (v / v0 if track_mt else 1.0)
        lam3 = 1.0 - root.d11 * root.n7 * math.log(MV / mu_star)
        a = 1.0 / (1.0 / a_star - lam3 / (12 * math.pi))
        return 1.0 / run(run(a, mu_star, mt, 6), mt, mu, 5)

    def slope(track_MV, track_mt, h=1e-4):
        return (inv_alpha(v0 * math.exp(h), track_MV, track_mt)
                - inv_alpha(v0 * math.exp(-h), track_MV, track_mt)) / (2 * h)
    top, both = slope(False, True), slope(True, True)
    out.update(kappa_g_top=top, kappa_g_both=both, kappa_g_tracking=both / top,
               alpha_s_MZ_chain=1.0 / inv_alpha(v0, True, True, mu=91.1876,
                                                run=cpl._run_SM_2loop))
    report(f"R7 d(1/alpha_s)/d ln v at 125 GeV, beta order 1: top {top:+.4f}, top and coset threshold {both:+.4f}: "
           f"kappa_g = {both / top:+.3f} if the 125 GeV state moved M_V (measured 0.92 +/- 0.08).")
    report("RESULT: the Einstein term and every scalar mass term share one moment of the "
           "regulator, so their ratio is free of it. For the Higgs that ratio must vanish at "
           "the regulator scale, and the Standard Model rows with the engine's couplings give "
           "C_V = -0.22 there. Dimension zero is removed by the traceless form. Dimension two "
           "is the Higgs: the bridge read through its quartic edge has the right sign and "
           "is more than eight times too large.")
    return out


if __name__ == "__main__":
    run()
