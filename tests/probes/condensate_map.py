"""The G2-breaking condensate map, executed.

The alpha_s chain needs one scale the WZW identity does not supply:
the G2 -> SU(3) breaking condensate v7 in M_V = g v7/sqrt(3).  The
coded value is v7^2 = E[v^2] v_EW^2 with E[v^2] = B^2/4 = 1/2 the
bridge retained second moment.  This probe executes
the field-map question behind that value and asserts three facts.

FACT 1 (the component identification is not a field map).  The engine
consumes v_EW in the real-canonical norm: G_F = 1/(sqrt(2) v_EW^2),
A0 = sqrt(pi alpha) v_EW, y_t = sqrt(2) m_t/v_EW.  Built here
explicitly: the SM doublet's four real components under SU(2)_L obey
the broken-generator sum rule

    Sum_a M_a^2 = g^2 C2(2) |h|^2,   C2(2) = 3/4,  |h| = v_EW,

with three equal masses M_W = g v_EW/2, in BOTH bookkeepings of the
same doublet (real-canonical and complex).  The 7 of G2 is a real
representation with the same canonical norm and the same trace
arithmetic C2 = dim(G) T(R)/dim(R), so a LINEAR identification of the
7-singlet with the condensing direction transports the real-norm
length: v7 = v_EW, M_V = 157.9 GeV, alpha_s(M_Z) = 0.1157, -2.5
sigma.  Reading v7 as "the neutral component of the complex doublet"
takes the complex-component norm of a real field: the factor
1/sqrt(2) is structural, not notational.

FACT 2 (the factor is the connector's per-mode contraction).  The
crystallized Z3 harmonic v_k = B cos(theta + 2 pi k/3) has exactly
two quadratic contractions (root.py):

    per-mode  E[v^2] = |v_hat_1|^2 = B^2/4 = 1/2   (what the bridge sees)
    total     P v^2 P = B^2/2 = 1                  (both channels)

The engine consumes the per-mode contraction at every quadratic
bridge reading (xi_bridge, delta_bridge).  Applied to the transported
condensate second moment it gives v7^2 = E[v^2] v_EW^2, which
reproduces the coded M_V = g v_EW/sqrt(6) = 111.648 GeV exactly, and
alpha_s(M_Z) separates the two contractions by 3.0 sigma: per-mode
+0.43, total -2.53.

FACT 3 (the scale menu closes on the electroweak condensate).  Every
scale in the theory is M_Pl e^{-S} for a derived action.  Scanning
the ladder against alpha_s(M_Z): the G2-face action scale M_R =
M_Pl e^{-2 pi^2} (the seesaw scale, the would-be independent G2-side
condensate) fails by -72 sigma, Lambda_conf by +96 sigma, and no
threshold at all by -6.7 sigma.  Only the electroweak row survives
(mu* = v_EW e^{15/512} sits 3 percent away and is numerically
admissible, but it is a gauge-matching endpoint, not a condensate).

FACT 4 (the reading closes: the closed-world sector argument).
Three lemmas discharge the reading assignment, the closed-world
genre of the EM-channel finality argument in root.py; the residual
is named below: the per-mode transfer to a
transported condensate is the engine's uniform contraction grammar,
not an independent transport computation, and the lemma chain uses
one empirical row-kill and one role exclusion.
  (a) NO INDEPENDENT SCALE: the derived actions each carry a
      sector already: S_lepton (the lepton knot), S_quark (the F4
      knot, v_EW), S_baryo = 2 pi^2 (the G2 face, M_R = M_Pl
      e^{-S_baryo}, consumed by the seesaw).  The 7-scalar carries
      none, so its condensate can only be a reading of the one
      condensing channel, and the unique route is the (7,26)
      connector.  Asserted: M_R = M_Pl e^{-S_baryo} exactly,
      and its alpha_s row is dead (FACT 3).
  (b) THE COEFFICIENT IS NOT FREE: the connector's normalized Schur
      trace over all 182 copies introduces no independent bridge
      amplitude, so the quadratic reading is E[v^2] with nothing to
      tune.
      Asserted: the cross-pin of the two independent structures,
      the octonionic CG amplitude B = sqrt(d10) (Fano cross
      product) against the MTC Schur trace, B^2/4 = E[v^2] exactly.
  (c) NO VENT DEFLATION: the (1 - h10) share law fires only at a
      confining interface, and the matching sits three decades
      above it.  Asserted: Lambda_conf << M_V < mu*.

LOOK-ELSEWHERE, stated: at the PDG width, five engine atoms fall
within 1 sigma as bare factors f^2 in v7^2 = f^2 v_EW^2 (1/2, 4/9,
3/5, 2/3, 9/13), so alpha_s corroborates the per-mode reading and
cannot alone pin it.  The pin is the structure above; the
neighbours separate from 1/2 by > 5 sigma at the watch-1 lattice
width (~0.03 percent), the kill.

Usage: python3 condensate_map.py
"""

import contextlib
import io
import math
import os
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
for _path in (os.path.dirname(_here),
              os.path.join(os.path.dirname(os.path.dirname(_here)),
                           "interference"),
              os.path.join(_here, "interference")):
    if _path not in sys.path:
        sys.path.insert(0, _path)
import _bootstrap  # noqa: E402,F401

from root import (d10, d11, n7, C2_7, C2_26, E_v2, dim_G2,   # noqa: E402
                  M_Pl_GeV, S_lepton, S_baryo, BA_ratio,
                  rk4_run, PDG_EW, WEB)

SIGMA_PDG = 0.0009


# ═══════════════════════════════════════════════════════════════════════
#  FACT 1: the SM broken-generator sum rule in both bookkeepings
# ═══════════════════════════════════════════════════════════════════════

_PAULI = [np.array([[0, 1], [1, 0]], dtype=complex),
          np.array([[0, -1j], [1j, 0]], dtype=complex),
          np.array([[1, 0], [0, -1]], dtype=complex)]


def _real_embedding(m_complex):
    """C^2 -> R^4 embedding of a complex 2x2 matrix, coordinates
    (Re H+, Im H+, Re H0, Im H0)."""
    r = np.zeros((4, 4))
    for i in range(2):
        for j in range(2):
            a, b = m_complex[i, j].real, m_complex[i, j].imag
            r[2 * i, 2 * j], r[2 * i, 2 * j + 1] = a, -b
            r[2 * i + 1, 2 * j], r[2 * i + 1, 2 * j + 1] = b, a
    return r


def sm_sum_rule(v_ew):
    """Vector masses of the broken SU(2)_L from the doublet VEV, in the
    real-canonical and the complex bookkeeping of the same field."""
    # real-canonical: h = sqrt(2) r, |h| = v_EW, kinetic (1/2)(dh)^2,
    # M^2_ab = g^2 (K^a h).(K^b h) with K^a the real embedding of
    # -i sigma^a/2
    K = [_real_embedding(-1j * 0.5 * s) for s in _PAULI]
    h = np.array([0.0, 0.0, v_ew, 0.0])
    masses_real = [float(np.dot(K[a] @ h, K[a] @ h)) for a in range(3)]
    # complex: H = (0, v/sqrt(2)), kinetic |dH|^2,
    # M^2_ab = 2 g^2 Re[<H>* T^a T^b <H>]
    H = np.array([0.0, v_ew / math.sqrt(2.0)], dtype=complex)
    masses_cplx = [float(2.0 * np.real(
        np.vdot(0.5 * s @ H, 0.5 * s @ H))) for s in _PAULI]
    casimir = 0.75  # C2(2) = dim(SU2) T(2)/dim(2) = 3 (1/2)/2
    return {
        "masses_real_over_g2": masses_real,
        "masses_complex_over_g2": masses_cplx,
        "casimir_sum_over_g2v2": sum(masses_real) / v_ew**2,
        "casimir_C2_doublet": casimir,
        "mw_over_gv": math.sqrt(masses_real[0]) / v_ew,
        "real_norm_length": float(np.linalg.norm(h)),
        "complex_norm_length": float(np.linalg.norm(H)),
    }


# ═══════════════════════════════════════════════════════════════════════
#  The alpha_s chain (couplings.py, with v7 as the scanned coordinate)
# ═══════════════════════════════════════════════════════════════════════

def _b0(nf):
    return (11 * d11 - 2 * nf) / (6.0 * math.pi)


def _b1(nf):
    c_a = d11
    c_f = (d11**2 - 1) / (2.0 * d11)
    return (34 * c_a**2 - (10 * c_a + 6 * c_f) * nf) / (24.0 * math.pi**2)


def _run2(a0, mu0, mu1, nf):
    return rk4_run(lambda a: -_b0(nf) * a**2 - _b1(nf) * a**3,
                   a0, math.log(mu0), math.log(mu1), 10000)


def alpha_s_of_v7(v7, mu_star, m_t, m_z):
    """M_V = g v7 sqrt(C2(7)/6), the Weinberg-Hall threshold at mu*,
    then the two-loop SM transport to M_Z."""
    a_mu = math.pi / 2**5
    g = math.sqrt(4.0 * math.pi * a_mu)
    n_broken = 2 * n7 - (d11**2 - 1)
    m_v = g * v7 * math.sqrt(float(C2_7) / n_broken)
    lam3 = 1.0 - float(d11 * n7) * math.log(m_v / mu_star)
    inv = 1.0 / a_mu - lam3 / (12.0 * math.pi)
    if inv <= 0.0 or 1.0 / inv > 1.0:
        return m_v, float("nan")
    a_mt = _run2(1.0 / inv, mu_star, m_t, 6)
    return m_v, _run2(a_mt, m_t, m_z, 5)


# ═══════════════════════════════════════════════════════════════════════
#  The probe
# ═══════════════════════════════════════════════════════════════════════

def run(report=print, engine_m_t_GeV=None):
    if engine_m_t_GeV is None:
        import masses as _masses
        import root as _root
        with contextlib.redirect_stdout(io.StringIO()):
            _R = _root.derive()
            _M = _masses.derive(_R)
        engine_m_t_GeV = _M["m_t"] / 1e3

    v_ew = M_Pl_GeV * math.exp(-WEB.state["S_quark"])
    mu_star = M_Pl_GeV * math.exp(-(S_lepton - float(C2_26)))
    m_z = PDG_EW["M_Z"]
    a_pdg = PDG_EW["alpha_s_MZ"]

    def pull(a):
        return (a - a_pdg) / SIGMA_PDG

    # ── FACT 1 ──
    sm = sm_sum_rule(v_ew)
    assert abs(sm["casimir_sum_over_g2v2"] - 0.75) < 1e-12
    assert max(sm["masses_real_over_g2"]) - min(
        sm["masses_real_over_g2"]) < 1e-9
    assert np.allclose(sm["masses_real_over_g2"],
                       sm["masses_complex_over_g2"], atol=1e-9)
    assert abs(sm["mw_over_gv"] - 0.5) < 1e-12
    assert abs(sm["real_norm_length"] / sm["complex_norm_length"]
               - math.sqrt(2.0)) < 1e-12
    # the same trace arithmetic on the G2 side: C2(7) = 14/7 = 2
    assert C2_7 == dim_G2 * 1 // n7 == d10

    report("FACT 1: SM sum rule Sum M^2 = g^2 (3/4) v_EW^2, three equal")
    report(f"  masses M_W = g v_EW/2, identical in both bookkeepings;")
    report(f"  real/complex norm ratio sqrt(2). C2(7) = 14/7 = {C2_7}.")

    # ── FACT 2: the two contractions ──
    mv_half, a_half = alpha_s_of_v7(
        v_ew * math.sqrt(float(E_v2)), mu_star, engine_m_t_GeV, m_z)
    mv_one, a_one = alpha_s_of_v7(v_ew, mu_star, engine_m_t_GeV, m_z)
    separation = (a_half - a_one) / SIGMA_PDG
    report(f"FACT 2: per-mode E[v^2] = {float(E_v2)}: "
           f"M_V = {mv_half:.3f}, alpha_s = {a_half:.5f} "
           f"({pull(a_half):+.2f} sigma)")
    report(f"        total 1:            M_V = {mv_one:.3f}, "
           f"alpha_s = {a_one:.5f} ({pull(a_one):+.2f} sigma)")
    report(f"        dichotomy separation {separation:.2f} sigma")

    # ── FACT 3: the scale menu ──
    lambda_conf = (1.0 / d10) * M_Pl_GeV * math.exp(-S_lepton)
    m_r = M_Pl_GeV * math.exp(-2.0 * math.pi**2)
    scale_menu = {}
    for sname, s in (("v_EW", v_ew), ("mu_star", mu_star),
                     ("M_R", m_r), ("Lambda_conf", lambda_conf)):
        for cname, c in (("per_mode", float(E_v2)), ("total", 1.0)):
            _, a = alpha_s_of_v7(s * math.sqrt(c), mu_star,
                                 engine_m_t_GeV, m_z)
            scale_menu[f"{sname}:{cname}"] = pull(a)
    a_mt = _run2(math.pi / 2**5, mu_star, engine_m_t_GeV, 6)
    a_none = _run2(a_mt, engine_m_t_GeV, m_z, 5)
    report("FACT 3: scale menu pulls (sigma): "
           + ", ".join(f"{k} {v:+.1f}" for k, v in scale_menu.items()))
    report(f"        no threshold: {pull(a_none):+.1f} sigma")

    # ── FACT 4: the closed-world lemmas ──
    # (a) the G2 face's own action scale is the seesaw scale, exactly
    assert abs(m_r / (M_Pl_GeV * math.exp(-S_baryo)) - 1.0) < 1e-12
    assert abs(scale_menu["M_R:per_mode"]) > 50
    # (b) the coefficient cross-pin: octonionic CG (B = sqrt(d10),
    #     Fano cross product) against the MTC Schur trace, exactly
    assert abs(BA_ratio**2 / 4.0 - float(E_v2)) < 1e-15
    # (c) no confining interface at the matching scale
    assert lambda_conf * 100.0 < mv_half < mu_star
    report("FACT 4: M_R = M_Pl e^(-S_baryo) (the G2-face action is the")
    report(f"  seesaw's, dead row); B^2/4 = E[v^2] cross-pin exact;")
    report(f"  Lambda_conf ({lambda_conf:.3f}) << M_V ({mv_half:.1f}) "
           f"< mu* ({mu_star:.1f}): no vent fires.")

    # ── look-elsewhere over the weight-class atoms ──
    atoms = {"E_v2 = 1/2": 0.5, "h10 = 2/9": 2.0 / 9.0,
             "1 - h10 = 7/9": 7.0 / 9.0, "Q0 = 2/3": 2.0 / 3.0,
             "Q0^2 = 4/9": 4.0 / 9.0, "h_7 = 2/5": 0.4,
             "h_26 = 3/5": 0.6, "1/d11 = 1/3": 1.0 / 3.0,
             "4/13": 4.0 / 13.0, "9/13": 9.0 / 13.0,
             "3/13": 3.0 / 13.0, "total = 1": 1.0}
    atom_pulls = {}
    for name, f2 in atoms.items():
        _, a = alpha_s_of_v7(v_ew * math.sqrt(f2), mu_star,
                             engine_m_t_GeV, m_z)
        atom_pulls[name] = pull(a)
    within = [k for k, v in atom_pulls.items() if abs(v) <= 1.0]
    report(f"  look-elsewhere: {len(within)} of {len(atoms)} atoms "
           f"within 1 sigma: {sorted(within)}")

    return {
        "sm_real_norm": sm,
        "closed_world": {
            "M_R_is_G2_face_action_scale": True,
            "coefficient_cross_pin": float(BA_ratio**2 / 4.0),
            "Lambda_conf_GeV": lambda_conf,
            "no_interface_at_matching": lambda_conf * 100.0 < mv_half,
        },
        "per_mode": {"M_V": mv_half, "alpha_s_MZ": a_half,
                     "pull": pull(a_half)},
        "linear_map": {"M_V": mv_one, "alpha_s_MZ": a_one,
                       "pull": pull(a_one)},
        "dichotomy_separation_sigma": separation,
        "scale_menu_pulls": scale_menu,
        "no_threshold_pull": pull(a_none),
        "atom_pulls": atom_pulls,
        "atoms_within_1sigma": sorted(within),
        "v_EW_GeV": v_ew,
    }


if __name__ == "__main__":
    run()
