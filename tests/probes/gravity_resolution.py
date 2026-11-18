"""How sharp gravity is, channel by channel, and where its error comes from.

Gravity has two readings in the engine.  The anchor obtains the
boundary M_b from the electron through the lepton transmutation and its
echo ladder (root.py).  The full Planck identification then reads M_b
as sqrt(hbar c/G).  The ledger reads the stiffness Sigma of the one
list (boundary.py, gravity.py).  This probe states the resolution of
each, attributes the anchor's residual, and checks which readings are
true readings and which are identities.

R1 (the ladder).  M_Pl against CODATA at web depth 1, 2 and 4:
   -1140, +21.5, +22.2 ppm.  Each layer is suppressed by a further
   power of alpha/2 pi.
R2 (attribution).  The anchor's inputs are m_e and alpha (measured)
   and S_lepton, C_EFF and the harmonic's Delta_e (exact).  The
   measured inputs move M_Pl by parts in 1e10 and 1e12.  The next
   echo of the same ladder is of order (alpha/2 pi)^3, near a part in
   1e9.  The +22 ppm residual is therefore neither input noise nor a
   missing echo of the ladder unless that echo's coefficient exceeds
   1e4.  It sits at the scale of G's own measurement, two standard
   deviations (CODATA 22 ppm on G, 11 ppm on M_Pl), or it comes from
   an identification. The probe does not decide which.  Written as an
   action shift it is delta S = 2.2e-5.
R3 (readings and identities).  A reading moves when an independent
   input moves.  The anchor moves with m_e one for one.  The ledger's
   G at the anchored regulator t0 = Sigma/(2 pi) M_b^-2 does not move
   when Sigma does: an identity.  The ledger's G at a regulator placed
   independently moves as 1/Sigma: a reading, but only of the ledger's
   coarse estimate (a fixed sign and a range of about two decades).
R4 (the coarse estimate).  The ledger's own resolution: Sigma over the
   vector contacts and the composite metric's own loop, and G over the
   Planckian band of regulator placements.  The sign is fixed in every
   cell. The magnitude spans about two decades.
R5 (where the frame is read).  The stiffness nodes of the web carry
   vent edges and one boundary reading (xi on the bridge row), all at
   web depth one.  The bridge carries the EM loop and receives no echo
   of it.  The common mode's own node reads at depth four.  It has two
   edges, the relative powers of the lepton circulant on its depth four
   rates, and their parent is the closed e<->q carrier (root.py, THE
   COMMON MODE).  The anchor node M_Pl_MeV sits at web depth four, read
   through the lepton loop.

RESULT. Computed: the anchor puts M_Pl at -1139.9, +21.5 and +22.2 ppm
from CODATA at web depth 1, 2 and 4, and the +22.2 ppm residual is 2.02
sigma of M_Pl's 11 ppm uncertainty (half of G's 22 ppm), far above the
moves of the measured inputs and the next echo of the ladder. The anchor
moves one for one with m_e, a reading, while the ledger's G at the
anchored regulator is an identity. At a placed regulator the ledger
gives an attractive sign in every cell over a G range of about two
decades.

Usage: python3 tests/probes/gravity_resolution.py     (about 20 s, runs the chain)
"""

import contextlib
import io
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

CODATA_M = r.M_PL_CODATA_GEV
SIGMA_G_REL = 22e-6                     # CODATA relative uncertainty on G
REL_M_E = 3.1e-10                       # CODATA relative uncertainty on m_e
INV_ALPHA_CS, INV_ALPHA_RB = 137.035999046, 137.035999206
ONE = B.rows(r.d10, r.d11, r.n7, r.n26)       # the one list: the bridge as amount (boundary.py)
CONTACT_SIGMA = tuple(B.sigma(ONE, 0.0, c, float(r.xi_bridge)) for c in B.CONTACTS)   # contacts +1/0/-1
assert all(abs(s - v) < 1e-9 for s, v in zip(CONTACT_SIGMA, (155 / 6, 263 / 6, 371 / 6)))
METRIC_LOOP = (0.0, -8.0 / 3.0, -23.0 / 3.0)  # none, on shell, de Donder off shell
BAND = (0.6, 15.0)                            # gate M: Planckian t0 M_b^2


def m_pl(m_e=None, inv_alpha=None, S=None, depth=4):
    """The anchor inversion, in GeV, through web depth 1, 2 or 4."""
    m_e = r.M_E_ANCHOR_MEV if m_e is None else m_e
    a = 1.0 / (r.inv_alpha_phys if inv_alpha is None else inv_alpha)
    e = a / (2 * math.pi)
    S = r.S_lepton if S is None else S
    De = r._Delta_e
    qed = 1.0 - e if depth >= 2 else 1.0
    vent = 1.0 - (r.d11**3 / r.d10) * e * e * De if depth >= 4 else 1.0
    return m_e / (r.C_EFF * math.exp(-S) * De**2 * qed * vent) / 1e3


def ppm(M):
    return 1e6 * (M / CODATA_M - 1.0)


def ladder():
    out = {d: ppm(m_pl(depth=d)) for d in (1, 2, 4)}
    assert abs(m_pl() / r.M_Pl_GeV - 1.0) < 1e-14, "the probe's inversion is the engine's"
    return out


def attribution():
    M0 = m_pl()
    e = 1.0 / (2 * math.pi * r.inv_alpha_phys)
    from_m_e = REL_M_E                                   # M is linear in m_e
    from_alpha = abs(m_pl(inv_alpha=INV_ALPHA_RB) / m_pl(inv_alpha=INV_ALPHA_CS) - 1.0)
    next_echo = e**3                                     # order of a depth 5 term
    step_4 = abs(m_pl(depth=4) / m_pl(depth=2) - 1.0)
    resid = M0 / CODATA_M - 1.0
    dS = math.log(M0 / CODATA_M)                         # M ~ e^S, so delta S = ln ratio
    sigma_M = SIGMA_G_REL / 2.0
    return dict(from_m_e=from_m_e, from_alpha=from_alpha, next_echo=next_echo,
                step_4=step_4, residual=resid, delta_S=dS,
                residual_over_sigma=resid / sigma_M,
                residual_over_next_echo=resid / next_echo)


def readings():
    """R3: perturb an independent input; a reading moves, an identity does not."""
    eps = 1e-6
    anchor = m_pl(m_e=r.M_E_ANCHOR_MEV * (1 + eps)) / m_pl() - 1.0
    G_anchored = lambda s: 2 * math.pi * (s / (2 * math.pi)) / s
    G_placed = lambda s, t0=1.0: 2 * math.pi * t0 / s
    s0 = CONTACT_SIGMA[0]
    return dict(
        anchor_response=anchor / eps,
        ledger_anchored_response=(G_anchored(2 * s0) / G_anchored(s0) - 1.0),
        ledger_placed_response=(G_placed(2 * s0) / G_placed(s0) - 1.0))


def silhouette():
    sig = [s + d for s in CONTACT_SIGMA for d in METRIC_LOOP]
    G = [2 * math.pi * t / s for s in sig for t in BAND]
    return dict(sigma_min=min(sig), sigma_max=max(sig),
                G_min=min(G), G_max=max(G),
                decades=math.log10(max(G) / min(G)),
                all_attractive=all(s > 0 for s in sig))


def _chain():
    from interference import couplings, gravity, masses, mixing
    with contextlib.redirect_stdout(io.StringIO()):
        R = r.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
        c = couplings.derive(R, m)
        gravity.derive(R, m, mx, c)


GRAVITY_NODES = ("a1_bridge", "Sigma_vent(UV)", "Sigma_vent(broken)", "G_ratio")


def echo_scan(build=True):
    if build and "Sigma_vent(UV)" not in r.WEB:
        _chain()
    nodes = {}
    for name in GRAVITY_NODES:
        led = r.WEB[name]
        nodes[name] = dict(web_depth=led.web_depth,
                           max_term_depth=max((t.depth for t in led.terms), default=0),
                           kinds=sorted({t.kind for t in led.terms}),
                           n_terms=len(led.terms))
    mpl, cm = r.WEB["M_Pl_MeV"], r.WEB["common_mode"]
    return dict(nodes=nodes, anchor_web_depth=mpl.web_depth, common_mode_web_depth=cm.web_depth,
                common_mode_parents=sorted({t.parent for t in cm.terms if t.parent}))


def run(report=print, build=True):
    out = {}
    lad = ladder()
    out["ladder_ppm"] = lad
    assert abs(lad[1] + 1139.9) < 0.5 and abs(lad[2] - 21.5) < 0.1 and abs(lad[4] - 22.2) < 0.1
    report(f"R1 ladder: M_Pl against CODATA at web depth 1, 2, 4: "
           f"{lad[1]:+.1f}, {lad[2]:+.1f}, {lad[4]:+.1f} ppm.")

    at = attribution()
    out["attribution"] = at
    assert at["from_m_e"] < 1e-9 and at["from_alpha"] < 1e-11
    assert at["next_echo"] < 2e-9 and at["residual_over_next_echo"] > 1e4
    assert 1.5 < at["residual_over_sigma"] < 2.5
    report(f"R2 attribution: m_e moves M_Pl by {at['from_m_e']:.1e}, the Cs/Rb alpha split by "
           f"{at['from_alpha']:.1e}. The depth 4 step was {at['step_4']:.1e} and a depth 5 term is "
           f"of order {at['next_echo']:.1e}.  The residual {1e6 * at['residual']:+.1f} ppm is "
           f"{at['residual_over_next_echo']:.1e} next echoes and "
           f"{at['residual_over_sigma']:.2f} sigma of G's measurement. As an action shift "
           f"delta S = {at['delta_S']:.2e}.")

    rd = readings()
    out["readings"] = rd
    assert abs(rd["anchor_response"] - 1.0) < 1e-6
    assert abs(rd["ledger_anchored_response"]) < 1e-15
    assert abs(rd["ledger_placed_response"] + 0.5) < 1e-12
    report(f"R3 readings: the anchor moves {rd['anchor_response']:.3f} per unit of m_e (a reading); "
           f"the ledger's G at the anchored t0 moves {rd['ledger_anchored_response']:.0e} when Sigma "
           f"doubles (an identity). At a placed t0 it moves {rd['ledger_placed_response']:+.2f} "
           f"(a reading of the coarse estimate).")

    sh = silhouette()
    out["silhouette"] = sh
    assert sh["all_attractive"] and 1.5 < sh["decades"] < 2.5
    report(f"R4 coarse estimate: Sigma {sh['sigma_min']:.2f} to {sh['sigma_max']:.2f} over contacts and "
           f"the metric's own loop, attractive in every cell; G {sh['G_min']:.3f} to "
           f"{sh['G_max']:.2f} M_b^-2 over the Planckian band, {sh['decades']:.2f} decades.")

    es = echo_scan(build)
    out["echoes"] = es
    for name, n in es["nodes"].items():
        assert n["web_depth"] <= 1 and n["max_term_depth"] <= 1, name
    assert es["anchor_web_depth"] == 4 and es["common_mode_web_depth"] == 4
    assert es["common_mode_parents"] == ["eq_loop"]
    report("R5 echoes: " + "; ".join(
        f"{k} depth {v['web_depth']}, terms {v['n_terms']} ({'/'.join(v['kinds']) or 'live base'}) "
        f"to depth {v['max_term_depth']}" for k, v in es["nodes"].items())
        + f".  The anchor node and the common mode sit at web depth {es['anchor_web_depth']} and "
          f"{es['common_mode_web_depth']}, the common mode reading {es['common_mode_parents'][0]}.")
    report("RESULT: one closed recursion (the anchor) gives gravity's number at 22 ppm, at the "
           "scale of G's own measurement and far above what the ladder or its measured inputs "
           "can supply. The ledger gives a coarse estimate spanning about two decades with a fixed sign. The "
           "stiffness nodes read at depth one and the common mode at depth four.")
    return out


if __name__ == "__main__":
    run()
