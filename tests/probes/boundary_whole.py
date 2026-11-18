"""The boundary reading c: a named channel reads a share, the common mode reads the whole.

THE QUESTION. The anchor computes M_b, the scale at which
alpha_G2 = 1/(24 pi) holds, and G = 1/(c M_b^2) with c = 1 taken from
G itself (tests/audit/g_path_trials.py G10). In terms of records,
G M_b^2 is the coefficient of the exchange between two closures at the
boundary rate. It is the fraction by which each closure moves the
other's clock, at the distance over which its own phase turns one
radian (gates_as_records.py U). The engine already
defines such a coefficient for a named channel. What does the same
definition give for the common mode, which has no name?

PREMISES. P1 the engine's channel coupling: "closure probability times
the normalized dimension one current response" (interference/wiring.py),
read at the boundary, where the map is read whole (root.py). P2 the
forgetting: the common mode enters at second order, through the second
moment of the harmonic (interference/searchlight/protected_forgetting.py).
P3 the one list and its gates S and M (interference/boundary.py), which
use no measured G.

STATEMENTS (the numbers below check each one).

  B1 (three shares, three wholes). Each factor a named channel reads
     is a share of a sum, and each sum is one.
       closure probability: the pointed orbit holds 3/36 = 1/12 of the
         S-matrix's first row, whose squares sum to 1.
       Response: the part that closes locally is R(1) = 1/(2 pi), the
         remainder 1 - 1/(2 pi), together 1.
       Second moment: one relative channel holds E[v^2] = 1/2, and the
         compression to the common mode, P v^2 P, is 1.
     The engine's named couplings are products of shares:
     alpha_G2 = (1/12)(1/(2 pi)), xi = alpha_G2 E[v^2] h, and the
     baryon fraction R(1).
  B2 (the reading). A name is what picks a share. The common mode
     carries none, so it reads each sum whole, and its coefficient at
     the boundary is 1 x 1 x 1: G M_b^2 = 1, c = 1.
  B3 (the twelve readings, no measured G). Of the readings
     c = 1/(probability x response x moment) over the shares and
     wholes of B1, gates S and M pass at every contact for four:
     c = 1, 1.189, 2 and 2.379. Every reading with the pointed
     probability fails S, and so does the local response alone at the
     covariant contact (c = 2 pi against 4.16). S's threshold is soft
     (g_path_trials.py G10). Across its six conventions, c = 1 is the
     only one of the twelve readings that passes under all six. The
     other three pass under five, four and four.
  B4 (the three left). The remainder alone, c = 1.189, gives no
     weight to any named closure, and the local share alone gives no
     weight to the amount. The dark ratio is then infinite or zero.
     When both carry weight it is 2 pi - 1, the reading dark_sector.py
     already makes. One relative
     channel, c = 2, is not a compression to the common mode: P v^2 P
     sums both channels, for every harmonic.
  B5 (the reading moves). c = 2/B^2, with B the harmonic's amplitude.
     The whole second moment is one because B^2 = d10, which is the
     statement Q = 2/3. A harmonic with B^2 off by a part in a
     thousand moves c by the same fraction, and moves m_mu/m_e by 2.4
     percent. The measured ratio is known to 2.2e-8, so the leptons
     fix c to 9e-10, where G alone left it unbounded
     (tests/audit/anchor_account.py A3).
  B6 (the same whole in the anchor). The lepton closure's common
     mode, the mean of its three rates, is M_b e^(-S) times
     C_EFF (1 + P v^2 P) = 1 at depth one. The measured mean is 0.12
     percent below it, the size of the depth two factor alpha/2 pi.
  B7 (the test). On the depth four web the common mode reads this
     whole as the count of the lepton circulant,
     G M_b^2 = 3 Q_l - 1 = 1 - 7.50 (alpha/2 pi)^2 (root.py, THE COMMON
     MODE; common_mode_count.py). The anchor's G then comes out 54.6 ppm
     below the recommended value, 2.43 of its standard deviations. The
     nearest of the other three readings, c = 1.189, comes out 16
     percent below.
  B8 (the channel map is the moment's whole). On the searchlight's
     circulant the common mode's count is P v^2 P = (B/A)^2/2 = 3 Q - 1
     for every B/A. So the channel map B/A = sqrt d10 (octonions.py),
     the Koide value Q = 2/3 and the whole second moment are one
     statement, an identification counted once on G's path. Three
     readers fix it:
     the amplitude gate with no data (B/A in [1.30541, 1.47367],
     amplitude_gate.py, G M_b^2 between 0.852 and 1.086), G at its
     standard deviation (B/A to 1.1e-5), and m_mu/m_e (B/A to 5e-10).
     The channel map's weight W = 1/d(fund), the quantum trace, is the
     whole reading in another form. The classical trace gives
     B/A = 1/sqrt 2, a count of 1/4, and a G off by a factor of four
     from the whole. The amplitude gate excludes it with no data.

RESULT. c = 1 is the statement that the common mode reads no share.
The menu is the engine's own: its shares and their sums. c = 1 is the
only one of the twelve readings that passes gate S under every
convention, with no measured G. The second moment is the only factor
that is a reading: it moves with the harmonic, and the leptons fix it
to 9e-10. So G's 54.6 ppm cannot sit in the whole, including the
count's depth four term (the carrier read through the three lepton
terms, -10.1 ppm). The distance belongs to the rule that combines the
measurements (g_measurements.py G6). That moment is the channel map
B/A read by the common mode (B8). The other two sums are one by
unitarity and by the definition of the shares. The gates use no
measured G. The menu of shares and sums is assembled with G's value
known.

Scope: the engine's coupling definition carried to the common mode at
the boundary, with the second moment as its third factor because the
common mode enters at second order. No dynamics of the exchange at the
boundary is computed.

Usage: python3 tests/probes/boundary_whole.py
"""

import itertools
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402
from interference import wiring  # noqa: E402
from interference.searchlight import protected_forgetting as pf  # noqa: E402

BAND = (0.6, 15.0)                           # gate M: t0 in M_b^-2 (gravity.py)
G_SIGMA = r.G_REL_SIGMA                      # relative uncertainty of the recommended G
MU_E_SIGMA = 2.2e-8                          # relative uncertainty of m_mu/m_e


def shares():
    """The three sums of B1: each named share and its whole."""
    d = wiring.su3_qdims(3)
    D2 = sum(v * v for v in d.values())
    pointed = sum(d[p] ** 2 for p in ((0, 0), (3, 0), (0, 3))) / D2
    local = wiring.normalized_current_response(1.0)
    _, total, mode = pf.moments(pf.harmonic())
    return dict(probability=dict(pointed=pointed, whole=sum(v * v for v in d.values()) / D2),
                response=dict(local=local, remainder=1.0 - local, whole=local + (1.0 - local)),
                moment=dict(mode=mode, whole=total))


def gates(c):
    """Gates S and M on the one list for G = 1/(c M_b^2), at every contact. No measured G."""
    one, xi = B.rows(r.d10, r.d11, r.n7, r.n26), float(r.xi_bridge)
    S = all(B.sigma(one, 1.0, ct, xi) >= 2 * math.pi * c for ct in B.CONTACTS)
    t0 = []
    for ct in B.CONTACTS:
        f = lambda t: B.sigma(one, t, ct, xi) - 2 * math.pi * c * t
        lo, hi = 1e-9, 200.0
        for _ in range(200):                              # sigma falls with t0: one root
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        t0.append(0.5 * (lo + hi))
    return S, all(BAND[0] < t < BAND[1] for t in t0)


def conventions_passed(c):
    """Of gate S's six conventions (g_path_trials.py G10), the number under which c passes at every contact."""
    one, xi = B.rows(r.d10, r.d11, r.n7, r.n26), float(r.xi_bridge)
    x_half = brentq(lambda x: expn(2, x) - 0.5, 1e-3, 5.0)            # a row at M_b at half weight
    conventions = [(th, F) for th in (x_half, 1.0, 1.0 / x_half) for F in (1.0, 2 * math.log(2))]
    return sum(all(B.sigma(one, th, ct, xi) * F >= 2 * math.pi * c * th for ct in B.CONTACTS)
               for th, F in conventions)


def lepton_ratio(B2):
    """m_mu/m_e of the harmonic with squared amplitude B2 at the engine's angle, depth one."""
    v = pf.harmonic(B=math.sqrt(B2), theta=float(r.h10))
    m = sorted((1.0 + v) ** 2)
    return m[1] / m[0]


def run(report=print):
    out = {}

    # B1: shares and wholes.
    s = shares()
    assert abs(s["probability"]["pointed"] - 1 / 12) < 1e-12 and abs(s["probability"]["whole"] - 1) < 1e-12
    assert abs(s["response"]["local"] - 1 / (2 * math.pi)) < 1e-15 and s["response"]["whole"] == 1.0
    assert abs(s["moment"]["mode"] - 0.5) < 1e-12 and abs(s["moment"]["whole"] - 1) < 1e-12
    alpha = s["probability"]["pointed"] * s["response"]["local"]
    xi = alpha * s["moment"]["mode"] * float(r.h_bridge)
    dark = s["response"]["remainder"] / s["response"]["local"]
    assert abs(alpha - r.alpha_G2_Pl) < 1e-15 and abs(xi - float(r.xi_bridge)) < 1e-15
    assert abs(dark - (2 * math.pi - 1)) < 1e-12
    out.update(shares=s, alpha_G2=alpha, xi=xi, dark_ratio=dark)
    report(f"B1 shares and wholes: probability {s['probability']['pointed']:.4f} of {s['probability']['whole']:.0f}; "
           f"response {s['response']['local']:.4f} and {s['response']['remainder']:.4f} of "
           f"{s['response']['whole']:.0f}. Moment {s['moment']['mode']:.1f} of {s['moment']['whole']:.0f}; "
           f"alpha_G2, xi and the dark ratio {dark:.4f} are the engine's")

    # B2: the common mode's coefficient.
    whole = s["probability"]["whole"] * s["response"]["whole"] * s["moment"]["whole"]
    assert abs(whole - 1) < 1e-12
    out.update(coefficient=whole, c=1 / whole)
    report(f"B2 the common mode reads each sum whole: G M_b^2 = {whole:.12f}, c = {1 / whole:.12f}")

    # B3: the twelve readings against S and M.
    menu = {}
    for (pn, p), (rn, rs), (mn, m) in itertools.product(*(s[k].items() for k in ("probability", "response", "moment"))):
        c = 1 / (p * rs * m)
        menu[(pn, rn, mn)] = (c, *gates(c))
    passing = sorted(c for c, S, M in menu.values() if S and M)
    assert len(menu) == 12 and len(passing) == 4
    assert np.allclose(passing, [1.0, 1 / (1 - 1 / (2 * math.pi)), 2.0, 2 / (1 - 1 / (2 * math.pi))])
    assert not any(S for (pn, _, _), (_, S, _) in menu.items() if pn == "pointed")
    assert not menu[("whole", "local", "whole")][1]
    one, xil = B.rows(r.d10, r.d11, r.n7, r.n26), float(r.xi_bridge)
    edge = B.sigma(one, 1.0, 1.0, xil) / (2 * math.pi)
    soft = {k: conventions_passed(v[0]) for k, v in menu.items()}
    assert [k for k, n in soft.items() if n == 6] == [("whole", "whole", "whole")]
    assert sorted(soft[k] for k, v in menu.items() if v[1] and v[2]) == [4, 4, 5, 6]
    out.update(menu={"/".join(k): v for k, v in menu.items()}, passing=passing, s_edge=edge,
               conventions={"/".join(k): n for k, n in soft.items()})
    report("B3 twelve readings, gates S and M with no measured G: pass for c = "
           + ", ".join(f"{c:.3f}" for c in passing) + f". The pointed probability fails S in all six, and the "
           f"local response alone gives c = {2 * math.pi:.2f} against S's {edge:.2f} at the covariant contact. "
           "Across S's six conventions c = 1 alone passes all six, the other three "
           + ", ".join(str(soft[k]) for k, v in sorted(menu.items(), key=lambda kv: kv[1][0])
                       if v[1] and v[2] and v[0] > 1.0))

    # B4: the three left.
    weights = {"both": (1.0, 1.0), "remainder alone": (0.0, 1.0), "local alone": (1.0, 0.0)}
    ratios = {k: (math.inf if wl == 0 else wr * s["response"]["remainder"] / (wl * s["response"]["local"]))
              for k, (wl, wr) in weights.items()}
    rng = np.random.default_rng(7)
    compress = max(abs(pf.moments(pf.harmonic(B=b, theta=t))[1] - b * b / 2)
                   for b, t in zip(rng.uniform(0.3, 2.0, 20), rng.uniform(0, 2 * math.pi, 20)))
    assert abs(ratios["both"] - dark) < 1e-12 and ratios["remainder alone"] == math.inf
    assert ratios["local alone"] == 0.0 and compress < 1e-12
    out.update(dark_readings=ratios, compression_error=compress)
    report(f"B4 dark ratio {ratios['both']:.4f} when both shares weigh, infinite for the remainder alone, zero for "
           f"the local share alone; P v^2 P = B^2/2, both channels, for twenty harmonics to {compress:.0e}")

    # B5: the reading moves with the harmonic, and the leptons fix it.
    eps = 1e-3
    c_of = lambda B2: 2.0 / B2
    moved_c = c_of(2 * (1 + eps)) / c_of(2.0) - 1
    slope = (lepton_ratio(2 * (1 + eps)) / lepton_ratio(2 * (1 - eps)) - 1) / (2 * eps)
    held = MU_E_SIGMA / abs(slope)
    Q = lambda B2: 1 / 3 + B2 / 6
    assert abs(moved_c + eps) < 1e-5 and abs(slope) > 1 and held < 1e-8 and Q(2.0) == 2 / 3
    assert abs(lepton_ratio(2.0) / (r.PDG_MASSES["mu"] / r.PDG_MASSES["e"]) - 1) < 2e-3
    out.update(c_moves=moved_c, ratio_slope=slope, c_held=held, Q=Q(2.0), Q_moved=Q(2 * (1 + eps)))
    report(f"B5 c = 2/B^2: B^2 up by {eps:.0e} moves c by {moved_c:+.1e}, Q to {Q(2 * (1 + eps)):.6f} and "
           f"m_mu/m_e by {slope * eps:+.1e}. The measured ratio fixes c to {held:.0e}")

    # B6: the lepton closure's common mode.
    unit = r.C_EFF * (1 + s["moment"]["whole"])
    mean = (r.PDG_MASSES["e"] + r.PDG_MASSES["mu"] + r.PDG_MASSES["tau"]) / 3
    frame = r.M_Pl_MeV * math.exp(-r.S_lepton)
    assert abs(unit - 1.0) < 1e-12 and abs(mean / frame - 1) < 2e-3
    out.update(common_mode_coefficient=unit, mean_over_frame=mean / frame - 1)
    report(f"B6 the lepton closure's common mode: C_EFF (1 + P v^2 P) = {unit:.0f}. The measured mean rate is "
           f"{100 * (mean / frame - 1):+.2f}% from M_b e^-S")

    # B7: the test.
    ppm = 1e6 * (r.G_PRED / r.G_CODATA - 1)
    nearest = passing[1]
    assert abs(ppm + 54.6) < 0.1 and abs(1 / nearest - 1) > 0.15
    out.update(g_ppm=ppm, g_sigma=ppm * 1e-6 / G_SIGMA, nearest=nearest, nearest_shift=1 / nearest - 1)
    report(f"B7 c = 1 read as the depth four count: G at {ppm:+.1f} ppm, {ppm * 1e-6 / G_SIGMA:+.1f} standard deviations. c = {nearest:.3f} "
           f"would put it at {100 * (1 / nearest - 1):+.0f}%")

    # B8: the channel map is the moment's whole.
    gate = (1.30541, 1.47367)                                            # amplitude_gate.py
    counted = {ba: pf.moments(pf.harmonic(B=ba, theta=float(r.h10)))[1] for ba in (*gate, r.BA_ratio)}
    assert all(abs(v - ba * ba / 2) < 1e-12 for ba, v in counted.items())
    assert abs(counted[r.BA_ratio] - 1) < 1e-12 and abs(r.WEB["common_mode"].value_at(dict(r.WEB.state, eq_loop=0.0)) - 1) < 1e-12
    held = {"amplitude gate": (gate[1] - gate[0]) / 2 / r.BA_ratio, "G": G_SIGMA / 2, "m_mu/m_e": held / 2}
    assert abs(counted[gate[0]] - 0.852) < 1e-3 and abs(counted[gate[1]] - 1.086) < 1e-3
    assert abs(held["G"] - 1.1e-5) < 1e-6 and held["m_mu/m_e"] < 1e-9
    classical = 1 / math.sqrt(2)                                       # octonions.py, the classical trace
    classical_count = pf.moments(pf.harmonic(B=classical, theta=float(r.h10)))[1]
    assert abs(classical_count - 0.25) < 1e-12 and not gate[0] <= classical <= gate[1]
    out.update(channel_map=counted, ba_held=held, classical_count=classical_count)
    report(f"B8 the count on the searchlight's circulant is (B/A)^2/2: {counted[gate[0]]:.3f} and {counted[gate[1]]:.3f} "
           f"at the amplitude gate's edges, 1 at sqrt 2, {classical_count:.2f} at the classical trace. B/A fixed by "
           + ", ".join(f"{k} to {v:.1e}" for k, v in held.items()))

    report("RESULT: c = 1 is the common mode reading no share. Its three sums are read by the running, the dark "
           "ratio and the lepton ratios, and the shares it rejects fail without G. The reading is made with G's "
           "value known.")
    return out


if __name__ == "__main__":
    run()
