"""The anchor's G against the measurements of G, one by one, and the precision block's joint statistic.

THE QUESTION.  The anchor inversion gives G = 6.673936e-11 when the
common mode is taken as the count of the lepton circulant on the depth
four web (root.py, THE COMMON MODE; common_mode_count.py).  This is 2.4
standard deviations below the recommended value 6.674 30(15)e-11
(root.py G_PRED).  That uncertainty is not the measurements' own: the
sixteen inputs disagree, and the adjustment multiplies every one of
their uncertainties by 3.9.  Where does the anchor's value sit in the
set itself, and what is the joint statistic of the precision block?

PREMISES.  P1 the sixteen input data of the adjustment, as printed in
Table XXX (30) of Mohr, Newell, Taylor and Tiesinga, "CODATA recommended
values of the fundamental physical constants: 2022" (arXiv:2409.03787),
with its three correlation coefficients.  P2 its Sec. XV.A: "the same
3.9 expansion factor is used in this adjustment", and "The factor 3.9
reduces all normalized residuals to 2 or less". Before the expansion
"the largest three normalized residuals were 7.7, 6.8, and 4.8 for
BIPM-14, JILA-18, and BIPM-01, respectively".  P3 the precision block of
root.py: G, m_mu, m_tau and v_EW against G_F.  The fine structure
constant is not in the block. Its references disagree with each other
and are stated on their own.  P4 two published combinations of the
same kind of data that add an unknown spread to each input in place of
multiplying its uncertainty: the consensus value 6.67408(24) of
Merkatas, Toman, Possolo and Schlamminger (arXiv:1905.09551), and
6.6740(15) of Rinaldi, Middleton, Del Pozzo and Gair
(arXiv:2209.07416).  Bodnar and Eriksson (arXiv:2104.01977) compare
the two models on these data and support the multiplicative one.  A
third rule, a median weighted by the stated uncertainties, gives
"G = 6.674486(48)" (Muravyov, Khudonogova and Shabramov 2025, as
reported by the IUPAP working group on G in 2025).

STATEMENTS (the numbers below check each one).

  G1 (the set).  The sixteen values give the generalized least squares
     mean 6.674300(38), chi-square 195.7 on 15 degrees of freedom: the
     set is not one population.  Multiplied by 3.9 the mean's
     uncertainty is 0.000148, the recommended 6.674 30(15).
  G2 (the anchor against the mean).  The anchor's value lies 9.6 of
     the mean's own standard deviations below it and 2.47 of the
     expanded ones, 0.000148 before rounding. Against the recommended
     value as printed, 6.67430(15), it lies 2.43 below.  The 2.43 is the
     reading with the expansion.
     Without it the distance is larger, not smaller.
  G3 (the anchor in the set).  Seven values lie below the anchor's and
     nine above: it sits between the seventh and eighth, MSL-03 and
     LANL-97, inside the scatter.  The individual pulls run from -5.3
     (JILA-18) to +10.0 (BIPM-14).  The anchor's value lies below the
     five most precise inputs (relative uncertainty below 2e-5), by 2.6
     to 7.1 of their own standard deviations.  The set neither singles the
     anchor's value out nor excludes it. Its precise half prefers a
     larger G.
  G4 (the block).  The four pulls are G -2.43, m_mu +0.12, m_tau -0.23
     and v_EW +0.08: chi-square 6.0 on 4, probability 0.20, and one
     pull of four beyond 2.43 has probability 0.06.  The block's
     chi-square is G's: the other three carry 0.07 of it.
  G5 (the fine structure constant, apart).  1/alpha sits +0.13
     standard deviations from the caesium recoil value, -14 from the
     rubidium one, -6.1 from the recommended blend and -7.8 from the
     electron moment value.  It enters no joint statistic here: the
     references differ by more than the pulls measure.
  G6 (the 2.43 comes from the combination rule).  The same sixteen values
     under other rules, each with the anchor's distance in that rule's
     standard deviations: an added spread fitted to chi-square = 15
     (Mandel and Paule), 6.67389(26), +0.2. The moment estimate of the
     spread (DerSimonian and Laird), 6.67399(16), -0.4. The unweighted
     mean, 6.67377(27), +0.6. The median, 6.67408(23), -0.6.  The
     published values of P4 put it at -0.6, 0.0 and -11.5.  The
     multiplicative rule keeps the stated weights and gives -2.43. Every
     rule that lets the scatter set the weights gives less than one,
     and the median weighted by the stated uncertainties gives eleven.
     The distance is a statement about the weights, and the literature
     has not settled them.
  G7 (what fails).  The same comparison rejects the anchor read through
     web depth one (G above every one of the sixteen, 103 expanded
     standard deviations) and the boundary read as 1/(4 pi G), the one
     rival gate S leaves (tests/audit/g_path_trials.py G10).

RESULT.  The 2.43 standard deviations are the anchor's distance from a
mean whose inputs scatter by 3.9 times their stated uncertainties.  They
belong to the rule that keeps the stated weights.  Under an added
spread, the median, or the published consensus values the distance is
below one.  In the set the anchor's value lies between the seventh and
eighth of the sixteen, and below the precise half.  A measurement at
1e-5 that the others confirm decides the question (root.py, THE
PRECISION BLOCK AND WHAT ENDS IT).  Nothing in the ladder can move the
value (gravity_resolution.py R2).

Usage: python3 tests/probes/g_measurements.py     (about 20 s, runs the chain)
"""

import contextlib
import io
import math
import os
import sys

import numpy as np
from scipy.stats import chi2, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import masses  # noqa: E402
from interference import root as r  # noqa: E402
from interference.reference import PDG_EW  # noqa: E402
import gravity_resolution as gr  # noqa: E402
import vertex_pinning as vp  # noqa: E402

# Table XXX (30) of arXiv:2409.03787, in 1e-11 m^3 kg^-1 s^-2: label, value, standard uncertainty.
TABLE = (
    ("NIST-82", 6.67248, 0.00043), ("TR&D-96", 6.6729, 0.0005), ("LANL-97", 6.67398, 0.00070),
    ("UWash-00", 6.674255, 0.000092), ("BIPM-01", 6.67559, 0.00027), ("UWup-02", 6.67422, 0.00098),
    ("MSL-03", 6.67387, 0.00027), ("HUST-05", 6.67222, 0.00087), ("UZur-06", 6.67425, 0.00012),
    ("HUST-09", 6.67349, 0.00018), ("BIPM-14", 6.67554, 0.00016), ("LENS-14", 6.67191, 0.00099),
    ("UCI-14", 6.67435, 0.00013), ("HUSTT-18", 6.674184, 0.000078), ("HUSTA-18", 6.674484, 0.000077),
    ("JILA-18", 6.67260, 0.00025),
)
CORRELATIONS = (("NIST-82", "LANL-97", 0.351), ("HUST-05", "HUST-09", 0.134), ("HUST-09", "HUSTT-18", 0.068))
EXPANSION = 3.9                                   # Sec. XV.A
RECOMMENDED = (6.67430, 0.00015)
PRECISE = 2.0e-5                                  # relative uncertainty below which the precise five lie
M_TAU_PDG = (1776.93, 0.09)                       # MeV, PDG 2024 (reference.py holds the value)
PUBLISHED = (("Merkatas et al.", 6.67408, 0.00024), ("Rinaldi et al.", 6.6740, 0.0015),
             ("Muravyov et al.", 6.674486, 0.000048))   # P4
INV_ALPHA = (("caesium recoil", 137.035999046, 27e-9), ("rubidium recoil", 137.035999206, 11e-9),
             ("recommended blend", 137.035999177, 21e-9), ("electron moment", 137.035999166, 15e-9))


def the_set():
    """Generalized least squares mean of the sixteen, its uncertainty, and the chi-square about a value."""
    labels = [row[0] for row in TABLE]
    x = np.array([row[1] for row in TABLE])
    s = np.array([row[2] for row in TABLE])
    C = np.diag(s**2)
    for a, b, rho in CORRELATIONS:
        i, j = labels.index(a), labels.index(b)
        C[i, j] = C[j, i] = rho * s[i] * s[j]
    W = np.linalg.inv(C)
    one = np.ones(len(x))
    var = 1.0 / (one @ W @ one)
    mean = var * (one @ W @ x)
    about = lambda g: float((x - g) @ W @ (x - g))
    return mean, math.sqrt(var), about


def combinations():
    """The sixteen values under other combination rules: name -> (value, standard uncertainty)."""
    x = np.array([row[1] for row in TABLE])
    s = np.array([row[2] for row in TABLE])
    n = len(x)

    def spread(tau):
        w = 1.0 / (s**2 + tau**2)
        m = float((w * x).sum() / w.sum())
        return m, float(w.sum()**-0.5), float((w * (x - m)**2).sum())

    lo, hi = 1e-7, 1e-2
    for _ in range(200):                                   # the spread at which chi-square equals n - 1
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if spread(mid)[2] > n - 1 else (lo, mid)
    w = 1.0 / s**2
    moment = math.sqrt(max(0.0, (spread(0.0)[2] - (n - 1)) / (w.sum() - (w**2).sum() / w.sum())))
    med = float(np.median(x))
    robust = 1.4826 * float(np.median(np.abs(x - med)))
    return {"added spread, chi-square = n - 1": spread(lo)[:2], "added spread, moment estimate": spread(moment)[:2],
            "unweighted mean": (float(x.mean()), float(x.std(ddof=1) / math.sqrt(n))),
            "median": (med, 1.2533 * robust / math.sqrt(n))}


def comparison_rules():
    """The three readings every G comparison quotes, in SI: the recommended value, the added spread, the median."""
    rules = combinations()
    return {"recommended": (r.G_CODATA, r.G_CODATA_SIGMA),
            "spread added": tuple(1e-11 * v for v in rules["added spread, chi-square = n - 1"][:2]),
            "median": tuple(1e-11 * v for v in rules["median"][:2])}


def compare(G):
    """A value of G against the set: pulls, distance from the mean, and its place in the scatter."""
    mean, u, about = the_set()
    pulls = {name: (v - G) / s for name, v, s in TABLE}
    below = sum(v < G for _, v, _ in TABLE)
    expanded = (G - mean) / (EXPANSION * u)
    return dict(pulls=pulls, raw=(G - mean) / u, expanded=expanded, below=below, above=len(TABLE) - below,
                chi2=about(G), passes=abs(expanded) < 3.0 and 0 < below < len(TABLE))


def references():
    """The precision block's measured values and standard uncertainties: G in units of the CODATA value, masses in MeV, v_EW in GeV."""
    v_meas = (math.sqrt(2) * PDG_EW["G_F"])**-0.5
    return {"G": (r.G_CODATA, r.G_CODATA_SIGMA), "m_mu": vp.M_MU_CODATA, "m_tau": M_TAU_PDG,
            "v_EW": (v_meas, v_meas * 0.5 * PDG_EW["G_F_err"] / PDG_EW["G_F"])}


def block(build=True):
    """Pulls of root.py's precision block, from the chain's own values."""
    if build:
        with contextlib.redirect_stdout(io.StringIO()):
            m = masses.derive(r.derive())
    else:
        m = dict(m_mu=r.WEB.state["m_mu"], m_tau=r.WEB.state["m_tau"], v_EW_pred_GeV=r.WEB.state["v_EW_GeV"])
    v_meas = (math.sqrt(2) * PDG_EW["G_F"])**-0.5
    return {
        "G": (r.G_PRED / r.G_CODATA - 1.0) / r.G_REL_SIGMA,
        "m_mu": (m["m_mu"] - vp.M_MU_CODATA[0]) / vp.M_MU_CODATA[1],
        "m_tau": (m["m_tau"] - M_TAU_PDG[0]) / M_TAU_PDG[1],
        "v_EW": (m["v_EW_pred_GeV"] / v_meas - 1.0) / (0.5 * PDG_EW["G_F_err"] / PDG_EW["G_F"]),
    }


def run(report=print, build=True):
    out = {}
    G = r.G_PRED * 1e11
    assert abs(G - 6.673936) < 1e-6

    # G1: the set against itself.
    mean, u, about = the_set()
    chi_set = about(mean)
    assert abs(mean - 6.674300) < 2e-6 and abs(u - 3.8e-5) < 1e-6 and abs(chi_set - 195.7) < 0.5
    assert abs(mean - RECOMMENDED[0]) < 5e-6 and abs(EXPANSION * u - RECOMMENDED[1]) < 5e-6
    residuals = sorted(((v - mean) / s for _, v, s in TABLE), key=abs, reverse=True)
    assert [round(abs(x), 1) for x in residuals[:3]] == [7.7, 6.8, 4.8]          # P2, as the adjustment prints them
    out.update(mean=mean, mean_unc=u, chi2_set=chi_set, dof=len(TABLE) - 1)
    report(f"G1 set: mean {mean:.6f}({1e6 * u:.0f}), chi-square {chi_set:.1f} on {len(TABLE) - 1}. Times "
           f"{EXPANSION} the uncertainty is {EXPANSION * u:.6f}, the recommended {RECOMMENDED[0]:.5f}"
           f"({1e5 * RECOMMENDED[1]:.0f}). Largest residuals {residuals[0]:+.1f}, {residuals[1]:+.1f}, "
           f"{residuals[2]:+.1f}")

    # G2, G3: the anchor's value against the mean and in the set.
    c = compare(G)
    ordered = sorted(TABLE, key=lambda row: row[1])
    neighbours = (ordered[c["below"] - 1][0], ordered[c["below"]][0])
    precise = {name: c["pulls"][name] for name, v, s in TABLE if s / v <= PRECISE}
    assert abs(c["raw"] + 9.6) < 0.1 and abs(c["expanded"] + 2.47) < 0.01 and c["passes"]
    assert (c["below"], c["above"]) == (7, 9) and neighbours == ("MSL-03", "LANL-97")
    assert abs(min(c["pulls"].values()) + 5.34) < 0.02 and abs(max(c["pulls"].values()) - 10.03) < 0.02
    assert len(precise) == 5 and 2.6 < min(precise.values()) and max(precise.values()) < 7.2
    assert c["chi2"] > chi_set + 50
    out.update(G=G, raw=c["raw"], expanded=c["expanded"], below=c["below"], above=c["above"], pulls=c["pulls"],
               neighbours=neighbours, precise=precise, chi2_anchor=c["chi2"])
    report(f"G2 anchor against the mean: {G:.6f} is {c['raw']:+.1f} of the mean's standard deviations, "
           f"{c['expanded']:+.2f} of the expanded ones. Chi-square of the set about it {c['chi2']:.1f}")
    report(f"G3 anchor in the set: {c['below']} below, {c['above']} above, between {neighbours[0]} and "
           f"{neighbours[1]}. Pulls from {min(c['pulls'].values()):+.2f} to {max(c['pulls'].values()):+.2f}; "
           "the precise half " + ", ".join(f"{k} {v:+.2f}" for k, v in precise.items()))

    # G4: the block's joint statistic.
    b = block(build)
    chi_block = sum(v * v for v in b.values())
    p_block = float(chi2.sf(chi_block, len(b)))
    worst = max(abs(v) for v in b.values())
    p_worst = 1.0 - (1.0 - 2.0 * float(norm.sf(worst)))**len(b)
    assert abs(b["G"] + 2.43) < 0.01 and abs(b["m_mu"] - 0.12) < 0.01 and abs(b["m_tau"] + 0.23) < 0.01
    assert abs(b["v_EW"] - 0.08) < 0.01
    assert abs(chi_block - 5.97) < 0.02 and 0.19 < p_block < 0.21 and 0.05 < p_worst < 0.07
    assert b["G"]**2 / chi_block > 0.98
    out.update(block=b, chi2_block=chi_block, p_block=p_block, p_worst=p_worst)
    report("G4 block: " + ", ".join(f"{k} {v:+.2f}" for k, v in b.items())
           + f". Chi-square {chi_block:.2f} on {len(b)}, probability {p_block:.2f}. One pull of {len(b)} beyond "
           f"{worst:.2f}: {p_worst:.2f}; G carries {b['G']**2 / chi_block:.2f} of it")

    # G5: the fine structure constant, apart from the block.
    alpha = {name: (r.inv_alpha_phys - v) / s for name, v, s in INV_ALPHA}
    assert [round(alpha[k], 1) for k, _, _ in INV_ALPHA] == [0.1, -14.2, -6.1, -7.8]
    split = (INV_ALPHA[1][1] - INV_ALPHA[0][1]) / math.hypot(INV_ALPHA[0][2], INV_ALPHA[1][2])
    assert split > 5.0
    out.update(alpha=alpha, alpha_split=split)
    report("G5 1/alpha, apart: " + ", ".join(f"{k} {v:+.2f}" for k, v in alpha.items())
           + f". The two recoil values differ by {split:.1f} of their combined standard deviations")

    # G6: the same set under other combination rules, and the published consensus values.
    rules = {k: (v, u, (G - v) / u) for k, (v, u) in combinations().items()}
    published = {k: (G - v) / u for k, v, u in PUBLISHED}
    assert [round(p, 1) for _, _, p in rules.values()] == [0.2, -0.4, 0.6, -0.6]
    assert max(abs(p) for _, _, p in rules.values()) < 1.0
    assert abs(published["Merkatas et al."]) < 0.7 and abs(published["Rinaldi et al."]) < 0.1
    assert abs(published["Muravyov et al."] + 11.46) < 0.01
    assert abs(rules["median"][0] - 6.674082) < 1e-6 and abs(c["expanded"]) > 2.4
    out.update(rules=rules, published=published)
    report("G6 combination rule: " + "; ".join(f"{k} {v:.5f}({1e5 * u:.0f}), {p:+.1f}" for k, (v, u, p) in rules.items())
           + ". Published " + ", ".join(f"{k} {p:+.1f}" for k, p in published.items())
           + f". The multiplicative rule, the recommended value {(G - RECOMMENDED[0]) / RECOMMENDED[1]:+.2f}")

    # G7: readings the comparison rejects.
    depth_one = r.G_CODATA * 1e11 * (r.M_PL_CODATA_GEV / gr.m_pl(depth=1))**2
    fails = {"web depth one": compare(depth_one), "boundary as 1/(4 pi G)": compare(G / (4 * math.pi))}
    assert not any(f["passes"] for f in fails.values())
    assert fails["web depth one"]["below"] == len(TABLE) and fails["web depth one"]["expanded"] > 100
    out.update(fails={k: (f["expanded"], f["below"]) for k, f in fails.items()}, depth_one=depth_one)
    report("G7 rejected: " + "; ".join(
        f"{k} {f['expanded']:+.3g} expanded standard deviations, {f['below']} of {len(TABLE)} below"
        for k, f in fails.items()))

    report(f"RESULT: the anchor's {G:.6f} is {(G - RECOMMENDED[0]) / RECOMMENDED[1]:+.2f} of the recommended value, "
           f"{c['expanded']:+.2f} expanded and {c['raw']:+.1f} unexpanded standard "
           f"deviations from the mean of a set with chi-square {chi_set:.0f} on {len(TABLE) - 1}. It lies inside the "
           "scatter, below the five most precise inputs. Under an added spread, the median or the "
           "added spread consensus it lies within one standard deviation, under the uncertainty weighted median "
           "eleven below. The block's chi-square is G's.")
    return out


if __name__ == "__main__":
    run()
