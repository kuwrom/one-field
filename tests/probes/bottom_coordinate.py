"""The bottom mass: the relation the engine uses and the entries it meets.

THE QUESTION. The engine computes m_b from the (c, b, t) Koide relation
with Q = Q0 + h11/K^3 = 289/432, using its own m_c and m_t. The result
is 4193.8 MeV, +1.5 standard deviations from the PDG m_b(m_b) =
4183(7). The triplet meets three kinds of table entry: m_c(m_c),
m_b(m_b), and the direct top mass. How much of the residual belongs to
the relation, how much to the inputs, and why these entries?

PREMISES. P1 the relation of interference/masses.py, the negative free
branch bracketed between m_c and m_t. P2 reference values: PDG 2024
m_c(m_c) = 1273.0(4.6), m_t = 172.57(29) GeV direct, m_b(m_b) =
4183(7); FLAG 2024 m_b(m_b) = 4200(14) for Nf = 2+1+1 and 4171(20) for
Nf = 2+1 (arXiv:2411.04268, Eqs. 61, 63). P3 a reading is the object
read against itself, with every echo of its channels included, and it
does not run (masses.py). A scale from outside the object would be a
second thing, which self-presence excludes (interference/__init__.py).
P4 the engine's own alpha_s(M_Z), G_F = 1/(sqrt 2 v_EW^2), M_W, M_Z
and m_t(m_t) (couplings, higgs_ew). P5 for B6 and B7 only, the loop
series of the strong channel on the comparison side: beta and gamma_m
to four loops (van Ritbergen, Vermaseren and Larin, PLB 400 (1997)
379; Vermaseren, Larin and van Ritbergen, PLB 405 (1997) 327), the
self entry series to four loops (Marquard et al., PRD 94 (2016)
074025). No reference mass enters B4 to B7.

STATEMENTS (the numbers below check each one).

  B1 (the lever). d ln m_b / d ln m_t = +1.61 and d ln m_b / d ln m_c =
     -0.61. One GeV of top moves m_b by 39 MeV, so the 8 to 11 MeV
     residual equals a quarter of a GeV of top.
  B2 (the relation fed with measured inputs). With the measured m_c and
     m_t the relation gives 4197.8 +- 14.6 MeV. With the spread of the
     inputs included, that is +0.9 standard deviations from 4183(7).
     This number tests the relation itself, independent of the
     engine's m_c and m_t.
  B3 (the references). The engine's 4193.8 is +1.5 sigma from PDG,
     -0.4 sigma from FLAG Nf = 2+1+1 and +1.1 sigma from FLAG Nf = 2+1.
  B4 (the entries an object fixes alone). By P3 a reading meets an
     entry the object fixes by itself. A table has two. The complete
     entry (the table's pole) holds all of the object's self
     interference and exists when the object is never bound into a
     hadron: colourless, or decaying before it binds (Gamma > Lambda).
     This is the measured object, not the word ladder, whose confining
     interface the top's word closes behind (words.py). The entry cut at the object's own mass
     (the table's m(m)) exists for every colourless object, and for a
     quark whose mass lies above confinement (m > Lambda). An entry at
     an outside scale (a shared mu, 2 GeV) is the third kind. Lambda
     is the Landau point of the engine's own alpha_s with thresholds at
     its own m_b and m_c: 0.20 GeV at one loop, 0.40 GeV at two. The
     top width from the engine's G_F and M_W is 1.48 GeV, 3.7 times the
     larger Lambda. The weak widths of the other quarks are below 2e-10
     GeV. m_c is 3.1 times the larger Lambda, m_s is 0.46 times the
     smaller. The classes hold for either Lambda.
  B5 (exhaustion). The reading is completed, so it meets the complete
     entry where one exists, the entry at its own mass where only that
     exists, and an outside scale only where the object fixes no entry
     of its own. Of the 3^9 = 19683 assignments of the nine fermions,
     324 use entries that exist and one follows this order: e, mu, tau
     and t complete, c and b at their own mass, u, d, s at an outside
     scale with ratios alone as content. It is the assignment
     masses.py uses. Each of the six possible orders leaves one
     assignment, and this order alone leaves that one. One shared
     scale for all reads six objects
     against an outside scale. The complete entry for all fails on c,
     b, u, d, s, which have none. Their own mass for all fails on u,
     d, s and cuts the leptons and the top. Complete leptons with
     every quark at its own mass cuts the top. Every rival fails on
     existence or order, never on a mass.
  B6 (the complete entry of a confined quark). The strong channel's
     self entry series at the object's own mass, with the engine's
     alpha_s, reaches about the confinement scale within its four
     known terms for every quark: the smallest is 0.20 GeV on t, 0.14
     GeV on b, 0.21 GeV on c. That is 0.12% of the top, inside which
     the engine's top
     sits (-0.036%), and 3.4% of the bottom and 16% of the charm,
     larger than any quark residual. So c and b fix no complete entry
     at the resolution of their readings.
  B7 (the same values at one outside scale). Carried to mu = m_t by
     beta order transport, which is not web depth, the engine's values
     give m_c(m_t) = 0.597 to 0.623 GeV, m_b(m_t) = 2.730 to 2.753 GeV
     and Q(c,b,t) = 0.7185 to 0.7205 at two to four loops, and m_b/m_c
     = 4.42 to 4.57 against 3.29 at their own masses. Tables at one
     shared scale quote Q near 0.72. That is the same spectrum read
     against an outside scale, which B5 excludes for c, b and t.

RESULT. The residual on m_b is inside the spread of the reference
itself. The relation fed with measured inputs is at +0.9 sigma, and the
two lattice averages bracket the engine's value. The entries the
triplet meets are forced: one order and an existence check computed
from the engine's own values leave one assignment, the one masses.py
uses. Q(c,b,t) = 289/432 holds at the entries the objects fix alone.
Read against one outside scale, the same values give 0.72.

Usage: python3 tests/probes/bottom_coordinate.py
"""

import contextlib
import io
import itertools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import couplings, higgs_ew, masses  # noqa: E402
from interference import root as r  # noqa: E402

M_C, M_T = (1273.0, 4.6), (172570.0, 290.0)
REFS = {"PDG m_b(m_b)": (4183.0, 7.0), "FLAG Nf=2+1+1": (4200.0, 14.0), "FLAG Nf=2+1": (4171.0, 20.0)}

FERMIONS = ("e", "mu", "tau", "u", "d", "s", "c", "b", "t")
COLOURED = {"u", "d", "s", "c", "b", "t"}
# B5 order: the complete entry, the entry at the object's own mass, an outside scale.
ENTRIES = ("complete", "own mass", "outside")
MASSES_PY = {"e": "complete", "mu": "complete", "tau": "complete", "t": "complete",
             "c": "own mass", "b": "own mass", "u": "outside", "d": "outside", "s": "outside"}
Z3, Z4, Z5 = 1.2020569031595942, math.pi**4 / 90, 1.0369277551433699


def engine():
    """The engine's masses in GeV, and its alpha_s(M_Z), G_F, M_W, M_Z, m_t(m_t)."""
    with contextlib.redirect_stdout(io.StringIO()):
        R = r.derive()
        m = masses.derive(R)
        c = couplings.derive(R, m)
        h = higgs_ew.derive_higgs_ew(m, c)
    key = {"e": "m_e", "mu": "m_mu", "tau": "m_tau"}
    return {"m": {f: m[key.get(f, "m_" + f)] / 1e3 for f in FERMIONS},
            "a_MZ": c["alpha_s_MZ_thresh"], "M_Z": h["M_Z_derived"], "M_W": h["M_W_derived"],
            "G_F": 1 / (math.sqrt(2) * m["v_EW_pred_GeV"] ** 2), "m_t_run": h["m_t_rge_coordinate_GeV"]}


def bottom(m_c, m_t, Q=float(r.Q0 + r.h11 / r.K**3)):
    sc, st = math.sqrt(m_c), math.sqrt(m_t)
    a, b = Q - 1, 2 * Q * (sc + st)
    c = (Q - 1) * (m_c + m_t) + 2 * Q * sc * st
    for sign in (1, -1):
        x = (-b + sign * math.sqrt(b * b - 4 * a * c)) / (2 * a)
        if m_c < x * x < m_t:
            return x * x
    raise ValueError("no bracketed root")


def confinement(E):
    """B4: the Landau point of the engine's own alpha_s below m_c, at one and two loops."""
    def run(a, mu0, mu1, nf):
        return couplings._rk4_run(lambda x: -couplings._b0(nf) * x * x - couplings._b1(nf) * x**3,
                                  a, math.log(mu0), math.log(mu1))
    m = E["m"]
    a_c = run(run(E["a_MZ"], E["M_Z"], m["b"], 5), m["b"], m["c"], 4)
    b0, b1 = couplings._b0(3), couplings._b1(3)
    one = m["c"] * math.exp(-1 / (b0 * a_c))
    two = m["c"] * math.exp(-1 / (b0 * a_c) - b1 / b0**2 * math.log(b1 * a_c / (b0 + b1 * a_c)))
    return one, two


def widths(E):
    """B4: the top width at leading order, and a bound on the other quarks' weak widths."""
    m, x = E["m"], (E["M_W"] / E["m"]["t"]) ** 2
    top = E["G_F"] * m["t"] ** 3 / (8 * math.pi * math.sqrt(2)) * (1 - x) ** 2 * (1 + 2 * x)
    return top, {q: 6 * E["G_F"] ** 2 * m[q] ** 5 / (192 * math.pi**3) for q in COLOURED - {"t"}}


def existing(E, lam):
    """B4: the entries each fermion fixes, for confinement scale lam."""
    top, bound = widths(E)
    gamma = {"t": top, **bound}
    out = {}
    for f in FERMIONS:
        have = {"outside"}
        if f not in COLOURED or gamma[f] > lam:
            have.add("complete")
        if f not in COLOURED or E["m"][f] > lam:
            have.add("own mass")
        out[f] = have
    return out


def ordered(have):
    """B5: the first entry, in ENTRIES order, that exists."""
    return next(k for k in ENTRIES if k in have)


def exhaust(have):
    """B5: every assignment, by existence, then by the order."""
    menu = list(itertools.product(ENTRIES, repeat=len(FERMIONS)))
    alive = [dict(zip(FERMIONS, a)) for a in menu if all(a[i] in have[f] for i, f in enumerate(FERMIONS))]
    best = {f: ordered(have[f]) for f in FERMIONS}
    return len(menu), len(alive), [a for a in alive if a == best]


def rivals(have):
    """B5: the named rivals, and where each fails existence or order."""
    named = {"one shared scale": {f: "outside" for f in FERMIONS},
             "complete for all": {f: "complete" for f in FERMIONS},
             "own mass for all": {f: "own mass" for f in FERMIONS},
             "complete leptons, quarks at own mass":
                 {f: ("own mass" if f in COLOURED else "complete") for f in FERMIONS}}
    return {name: ([f for f in FERMIONS if a[f] not in have[f]],
                   [f for f in FERMIONS if a[f] in have[f] and a[f] != ordered(have[f])])
            for name, a in named.items()}


def _series(nf):
    """beta_i, gamma_i in a/pi with mu^2 d/d mu^2, to four loops (P5)."""
    beta = [(11 - 2 * nf / 3) / 4, (102 - 38 * nf / 3) / 16,
            (2857 / 2 - 5033 / 18 * nf + 325 / 54 * nf**2) / 64,
            ((149753 / 6 + 3564 * Z3) - (1078361 / 162 + 6508 / 27 * Z3) * nf
             + (50065 / 162 + 6472 / 81 * Z3) * nf**2 + 1093 / 729 * nf**3) / 256]
    gamma = [1.0, (202 / 3 - 20 * nf / 9) / 16,
             (1249 + (-2216 / 27 - 160 / 3 * Z3) * nf - 140 / 81 * nf**2) / 64,
             (4603055 / 162 + 135680 / 27 * Z3 - 8800 * Z5
              + (-91723 / 27 - 34192 / 9 * Z3 + 880 * Z4 + 18400 / 9 * Z5) * nf
              + (5242 / 243 + 800 / 9 * Z3 - 160 / 3 * Z4) * nf**2
              + (-332 / 243 + 64 / 27 * Z3) * nf**3) / 256]
    return beta, gamma


def transport(a0, m0, mu0, mu1, nf, loops, n=4000):
    """alpha_s and a mass carried from mu0 to mu1 by beta order transport (RK4 in ln mu^2)."""
    beta, gamma = (s[:loops] for s in _series(nf))

    def f(x):
        return (-sum(b * x ** (i + 2) for i, b in enumerate(beta)),
                -sum(g * x ** (i + 1) for i, g in enumerate(gamma)))
    x, lm = a0 / math.pi, math.log(m0)
    h = (2 * math.log(mu1) - 2 * math.log(mu0)) / n
    for _ in range(n):
        k1 = f(x)
        k2 = f(x + h / 2 * k1[0])
        k3 = f(x + h / 2 * k2[0])
        k4 = f(x + h * k3[0])
        x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        lm += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return x * math.pi, math.exp(lm)


def self_entry_floor(E):
    """B6: the smallest term of the strong self entry series at each quark's own mass, in GeV."""
    m = E["m"]
    a_t, _ = transport(E["a_MZ"], 1.0, E["M_Z"], E["m_t_run"], 5, 4)
    a_b, _ = transport(E["a_MZ"], 1.0, E["M_Z"], m["b"], 5, 4)
    a_c, _ = transport(a_b, 1.0, m["b"], m["c"], 4, 4)
    out = {}
    for q, a, nl, mq in (("t", a_t, 5, E["m_t_run"]), ("b", a_b, 4, m["b"]), ("c", a_c, 3, m["c"])):
        k = (4 / 3, 13.4434 - 1.0414 * nl, 190.595 - 26.655 * nl + 0.6527 * nl**2,
             3567.60 - 745.721 * nl + 43.3963 * nl**2 - 0.678141 * nl**3)
        out[q] = min(c * (a / math.pi) ** (i + 1) for i, c in enumerate(k)) * mq
    return out


def koide(x, y, z):
    return (x + y + z) / (math.sqrt(x) + math.sqrt(y) + math.sqrt(z)) ** 2


def outside_scale(E, loops):
    """B7: the engine's values carried to mu = m_t."""
    m = E["m"]
    a_b, _ = transport(E["a_MZ"], 1.0, E["M_Z"], m["b"], 5, loops)
    a_c, _ = transport(a_b, 1.0, m["b"], m["c"], 4, loops)
    _, c_at_b = transport(a_c, m["c"], m["c"], m["b"], 4, loops)
    _, c_at_t = transport(a_b, c_at_b, m["b"], E["m_t_run"], 5, loops)
    _, b_at_t = transport(a_b, m["b"], m["b"], E["m_t_run"], 5, loops)
    return {"m_c": c_at_t, "m_b": b_at_t, "Q": koide(c_at_t, b_at_t, E["m_t_run"]), "b/c": m["b"] / c_at_b}


def run(report=print):
    out = {}
    E = engine()
    ENGINE = {k: E["m"][k[2:]] * 1e3 for k in ("m_c", "m_t", "m_b")}
    assert abs(ENGINE["m_b"] - 4193.8) < 0.05
    assert abs(bottom(ENGINE["m_c"], ENGINE["m_t"]) - ENGINE["m_b"]) < 1e-6

    # B1: the lever.
    eps = 1e-6
    b0 = bottom(ENGINE["m_c"], ENGINE["m_t"])
    dt = (bottom(ENGINE["m_c"], ENGINE["m_t"] * (1 + eps)) / b0 - 1) / eps
    dc = (bottom(ENGINE["m_c"] * (1 + eps), ENGINE["m_t"]) / b0 - 1) / eps
    per_gev = dt * b0 / ENGINE["m_t"] * 1e3
    assert abs(dt - 1.608) < 0.005 and abs(dc + 0.608) < 0.005 and abs(per_gev - 39.1) < 0.3
    out.update(dlnt=dt, dlnc=dc, per_gev=per_gev)
    report(f"B1 d ln m_b/d ln m_t = {dt:+.3f}, d ln m_b/d ln m_c = {dc:+.3f}: {per_gev:.1f} MeV per GeV of top")

    # B2: the relation fed with measured inputs.
    fed = bottom(M_C[0], M_T[0])
    spread = math.hypot(bottom(M_C[0] + M_C[1], M_T[0]) - fed, bottom(M_C[0], M_T[0] + M_T[1]) - fed)
    pull = (fed - REFS["PDG m_b(m_b)"][0]) / math.hypot(REFS["PDG m_b(m_b)"][1], spread)
    assert abs(fed - 4197.8) < 0.1 and abs(spread - 14.6) < 0.1 and abs(pull - 0.91) < 0.01
    out.update(fed=fed, spread=spread, fed_pull=pull)
    report(f"B2 the relation on measured m_c, m_t: {fed:.1f} +- {spread:.1f} MeV, {pull:+.2f} sigma of PDG")

    # B3: the references.
    pulls = {k: (ENGINE["m_b"] - v) / s for k, (v, s) in REFS.items()}
    assert abs(pulls["PDG m_b(m_b)"] - 1.54) < 0.01 and abs(pulls["FLAG Nf=2+1+1"] + 0.44) < 0.01
    out.update(pulls=pulls)
    report("B3 the engine's 4193.8 MeV: " + ", ".join(f"{k} {p:+.2f}" for k, p in pulls.items()))

    # B4: the entries each object fixes, from the engine's own values.
    lam = confinement(E)
    top, bound = widths(E)
    assert 0.15 < lam[0] < lam[1] < 0.5 and 1.3 < top < 1.6 and max(bound.values()) < 2e-10
    classes = [existing(E, x) for x in lam]
    assert classes[0] == classes[1]
    have = classes[0]
    margins = {"Gamma_t/Lambda": top / lam[1], "m_c/Lambda": E["m"]["c"] / lam[1],
               "m_s/Lambda": E["m"]["s"] / lam[0]}
    assert margins["Gamma_t/Lambda"] > 3 and margins["m_c/Lambda"] > 3 and margins["m_s/Lambda"] < 0.5
    out.update(Lambda=lam, Gamma_t=top, margins=margins)
    report(f"B4 Lambda {lam[0]:.3f} to {lam[1]:.3f} GeV, Gamma_t {top:.3f} GeV, "
           + ", ".join(f"{k} {v:.2f}" for k, v in margins.items()))

    # B5: the exhaustion.
    total, alive, survivors = exhaust(have)
    assert (total, alive, len(survivors)) == (19683, 324, 1) and survivors[0] == MASSES_PY
    by_order = {o: {f: next(k for k in o if k in have[f]) for f in FERMIONS} for o in itertools.permutations(ENTRIES)}
    assert [o for o, a in by_order.items() if a == MASSES_PY] == [ENTRIES]
    fails = rivals(have)
    assert fails["one shared scale"] == ([], ["e", "mu", "tau", "c", "b", "t"])
    assert fails["complete for all"][0] == ["u", "d", "s", "c", "b"]
    assert fails["own mass for all"] == (["u", "d", "s"], ["e", "mu", "tau", "t"])
    assert fails["complete leptons, quarks at own mass"] == (["u", "d", "s"], ["t"])
    out.update(assignments=(total, alive, len(survivors)), survivor=survivors[0], rivals=fails)
    report(f"B5 {total} assignments, {alive} exist, {len(survivors)} in order: the entries masses.py uses")
    for name, (e, o) in fails.items():
        report(f"   {name}: no such entry on {e or 'none'}, out of order on {o or 'none'}")

    # B6: the complete entry of a confined quark.
    floor = self_entry_floor(E)
    rel = {q: floor[q] / (E["m_t_run"] if q == "t" else E["m"][q]) for q in floor}
    assert all(0.1 < floor[q] < 0.25 for q in floor)
    assert rel["t"] < 0.002 and rel["b"] > 0.03 and rel["c"] > 0.1
    out.update(floor=floor, floor_rel=rel)
    report("B6 smallest self entry term: " + ", ".join(f"{q} {floor[q]:.3f} GeV ({100 * rel[q]:.2f}%)" for q in floor))

    # B7: the same values at one outside scale.
    rows = {L: outside_scale(E, L) for L in (2, 3, 4)}
    own = E["m"]["b"] / E["m"]["c"]
    for x in rows.values():
        assert 0.59 < x["m_c"] < 0.63 and 2.72 < x["m_b"] < 2.76 and 0.717 < x["Q"] < 0.722 and 4.40 < x["b/c"] < 4.60
    assert abs(koide(E["m"]["c"], E["m"]["b"], E["m"]["t"]) - 289 / 432) < 1e-9 and abs(own - 3.292) < 0.002
    out.update(outside=rows, own_ratio=own)
    for L, x in rows.items():
        report(f"B7 {L} loops at m_t: m_c {x['m_c']:.3f}, m_b {x['m_b']:.3f} GeV, Q {x['Q']:.4f}, "
               f"m_b/m_c {x['b/c']:.3f} (own masses {own:.3f}, Q 289/432)")
    return out


if __name__ == "__main__":
    run()
