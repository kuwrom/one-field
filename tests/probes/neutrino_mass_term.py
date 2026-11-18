"""The neutrino mass term: which one the two F4 singlets close on.

THE QUESTION. The 26 of F4 splits as 8_v + 8_s + 8_c + 1 + 1
(masses.py). Beside the three lepton doublets, the two singlets are the
only neutral rows the engine holds. Three mass terms can be written on
that content, so there are eight ways to switch them on. Which does the
engine's label theorem admit, and which do the data admit?

PREMISES. P1 the content: three lepton doublets carrying lepton label
1, the label of the charged leptons, and the two singlets, with a
label l that is not assumed. P2 the terms. a is a mass of the singlets with themselves
and changes the label by 2 l. b joins a singlet to a doublet and
changes it by l - 1. c is a mass of the doublets with themselves and
changes it by 2. P3 the label theorem of root.EchoTerm (WHY THE PROTON
IS STABLE): "no edge kind changes a knot's species.  So baryon and
lepton number are conserved node labels of the graph." A term that
changes the label is an identity changing edge, and the engine holds
none. P4 two measured splittings (reference.NUFIT_OSCILLATION, NuFIT
6.0), so at least two neutrinos have mass. P5 comparison side only:
neutrinoless double beta decay m_bb < 28 to 122 meV at 90% CL
(KamLAND-Zen, PRL 135, 262501 (2025)), the kinematic mass m_beta < 0.45 eV
at 90% CL (KATRIN, Science 388, 180 (2025)), the sum of masses
< 0.12 eV at 95% CL and N_eff = 2.99 +- 0.17 (Planck 2018, A&A 641,
A6 (2020)), and the share of N_eff carried by right handed halves that
couple through their mass alone, 7.5e-12 (m/0.1 eV)^2 (Luo, Rodejohann
and Xu, JCAP 03 (2021) 082).

STATEMENTS (the numbers below check each one).

  N1 (the menu). 2^3 = 8 structures. For each, a generic mass matrix on
     the five neutral states gives the number of massive neutrinos
     (the rank of its doublet columns) and their kind: massive states
     in equal pairs are full fermions (Dirac), single ones are their
     own conjugates (Majorana). b alone gives two Dirac neutrinos and
     one massless. a with b, the type-I seesaw, gives two Majorana
     neutrinos and one massless. Every structure with c gives three
     Majorana neutrinos. The empty structure and a alone give none.
  N2 (the label gate, blind). A structure passes for a singlet label l
     when each of its terms leaves the label unchanged. b passes for
     l = 1 alone, a for l = 0 alone, c for no l. Three structures pass
     for some l: the empty one, a alone, b alone. a and b never pass
     together, so the seesaw fails for every label the singlets could
     carry. root.EchoTerm rejects any edge kind beyond ratio, echo and
     vent, and every edge of the solved web is one of the three.
  N3 (the data gate). The two splittings stand at 37 and 93 standard
     deviations, so a structure needs two massive neutrinos. Six pass.
     The empty structure and a alone fail.
  N4 (both gates). One structure of eight passes both: b alone, with
     l = 1. The singlets are the right handed halves of two Dirac
     neutrinos. Neither gate decides alone: the label gate leaves 3,
     the data gate leaves 6.
  N5 (what the survivor says). The mass matrix has rank 2, so the
     lightest neutrino is massless (normal ordering m1 = 0, inverted
     m3 = 0). No term changes the label by 2, so neutrinoless double
     beta decay is absent: m_bb = 0. For normal ordering the measured
     splittings and the engine's angles give a sum of 58.78 meV
     against < 120, and m_beta = 8.81 meV against < 450.
  N6 (what the data separate today). The seesaw row has the same
     masses and angles and differs in m_bb alone: 1.51 to 3.70 meV for
     normal ordering, 7.6 times below the strongest limit at its
     strictest. The data pass both rows. The label gate leaves one.
  N7 (the right handed halves stay out of the early bath). Their only
     coupling is the mass, y = sqrt 2 m3/v_EW = 2.9e-13. With every
     coupling and phase space factor set to 1, the rate y^2 T over the
     expansion rate at T = v_EW is 2.4e-10. They never share the bath.
     The published share is 1.9e-12 at m3, eleven orders inside the
     width of the measured N_eff.

RESULT. Of the eight mass structures on the engine's neutral content,
one conserves the lepton label and gives the two measured splittings:
the singlet doublet term with the singlets carrying the lepton label.
The neutrinos are Dirac, the lightest is massless, and m_bb = 0. The
type-I seesaw fails the label gate for every label,
since its two terms need l = 0 and l = 1 at once. The data do not
separate the two rows today. An observed neutrinoless double beta decay kills the reading,
and so do three established nonzero masses.

Usage: python3 tests/probes/neutrino_mass_term.py
"""

import contextlib
import io
import itertools
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import masses, mixing  # noqa: E402
from interference import root as r  # noqa: E402

TERMS = ("a", "b", "c")
# P2: the change of the lepton label under each term, for singlet label l.
CHANGE = {"a": lambda l: 2 * l, "b": lambda l: l - 1, "c": lambda l: 2}
LABELS = range(-3, 4)
N_DOUBLETS, N_SINGLETS = r.d11, r.n26 - 3 * (r.d11**2 - 1)

# P5, comparison side, in meV unless stated.
M_BB_LIMIT = (28.0, 122.0)      # KamLAND-Zen, 90% CL, over the nuclear matrix elements
M_BETA_LIMIT = 450.0            # KATRIN, 90% CL
SUM_LIMIT = 120.0               # Planck 2018, 95% CL
N_EFF = (2.99, 0.17)            # Planck 2018, 68% CL
N_EFF_DIRAC = 7.5e-12           # Luo, Rodejohann and Xu, at a mass of 100 meV, scaling as the mass squared
G_STAR = 106.75                 # Standard Model degrees of freedom above the electroweak scale


def structures():
    """N1: the eight structures, as the tuple of terms switched on."""
    return [tuple(t for t, x in zip(TERMS, bits) if x) for bits in itertools.product((0, 1), repeat=3)]


def mass_matrix(on, rng):
    """A generic complex symmetric mass matrix on (three doublets, two singlets)."""
    def block(rows, cols):
        return rng.normal(size=(rows, cols)) + 1j * rng.normal(size=(rows, cols))
    n, k = N_DOUBLETS, N_SINGLETS
    M = np.zeros((n + k, n + k), complex)
    if "c" in on:
        x = block(n, n)
        M[:n, :n] = x + x.T
    if "a" in on:
        x = block(k, k)
        M[n:, n:] = x + x.T
    if "b" in on:
        x = block(n, k)
        M[:n, n:], M[n:, :n] = x, x.T
    return M


def spectrum(on, seed=7):
    """N1: massive neutrinos (rank of the doublet columns) and their kind."""
    M = mass_matrix(on, np.random.default_rng(seed))
    massive = int(np.linalg.matrix_rank(M[:, :N_DOUBLETS], tol=1e-9))
    values = [s for s in np.linalg.svd(M, compute_uv=False) if s > 1e-9]
    paired = all(any(i != j and abs(s - t) < 1e-9 * s for j, t in enumerate(values))
                 for i, s in enumerate(values))
    kind = "none" if not massive else "Dirac" if paired else "Majorana"
    return {"massive": massive, "massless": N_DOUBLETS - massive, "kind": kind}


def conserving_labels(on):
    """N2: the singlet labels for which every term of the structure conserves the label."""
    return [l for l in LABELS if all(CHANGE[t](l) == 0 for t in on)]


def engine_kinds():
    """N2: the edge taxonomy read from the engine."""
    try:
        r.EchoTerm(["nu", "nu bar"], 1.0, 1, kind="majorana")
        rejected = False
    except AssertionError:
        rejected = True
    kinds = {t.kind for led in r.WEB.values() if isinstance(led, r.Ledger) for t in led.terms}
    return rejected, kinds


def benchmark():
    """N5 to N7: the normal ordering numbers from the engine's angles and the measured splittings."""
    with contextlib.redirect_stdout(io.StringIO()):
        R = r.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
    s12, s13 = mx["sin2_12"], mx["sin2_13"]
    m2, m3 = (math.sqrt(r.NUFIT_OSCILLATION[k][0]) * 1e3 for k in ("dm2_21", "dm2_31"))
    t2, t3 = s12 * (1 - s13) * m2, s13 * m3
    y = math.sqrt(2) * m3 * 1e-12 / m["v_EW_pred_GeV"]
    bath = y**2 * r.M_Pl_GeV / (1.66 * math.sqrt(G_STAR) * m["v_EW_pred_GeV"])
    return {"sum": m2 + m3, "m_beta": math.sqrt((1 - s13) * s12 * m2**2 + s13 * m3**2),
            "m_bb_seesaw": (abs(t2 - t3), t2 + t3), "yukawa": y, "bath": bath,
            "n_eff": N_EFF_DIRAC * (m3 / 100.0) ** 2}


def run(report=print):
    out = {}

    # N1: the menu.
    menu = structures()
    rows = {on: spectrum(on) for on in menu}
    assert len(menu) == 8 and (N_DOUBLETS, N_SINGLETS) == (3, 2)
    assert rows[("b",)] == {"massive": 2, "massless": 1, "kind": "Dirac"}
    assert rows[("a", "b")] == {"massive": 2, "massless": 1, "kind": "Majorana"}
    assert all(rows[on] == {"massive": 3, "massless": 0, "kind": "Majorana"} for on in menu if "c" in on)
    assert rows[()]["massive"] == rows[("a",)]["massive"] == 0
    assert all(spectrum(on, seed) == rows[on] for on in menu for seed in (1, 2, 3))
    for on in menu:
        x = rows[on]
        report(f"N1 {'+'.join(on) or 'none':6s} massive {x['massive']}, massless {x['massless']}, {x['kind']}")

    # N2: the label gate.
    labels = {on: conserving_labels(on) for on in menu}
    label_pass = [on for on in menu if labels[on]]
    rejected, kinds = engine_kinds()
    assert rejected and kinds <= {"ratio", "echo", "vent"}
    assert label_pass == [(), ("b",), ("a",)]
    assert labels[("b",)] == [1] and labels[("a",)] == [0] and labels[("a", "b")] == []
    assert all(not labels[on] for on in menu if "c" in on)
    report("N2 label gate passes: " + ", ".join(f"{'+'.join(on) or 'none'} (l = {labels[on][0] if on else 'any'})"
                                                for on in label_pass) + ". The seesaw a+b passes for no l.")

    # N3: the data gate.
    sig = [r.NUFIT_OSCILLATION[k][0] / r.NUFIT_OSCILLATION[k][1] for k in ("dm2_21", "dm2_31")]
    assert min(sig) > 30
    data_pass = [on for on in menu if rows[on]["massive"] >= 2]
    assert len(data_pass) == 6 and () not in data_pass and ("a",) not in data_pass
    report(f"N3 splittings at {sig[0]:.0f} and {sig[1]:.0f} sigma: {len(data_pass)} structures have two massive neutrinos")

    # N4: both gates.
    both = [on for on in menu if on in label_pass and on in data_pass]
    assert both == [("b",)]
    counts = {"menu": len(menu), "label": len(label_pass), "data": len(data_pass), "both": len(both)}
    survivor = {"terms": "+".join(both[0]), "label": labels[both[0]][0], **rows[both[0]]}
    assert counts == {"menu": 8, "label": 3, "data": 6, "both": 1}
    out.update(counts=counts, survivor=survivor, seesaw_labels=labels[("a", "b")])
    report(f"N4 {counts['menu']} structures, {counts['label']} pass the label gate, {counts['data']} the data gate, "
           f"{counts['both']} both: {survivor['terms']} with l = {survivor['label']}, {survivor['kind']}")

    # N5, N6: the survivor's numbers and the seesaw row's.
    b = benchmark()
    lo, hi = b["m_bb_seesaw"]
    margin = M_BB_LIMIT[0] / hi
    assert abs(b["sum"] - 58.78) < 0.05 and abs(b["m_beta"] - 8.81) < 0.05
    assert abs(lo - 1.51) < 0.05 and abs(hi - 3.70) < 0.05
    assert b["sum"] < SUM_LIMIT and b["m_beta"] < M_BETA_LIMIT and 7.0 < margin < 8.0
    out.update(m_bb=0.0, m_bb_seesaw=(lo, hi), m_bb_margin=margin, sum=b["sum"], m_beta=b["m_beta"])
    report(f"N5 survivor: m_bb = 0, sum {b['sum']:.2f} meV (< {SUM_LIMIT:.0f}), m_beta {b['m_beta']:.2f} meV (< {M_BETA_LIMIT:.0f})")
    report(f"N6 seesaw row: m_bb {lo:.2f} to {hi:.2f} meV, {margin:.1f} times below the limit {M_BB_LIMIT[0]:.0f} meV")

    # N7: the right handed halves and the early bath.
    assert 2.8e-13 < b["yukawa"] < 3.0e-13 and 2.0e-10 < b["bath"] < 3.0e-10
    assert 1.8e-12 < b["n_eff"] < 2.0e-12 and b["n_eff"] < 1e-10 * N_EFF[1]
    out.update(yukawa=b["yukawa"], bath=b["bath"], n_eff=b["n_eff"])
    report(f"N7 y = {b['yukawa']:.2e}, rate over expansion at v_EW {b['bath']:.1e}. "
           f"Published share of N_eff {b['n_eff']:.1e} against a width of {N_EFF[1]}")
    report("RESULT: one mass structure of eight conserves the lepton label and gives the two measured splittings. "
           "The neutrinos are Dirac, the lightest is massless, and m_bb = 0.")
    return out


if __name__ == "__main__":
    run()
