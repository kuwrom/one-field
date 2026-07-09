"""
surplus_edges.py: web completeness, run in the SURPLUS direction
(derivation program D4 of the gaps audit).

The scorecard runs the grammar TOWARD the PDG: nine targets, nine
hits.  Completeness demands the opposite scan too: enumerate every
admissible edge the grammar permits and ask what states it generates
that are NOT in the PDG.  Each such state is one of exactly two
things, and the probe labels which:

  PREDICTION: an admissible edge whose state is absent from the PDG,
      a registered exposure (a new state the framework says exists,
      or a mass relation not yet measured).
  FALSIFIER: an edge the grammar EXCLUDES by lemma whose state is
      also absent: it must stay absent, or the lemma that killed it
      is wrong (the 64 GeV and 139 GeV terminus ghosts are the
      canonical examples: already dead against the top).

The menus come from tests/edge_grammar.py (the committed D1 spec).
The selection lemmas from interference/words.py (N-ality
superselection: steps must be triality-neutral. Conversion lemma:
bases are pure-tower integer ranks, mixed words are vent-side).

Usage: python3 surplus_edges.py    (fast: pure linear algebra)
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "interference"))

import edge_grammar as eg  # noqa: E402
from root import PDG_MASSES  # noqa: E402

TOL = 0.01  # the scorecard's own band


def _pdg_match(mass_mev):
    return next((n for n, v in PDG_MASSES.items()
                 if abs(mass_mev / v - 1.0) < TOL), None)


def run(report=print):
    report("SURPLUS EDGES: the grammar enumerated AWAY from the data")
    report("=" * 70)

    m_lep = {1: PDG_MASSES["e"], 2: PDG_MASSES["mu"], 3: PDG_MASSES["tau"]}
    used_bases = {(1, "ff"), (1, "aa"), (2, "ffa"), (3, "f4+a4")}

    rows = []      # (gen, word, count, admissible?, reason, state MeV)
    # generations 1 and 2: words of length n+1 over {f, a}
    for gen, length in ((1, 2), (2, 3)):
        for w, cnt in sorted(eg.word_menu(length).items()):
            # conversion lemma: a BASE is a pure tower at gen 1 /
            # seed+neutral-steps at gen 2. Any in-word channel switch
            # beyond the seed+step grammar imports an amplitude
            if gen == 1:
                adm = w in ("ff", "aa")
                why = "pure seed" if adm else "mixed seed: conv. lemma"
            else:
                adm = w in ("ffa", "aaa")   # seed + one NEUTRAL step
                why = ("seed+adjoint step" if adm
                       else "charged or in-tower step: N-ality/conv.")
            rows.append((gen, w, cnt, adm, why, cnt * m_lep[gen]))
    # generation 3: the committed Z2-symmetric terminus menu
    for name, cnt in sorted(eg.terminus_menu().items()):
        adm = name == "f4+a4"
        why = ("pure-tower Z2 orbit sum" if adm
               else "mixed word: conversion lemma (vent-side)")
        rows.append((3, name, cnt, adm, why, cnt * m_lep[3]))

    surplus, falsifiers = [], []
    report(f"  {'gen':>3s} {'word':<8s} {'count':>5s} {'status':<11s} "
           f"{'state':>12s}  disposition")
    for gen, w, cnt, adm, why, mev in rows:
        match = _pdg_match(mev)
        used = (gen, w) in used_bases
        if used:
            disp = f"USED -> {match} (canon)"
        elif adm and match is None:
            disp = "SURPLUS -> register as prediction or derive exclusion"
            surplus.append((gen, w, mev))
        elif adm and match is not None:
            disp = f"admissible, lands on {match} (report)"
            surplus.append((gen, w, mev))
        elif match is not None:
            # excluded WORD whose count is degenerate with a used one:
            # the lemma excludes the word, not the value. Harmless
            disp = f"excluded word ({why}). Value degenerate with {match}"
        else:
            disp = f"excluded ({why}). Value must stay absent"
            falsifiers.append((gen, w, mev))
        report(f"  {gen:>3d} {w:<8s} {cnt:>5d} "
               f"{'ADMISSIBLE' if adm else 'excluded':<11s} "
               f"{mev/1e3:>9.2f} GeV  {disp}")

    report("-" * 70)
    report(f"  used bases:  {sorted(used_bases)}, all admissible")
    report(f"  surplus (the free prediction channel): "
           f"{[(g, w, round(m/1e3, 2)) for g, w, m in surplus]}")
    report(f"  falsifiers (excluded, must stay absent): "
           f"{[(g, w, round(m/1e3, 2)) for g, w, m in falsifiers]}")

    # ── frozen assertions ────────────────────────────────────────────
    # 1. Every USED base is admissible under the lemmas
    adm_set = {(g, w) for g, w, c, a, y, m in rows if a}
    assert used_bases <= adm_set, "a canon base violates its own grammar"
    # 2. The terminus ghosts are dead and must stay dead: the excluded
    #    mixed words land at ~64 and ~139 GeV, matching no PDG state
    ghosts = sorted(round(m / 1e3) for g, w, m in falsifiers if g == 3)
    assert ghosts == [64, 139], f"terminus ghosts moved: {ghosts}"
    assert all(_pdg_match(m) is None for g, w, m in falsifiers if g == 3)
    # 3. The surplus channel is FINITE and enumerated: the grammar
    #    cannot quietly grow new states (menu sizes frozen)
    assert len(rows) == 4 + 8 + 3
    return dict(surplus=surplus, falsifiers=falsifiers, rows=rows)


if __name__ == "__main__":
    run()
