"""
edge_grammar.py: the machine-readable admissible-edge specification.

WHY THIS FILE EXISTS (derivation program D1 of the gaps audit): the
framework's central negative claim, "zero dimensionless free
parameters", is only checkable if "admissible edge" is a formal,
finite, machine-readable definition rather than an impression left by
the code.  This module IS that definition.  It generalises the depth
pre-registration rule: the complete menus below are committed here,
in code, independent of any evaluation against data.  The null-model
test (test_coverage.py::test_grammar_null_model_scorecard_not_cheap)
Monte-Carlos over exactly these menus. The surplus-edge probe
(tests/probes/surplus_edges.py) enumerates them in the forward
direction.  Neither can widen the grammar without editing this file,
and editing this file is a registered event (the menus are frozen by
tests).

THE SPEC (data first, helpers below):

  GRAMMAR["integers"]      the four irreducible integers
  GRAMMAR["depth"]         Hurwitz gate: composite depth <= 3.
                           a new depth requires a new derived edge
  GRAMMAR["kinds"]         edge taxonomy of the one graph
  GRAMMAR["registration"]  the depth pre-registration rule, verbatim
  GRAMMAR["operations"]    the admissible edge-factor constructors

Every constructor draws from FINITE menus computed from the four
integers (weights of the ten SU(3)_3 primaries, boundary-walk counts
of channel words on the D6 nimrep, Koide structures Q0 + h/K^p,
small dictionary rationals).  Nothing here reads a measured value.
"""

import itertools
import math
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "interference")))

from root import d10, d11, n7, n26, K, h10, h11, delta, Q0  # noqa: E402
from words import N_FUND  # noqa: E402


# ═══════════════════════════════════════════════════════════════════════
#  The machine-readable spec
# ═══════════════════════════════════════════════════════════════════════

GRAMMAR = dict(
    integers=dict(d10=int(d10), d11=int(d11), n7=int(n7), n26=int(n26)),
    depth=dict(
        max_depth=3,
        rule="HURWITZ GATE: composite depth <= 3 (phase of phase of "
             "phase). A deeper path requires a NEW derived edge (a "
             "theorem), never a deeper ledger (root.EchoTerm raises)."),
    kinds=("ratio", "echo", "vent"),
    registration=(
        "DEPTH PRE-REGISTRATION RULE: before any new echo depth is "
        "evaluated in any channel, the complete grammar-permitted edge "
        "menu at that depth, with multiplicities, is committed first. "
        "A term identified after a residual is known, off-menu, is "
        "inadmissible.  This module is the committed menu."),
    operations=dict(
        circulant="sqrt(m_k) ~ 1 + (B/A) cos(theta + 2 pi k/3). "
                  "B/A = sqrt(d10) (Fano CG, fixed). theta from "
                  "theta_menu() (weights of the ten primaries and "
                  "their standard multiples), one theta per closure",
        word_multiplier="m_q(n) = [count(w) + corr] * m_lep(n). w a "
                        "channel word over {f, a} with the length rule "
                        "|w_n| = n + 1. count(w) = boundary-walk count "
                        "on the D6 nimrep (word_menu)",
        koide_structure="Q = Q0 + h/K^p, h in {h(fund), h(adj)} = "
                        "{2/9, 1/2}, p in 1..4 (koide_menu). Inverted "
                        "for the held-out mass at fixed partners",
        correction="one multiplicative vent (1 + c), c from "
                   "correction_menu() (small dictionary rationals), "
                   "or none",
    ),
)


# ═══════════════════════════════════════════════════════════════════════
#  Menu enumerators (each is complete, finite, and data-blind)
# ═══════════════════════════════════════════════════════════════════════

def primary_weights():
    """Conformal weights h = C2/(k+3) of all ten SU(3)_3 primaries."""
    prims = [(a, b) for a in range(4) for b in range(4) if a + b <= 3]
    return {(a, b): (a*a + b*b + a*b + 3*a + 3*b) / (3.0 * (3 + 3))
            for (a, b) in prims}


def theta_menu():
    """The natural angle menu of the MTC: every nonzero primary weight
    h and its standard multiples (h, pi*h, 2*pi*h).  Mirrors
    tests/probes/theta_menu.py. ~27 candidates before the positive-
    branch cut."""
    menu = set()
    for h in primary_weights().values():
        if h == 0.0:
            continue
        menu.update((h, math.pi * h, 2 * math.pi * h))
    return sorted(menu)


def word_menu(length):
    """ALL channel words over {f, a} of the given length, with their
    boundary-walk counts on the D6 nimrep.  Returns {word: count}.
    The grammar's selection rules (neutral steps, purity at the
    terminus) act downstream. The menu itself is the full cartesian
    set, so the null model cannot be accused of pre-narrowing."""
    N = N_FUND
    A = N @ N.T - np.eye(6)
    mat = {"f": N, "a": A}
    out = {}
    for letters in itertools.product("fa", repeat=length):
        M = np.eye(6)
        for c in letters:
            M = M @ mat[c]
        out["".join(letters)] = int(round(M.sum(axis=1)[0]))
    return out


def terminus_menu():
    """Length-4 Z2-symmetric terminus candidates (the committed menu
    of words.py part 3): pure orbit sum, and the two mixed classes."""
    w = word_menu(4)
    return {
        "f4+a4": w["ffff"] + w["aaaa"],           # 97 (selected)
        "f2a2": w["ffaa"],                        # mixed, excluded
        "f3a+fa3": w["fffa"] + w["faaa"],         # mixed, excluded
    }


def koide_menu():
    """The eight Koide structures Q0 + h/K^p (test_coverage._koide_menu
    in Q-space)."""
    out = {}
    for hname, h in (("fund", h10), ("adj", h11)):
        for p in range(1, 5):
            out[(hname, p)] = float(Q0 + Fraction(h) / K**p)
    return out


def correction_menu():
    """The committed corr ladder of the mass sector (masses.py:
    h10 -> delta -> 1/2K) plus the altitude family 1/K^p, p in 1..3.
    Exactly these. A wider menu would have to be committed HERE
    first (the pre-registration rule)."""
    return {
        "h10": float(h10), "delta": float(delta),
        "1/2K": 1.0 / (2 * K),
        "1/K": 1.0 / K, "1/K2": 1.0 / K**2, "1/K3": 1.0 / K**3,
    }


if __name__ == "__main__":
    print("EDGE GRAMMAR (machine-readable admissible-edge spec)")
    print("=" * 64)
    print(f"  integers: {GRAMMAR['integers']}")
    print(f"  depth:    <= {GRAMMAR['depth']['max_depth']} (Hurwitz)")
    print(f"  kinds:    {GRAMMAR['kinds']}")
    print(f"  theta menu: {len(theta_menu())} candidates")
    for L in (2, 3, 4):
        print(f"  word menu |w| = {L}: {word_menu(L)}")
    print(f"  terminus menu: {terminus_menu()}")
    print(f"  koide menu: { {k: round(v, 6) for k, v in koide_menu().items()} }")
    print(f"  correction menu: "
          f"{ {k: round(v, 6) for k, v in correction_menu().items()} }")
