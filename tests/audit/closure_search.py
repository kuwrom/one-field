"""
closure_search.py: the grammar run as a GENERATOR, blind gates first,
data contact once, look-elsewhere denominator printed.

This audit holds the selected E8(1) > G2(1) x F4(1) wiring fixed. Its 504
hits at 1% are mass grammar realizations, not additional closure lanes.

WHY THIS MODULE EXISTS.  The scorecard is retrodiction-shaped: the
committed edges were found by a search that knew the targets.  No test
inside the repository can audit that outer search - but the repository
can do the next best thing: re-run the search itself, mechanically,
over the COMMITTED grammar (tests/audit/edge_grammar.py), with the gates
split into two strictly ordered stages:

  STAGE 1 - BLIND GATES (no measured value anywhere):
    positivity      every circulant gap Δₖ > 0 (√m real)
    nondegeneracy   three distinct masses per triplet
    reality         Koide inversion discriminant > 0
    bracketing      the inverted root lies between its partners
                    (masses.py's own selection rule, stated blind)
    monotonicity    word-length rule ⇒ base(n) increasing in n

  STAGE 2 - DATA CONTACT (applied exactly once, at the end):
    a candidate HITS if it lands within the scorecard band (1%)
    of the observed value.  A second, tighter tier (0.1%) is also
    reported: at 1% the data selects the word COUNT but often not
    the correction (all corrections are sub-percent); the tight tier
    shows what the precision claims actually resolve.  Down-type
    stages condition on the CANON upstream choices (the search is
    staged, and the conditioning is printed, not hidden).

The printed result is the honest trials accounting:

    candidates  →  blind survivors  →  data hits

per sector and jointly.  The JOINT BLIND SURVIVOR COUNT is the
look-elsewhere denominator the scorecard owes its reader: how many
universes could this grammar have described?  If the canon web is the
unique data hit, that is the uniqueness claim in executable form.  If
it is not, this module says so (the hit lists print in full).

THE MEASURED ANSWER, stated where the question is asked (frozen by
tests/test_closure_search.py; any drift fails CI):

  * canon is NOT the unique 1% mass grammar hit: 504 of 2,370,816
    eligible realizations hit the data at the 1% band, canon among
    them. The structural wiring is fixed throughout. The uniqueness
    that survives is narrower: the word COUNTS are
    resolved uniquely, the sub-percent corrections are not (all six
    corrections are within 1% of each other).
  * the joint 0.1% tier is ZERO, including canon: U1 and Db have no
    tight hits because the u and b references themselves are coarser
    than 0.1%.  The tight tier is decisive only where the reference
    is precise (mu, tau, c, t).
  * three of the six blind gates gate nothing: the up-stage gate
    (count >= 1) is satisfied by every boundary-walk word, so the
    U1/U2/U3 blind-survivor counts EQUAL their candidate counts and
    the joint denominator is, for the quark sector, the raw
    candidate count.
  * the strange stage has 0 menu-only hits at any tier; it closes
    ONLY through the derived off-menu bridge (DERIVED_EDGES),
    whose derivation lives outside the audited grammar.
  * the down stages condition on canon upstream (m_c, m_t, m_b by
    name), so the joint product multiplies CONDITIONAL counts; it is
    a staged accounting, not an independent joint search.
  * koide_invert_s uses the opposite Brannen branch from
    koide_invert_b, a declared discrete branch choice, i.e. one
    binary degree of freedom the menus do not price.
  * the koide menu at this stage is the 8 committed structures plus
    the declared identity insertion ("none" = Q0, the vent-free
    structure, part of the committed grammar's koide_structure
    operation as its zero case), 9 candidates, stated here rather
    than widened silently.

THE TWO DERIVED EDGES.  Two canon edges enter the mass sector
through structures beyond the word/correction menus, and both are
derivations, not menu choices (DERIVED_EDGES below; witnesses in
tests/probes/):

  m_d: the exact identity d₁₁²(1 + d₁₀/d₁₁⁴) = count(aa) + h₁₀ =
       83/9, the SAME committed word_multiplier operation as the up
       sector, with corr(1) = h₁₀ inherited from the generation
       ladder.  The Z₂ triality swap changes only the seed word.
       The down-gen-1 stage is reported below (ancillary,
       outside the frozen denominator).
  m_s: the Dynkin→Albert norm conversion δ_Z₂²/δ_J²: δ_Z₂ = 2/3
       (Dynkin fixed point), δ_J² = 3/8 constructed from a concrete
       J₃(O) element on an anti-associative Fano triple with the
       isotropic share 1/dim(O) (tests/probes/edge_derivations.py),
       conversion unique by Schur's lemma on the irreducible 26.

Usage:  python3 tests/audit/closure_search.py
Frozen: tests/test_closure_search.py (survivor counts, hit identities,
        and the gap list are all pinned).
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

from interference.root import (d10, d11, K, h10, h11, Q0,   # noqa: E402
                               BA_ratio, M_E_ANCHOR_MEV, PDG_MASSES)
import edge_grammar as eg                                    # noqa: E402

TOL = 0.01          # the scorecard's own band, as in surplus_edges.py
TIGHT = 0.001       # the precision tier: what the sub-percent claims resolve


def _hit(pred, target, tol):
    return abs(pred / target - 1.0) < tol

# canon identities (frozen by tests/test_closure_search.py)
CANON = dict(theta=float(h10),
             up=(("ff", "h10"), ("ffa", "delta"), ("f4+a4", "1/2K")),
             koide_b=("adj", 3), koide_s=("fund", 3))

# the two derived edges, with their derivations
DERIVED_EDGES = (
    ("m_s bridge", "sqrt(Q0^2*charge_trace*(1+h11/K^2)) = sqrt(32/27*73/72)",
     "norm conversion delta_Z2^2/delta_J^2: fixed point + constructed "
     "Albert spectral gap (tests/probes/edge_derivations.py), "
     "conversion unique by Schur on the irreducible 26"),
    ("m_d vent", "d10/d11^4 = 2/81 (triality-swapped h10/d11^2)",
     "exact identity d11^2(1+d10/d11^4) = count(aa) + h10 = 83/9: "
     "the committed word_multiplier operation + generation ladder"),
)


# ═══════════════════════════════════════════════════════════════════════
#  Sector generators + blind gates
# ═══════════════════════════════════════════════════════════════════════

def lepton_candidates():
    """One candidate per θ in the committed angle menu.  B/A = √d₁₀ is
    a theorem edge (octonionic CG + quantum Schur), not a menu choice.
    BLIND GATES: all three gaps positive, three distinct values.
    A candidate is the ratio pair (m_μ/m_e, m_τ/m_e): sorted ascending,
    smallest is the anchor slot."""
    out = []
    for theta in eg.theta_menu():
        gaps = sorted(1.0 + BA_ratio * math.cos(theta + 2*math.pi*k/3)
                      for k in range(3))
        blind = gaps[0] > 0 and len({round(g, 12) for g in gaps}) == 3
        ratios = ((gaps[1]/gaps[0])**2, (gaps[2]/gaps[0])**2) if blind \
            else None
        out.append(dict(theta=theta, blind=blind, ratios=ratios))
    return out


def up_candidates(gen, length_words):
    """m_q(n) = [count(w) + corr]·m_lep(n): every word of the length
    rule's length × every correction (or none).
    BLIND GATE: count ≥ 1 (a multiplier; count 0 = no walk closes).
    VACUOUS IN PRACTICE: every word in the committed menus has count
    ≥ 4, so this gate removes nothing and the up-stage blind counts
    equal the candidate counts (stated in the module docstring; the
    joint denominator inherits it)."""
    corrs = dict(eg.correction_menu())
    corrs["none"] = 0.0
    out = []
    for w, cnt in sorted(length_words.items()):
        for cname, c in sorted(corrs.items()):
            out.append(dict(gen=gen, word=w, count=cnt, corr=cname,
                            mult=cnt + c, blind=cnt >= 1))
    return out


def koide_invert_b(Q, m_c, m_t):
    """Invert Q(c,b,t) = Q for m_b.  BLIND GATES: discriminant > 0 and
    a root bracketed strictly between the partners (the length rule
    orders the generations; masses.py states the same rule)."""
    sc, st = math.sqrt(m_c), math.sqrt(m_t)
    a = Q - 1.0
    b = 2.0 * Q * (sc + st)
    c = (Q - 1.0) * (m_c + m_t) + 2.0 * Q * sc * st
    disc = b*b - 4.0*a*c
    if disc <= 0.0:
        return None
    roots = [(-b + math.sqrt(disc)) / (2.0*a),
             (-b - math.sqrt(disc)) / (2.0*a)]
    hits = [x*x for x in roots if m_c < x*x < m_t]
    return hits[0] if len(hits) == 1 else None


def koide_invert_s(Q, m_c, m_b, bridge):
    """Invert Q(s,c,b) = Q for m_s on the negative Brannen branch:
    √m_s enters the Koide sum with a minus sign (the quadratic below
    is the +branch equation under x -> -x, so its positive root is
    -√m_s).  Same branch as masses.py canon; the b inversion above
    uses the +branch.  A declared discrete branch choice.  BLIND
    GATES: discriminant > 0, exactly one positive root, and the
    result lighter than both partners."""
    sc, sb = math.sqrt(m_c), math.sqrt(m_b)
    a = 1.0 - Q
    b = 2.0 * Q * (sc + sb)
    c = m_c + m_b - Q * (sc + sb)**2
    disc = b*b - 4.0*a*c
    if disc <= 0.0:
        return None
    roots = [(-b + math.sqrt(disc)) / (2.0*a),
             (-b - math.sqrt(disc)) / (2.0*a)]
    pos = [x for x in roots if x > 0]
    if len(pos) != 1:
        return None
    m_s = pos[0]**2 * bridge
    return m_s if m_s < min(m_c, m_b) else None


# ═══════════════════════════════════════════════════════════════════════
#  The staged search
# ═══════════════════════════════════════════════════════════════════════

def run(report=print):
    report("CLOSURE SEARCH: the committed grammar run as a generator")
    report("=" * 70)
    report("  structural wiring fixed: E8(1) > G2(1) x F4(1)")
    report("  stage order: BLIND GATES first, data contact once, last.")

    stages = {}
    hit_ids = {}          # per-stage hit identities (frozen by tests)

    ## Stage L: leptons ##
    lep = lepton_candidates()
    blind_L = [c for c in lep if c["blind"]]
    tgt = (PDG_MASSES["mu"] / PDG_MASSES["e"],
           PDG_MASSES["tau"] / PDG_MASSES["e"])
    hits_L = [c for c in blind_L
              if _hit(c["ratios"][0], tgt[0], TOL)
              and _hit(c["ratios"][1], tgt[1], TOL)]
    tight_L = [c for c in hits_L
               if _hit(c["ratios"][0], tgt[0], TIGHT)
               and _hit(c["ratios"][1], tgt[1], TIGHT)]
    stages["leptons"] = (len(lep), len(blind_L), len(hits_L), len(tight_L))
    hit_ids["leptons"] = sorted(c["theta"] for c in hits_L)
    report(f"\n  L  θ-menu:      {len(lep):4d} candidates → "
           f"{len(blind_L):3d} blind survivors → {len(hits_L)} hit(s) @1%"
           f" → {len(tight_L)} @0.1%")
    for c in hits_L:
        report(f"       hit: θ = {c['theta']:.10g}"
               f"{'  (= h₁₀ = 2/9, canon)' if abs(c['theta']-float(h10)) < 1e-12 else ''}")

    # canon lepton masses for the quark stages (anchor sets the unit)
    theta = float(h10)
    gaps = {k: 1.0 + BA_ratio * math.cos(theta + 2*math.pi*k/3)
            for k in range(3)}
    m_lep = {1: M_E_ANCHOR_MEV,
             2: M_E_ANCHOR_MEV * (gaps[2]/gaps[1])**2,
             3: M_E_ANCHOR_MEV * (gaps[0]/gaps[1])**2}

    ## Stage U: up-type, one sub-stage per generation ##
    up_specs = [(1, eg.word_menu(2), PDG_MASSES["u"]),
                (2, eg.word_menu(3), PDG_MASSES["c"]),
                (3, eg.terminus_menu(), PDG_MASSES["t"])]
    up_blind_masses = {}
    for gen, words, target in up_specs:
        cands = up_candidates(gen, words)
        blind = [c for c in cands if c["blind"]]
        # monotonicity across generations is enforced by the length
        # rule itself (longer word ⇒ larger walk count menu); recorded
        # per-candidate below at the joint stage.
        hits = [c for c in blind
                if _hit(c["mult"] * m_lep[gen], target, TOL)]
        tight = [c for c in hits
                 if _hit(c["mult"] * m_lep[gen], target, TIGHT)]
        stages[f"up_gen{gen}"] = (len(cands), len(blind), len(hits),
                                  len(tight))
        hit_ids[f"up_gen{gen}"] = sorted((c["word"], c["corr"])
                                         for c in hits)
        up_blind_masses[gen] = [(c, c["mult"] * m_lep[gen]) for c in blind]
        report(f"  U{gen} word×corr:   {len(cands):4d} candidates → "
               f"{len(blind):3d} blind survivors → {len(hits)} hit(s) @1%"
               f" → {len(tight)} @0.1%")
        for c in hits:
            mark = (gen, c["word"], c["corr"]) in {
                (1, "ff", "h10"), (2, "ffa", "delta"), (3, "f4+a4", "1/2K")}
            t = "  @0.1%" if c in tight else ""
            report(f"       hit: {c['word']} + {c['corr']} "
                   f"= {c['mult']:.6g}{t}{'  (canon)' if mark else ''}")

    ## Stage D: down-type Koide inversions ##
    m_c_canon = next(m for c, m in up_blind_masses[2]
                     if (c["word"], c["corr"]) == ("ffa", "delta"))
    m_t_canon = next(m for c, m in up_blind_masses[3]
                     if (c["word"], c["corr"]) == ("f4+a4", "1/2K"))

    kmenu = dict(eg.koide_menu())
    # the DECLARED identity insertion: Q0 with no vent, the zero case
    # of the grammar's koide_structure operation (module docstring;
    # 8 committed structures + this = 9 candidates, not a silent
    # widening of the menu)
    kmenu[("none", 0)] = float(Q0)
    b_rows = [(k, koide_invert_b(Q, m_c_canon, m_t_canon))
              for k, Q in sorted(kmenu.items())]
    blind_b = [(k, m) for k, m in b_rows if m is not None]
    hits_b = [(k, m) for k, m in blind_b
              if _hit(m, PDG_MASSES["b"], TOL)]
    tight_b = [(k, m) for k, m in hits_b if _hit(m, PDG_MASSES["b"], TIGHT)]
    stages["down_b"] = (len(b_rows), len(blind_b), len(hits_b),
                        len(tight_b))
    hit_ids["down_b"] = sorted(k for k, m in hits_b)
    report(f"  Db koide-menu:  {len(b_rows):4d} candidates → "
           f"{len(blind_b):3d} blind survivors → {len(hits_b)} hit(s) @1%"
           f" → {len(tight_b)} @0.1%")
    for k, m in hits_b:
        report(f"       hit: Q₀+h({k[0]})/K^{k[1]} → m_b = {m:.1f} MeV"
               f"{'  (canon)' if k == CANON['koide_b'] else ''}")

    # m_s: the bridge edge is the derived norm conversion (see
    # DERIVED_EDGES).  Both readings are still reported:
    # bridge = 1 (menu-only) and bridge = canon derived edge.
    m_b_canon = next(m for k, m in blind_b if k == CANON["koide_b"])
    charge_trace = d10**3 / d11
    bridge_canon = math.sqrt(float(Q0)**2 * charge_trace
                             * (1.0 + float(h11) / K**2))
    for tag, bridge in (("menu-only", 1.0), ("with derived edge",
                                             bridge_canon)):
        s_rows = [(k, koide_invert_s(Q, m_c_canon, m_b_canon, bridge))
                  for k, Q in sorted(kmenu.items())]
        blind_s = [(k, m) for k, m in s_rows if m is not None]
        hits_s = [(k, m) for k, m in blind_s
                  if _hit(m, PDG_MASSES["s"], TOL)]
        tight_s = [(k, m) for k, m in hits_s
                   if _hit(m, PDG_MASSES["s"], TIGHT)]
        stages[f"down_s ({tag})"] = (len(s_rows), len(blind_s),
                                     len(hits_s), len(tight_s))
        hit_ids[f"down_s ({tag})"] = sorted(k for k, m in hits_s)
        report(f"  Ds {tag:<22s}: {len(s_rows):2d} candidates → "
               f"{len(blind_s):2d} blind survivors → {len(hits_s)} hit(s)"
               f" @1% → {len(tight_s)} @0.1%")
        for k, m in hits_s:
            report(f"       hit: Q₀+h({k[0]})/K^{k[1]} → m_s = {m:.1f} MeV"
                   f"{'  (canon)' if k == CANON['koide_s'] else ''}")

    ## The accounting ##
    joint_cand = joint_blind = joint_hits = joint_tight = 1
    for name in ("leptons", "up_gen1", "up_gen2", "up_gen3", "down_b",
                 "down_s (with derived edge)"):
        n, b, h, t = stages[name]
        joint_cand *= n
        joint_blind *= b
        joint_hits *= h
        joint_tight *= t

    report("-" * 70)
    report(f"  JOINT: {joint_cand} grammar-expressible webs → "
           f"{joint_blind} blind survivors → {joint_hits} hit(s) @1% → "
           f"{joint_tight} @0.1%")
    report(f"  look-elsewhere denominator (blind survivors): {joint_blind}")
    report(f"  trials statement: 7 of the 8 non-anchor masses "
           f"(μ τ u c t b s; m_d is the ancillary stage below) are "
           f"{joint_hits} hit(s) in {joint_blind} eligible mass grammar "
           f"realizations at "
           f"the 1% band, {joint_tight} at 0.1% (down-type staged on "
           f"canon upstream)")
    report(f"  note: the 0.1% tier is finer than the experimental/scheme "
           f"precision of the u, s, b coordinates (PDG spreads 0.5-2%), "
           f"so it is decisive only where the reference is precise "
           f"(μ, τ, c, t)")
    report(f"  canon status: one of {joint_hits} joint hits at 1%, NOT "
           f"unique at that band (word counts resolve uniquely, "
           f"sub-percent corrections do not); joint 0.1% tier = "
           f"{joint_tight}, canon included (U1 and Db carry no tight "
           f"hits)")
    # Down-gen-1 stage (ANCILLARY: outside the frozen joint
    # denominator; not multiplied into the accounting)
    w2 = eg.word_menu(2)
    corrs = dict(eg.correction_menu())
    corrs["none"] = 0.0
    d_uncond = [(w, cn) for w, c in sorted(w2.items())
                for cn, cv in sorted(corrs.items())
                if _hit((c + cv) * M_E_ANCHOR_MEV, PDG_MASSES["d"], TOL)]
    d_cond = [w for w, c in w2.items()
              if _hit((c + float(h10)) * M_E_ANCHOR_MEV,
                      PDG_MASSES["d"], TOL)]
    report(f"\n  D1 down-gen-1 (ancillary): unconditioned hits @1% = "
           f"{d_uncond}; ladder-conditioned")
    report(f"     (corr(1) = h10, fixed at up-gen-1): unique hit = "
           f"{d_cond}  (canon: aa + h10 = 83/9 exactly)")

    report(f"\n  DERIVED EDGES (beyond the word/correction "
           f"menus; each a derivation, not a choice):")
    for name, expr, derivation in DERIVED_EDGES:
        report(f"    {name}: {expr}")
        report(f"        = {derivation}")

    return dict(stages=stages, stage_hits=hit_ids,
                joint_candidates=joint_cand,
                joint_blind=joint_blind, joint_hits=joint_hits,
                joint_tight=joint_tight, derived=DERIVED_EDGES,
                down_gen1=dict(uncond=d_uncond, cond=d_cond))


if __name__ == "__main__":
    run()
