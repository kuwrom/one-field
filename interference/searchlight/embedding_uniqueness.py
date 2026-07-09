"""
E₈ conformal embedding uniqueness: the coherence gate (gravity paper, Sec. IX).

CLAIM.  Among the rank-8 parents and factor pairs in the invariant table,
E₈(1) ⊃ G₂(1) × F₄(1) is the UNIQUE unordered survivor satisfying all six
coherence gates.  The surviving embedding is independently known; this
numerical scan does not compute branching rules or embedding indices.

  (a) c_coset = c_parent − c_G − c_F = 0  (no coset Virasoro sector),
  (b) h_bridge = h(fund_G) + h(fund_F) = 1  (dimension-one extension),
  (c) dim(parent) = dim(G) + dim(F) + dim(fund_G)·dim(fund_F)
      (every mode assigned: gauge + matter + bridge, no leftovers),
  (d) bridge dimension n_G·n_F > 0  (non-trivial matter sector),
  (e) parent is self-adjoint: fund = adjoint  (self-encoding: the parent
      theory IS its own deepest representation, no external structure),
  (f) both factors have trivial centre  (all quantum numbers emergent.
      a non-trivial centre would pre-inscribe charges onto the vacuum).

Gates (a)-(d) are the numerical conformal embedding conditions.
Gates (e)-(f) are the physics requirements from the framework's
ontology: the substrate is self-present (e), and quantum numbers
crystallise upon breaking, they are not pre-inscribed (f).

WHERE THE CLOSURE ACTUALLY HAPPENS (measured by this scan, stated
up front rather than left to the reader): gate (e) alone eliminates
every non-E₈ parent (A₈, B₈, C₈, D₈ are not self-adjoint), and gate
(f) alone eliminates the one genuine contender, E₈ ⊃ E₇ × A₁, which
passes (a)-(e) in full (c_coset = 0, h_bridge = 1, 133 + 3 + 112 =
248, bridge 56×2 > 0) and fails only on centres (Z₂ × Z₂).  The
conformal gates (a)-(d) narrow the field; the two framework closure gates do
the decisive work.  The scan prints the contender so this is auditable, and
derive() returns it.

SCOPE OF THE SCAN (what "exhaustive" means here): the scan is
exhaustive over the DATA TABLE below, all pairs of simple algebras
with rank(G) + rank(F) ≤ 8 tested against the six gates by their
numerical invariants.  It does NOT verify that a candidate pair is
actually a subalgebra of the parent (no branching rules, no
embedding indices are computed); a pair passing the numerical gates
without being an embedding would survive the scan.  For the one
survivor this is settled by the literature (E₈ ⊃ G₂ × F₄ is a known
maximal conformal embedding), and the scan's role is to show that
no OTHER candidate even passes the numerical gates.  Isomorphic
low-rank duplicates in the table (D₃ ≅ A₃, C₂ ≅ B₂) are skipped in
the scan so the candidate count is over distinct algebras.

═══════════════════════════════════════════════════════════════════════
THE FIVE RANK-8 CANDIDATES
═══════════════════════════════════════════════════════════════════════

At rank 8 there are exactly five simple Lie algebras:
    A₈ = SU(9),  B₈ = SO(17),  C₈ = Sp(16),  D₈ = SO(16),  E₈

For each, we enumerate all sub-pairs (G, F) of table entries at level
1 whose central charges saturate the parent's (c_coset = 0).  That
condition alone does not close the field (the scan finds a hundred
such pairs); it is the conjunction of all six gates that closes on
one embedding up to ordering, with gates (e) and (f) doing the
decisive work (see above).

═══════════════════════════════════════════════════════════════════════
WHY THESE GATES MATTER
═══════════════════════════════════════════════════════════════════════

c_coset = 0:  in the unitary embedding, the coset between the parent and
  sub-pair has no independent Virasoro sector. This identity alone does
  not determine the heat-kernel spin assignment; the separate
  non-adjoint representation routing assigns the bridge as a scalar-like
  coherence channel rather than a dynamical gauge field.

h_bridge = 1:  the bridge is a dimension-one extension current; that fact
  alone is not a proof of an exactly marginal (1,1) deformation.  Under the
  framework's single normalized-current definition its factor is 1/(2π),
  giving α(0) and Ω_DM/Ω_b = 2π − 1.  In the unitary embedding,
  c_coset = 0 removes an independent coset Virasoro sector; registered Web
  edges control later responses.  Any h ≠ 1 breaks both predictions.

Reference:
    Kahsay, Kibrom Kidane (2026). Trilogy.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

from fractions import Fraction


# ═══════════════════════════════════════════════════════════════════════
#  Lie algebra data (dim, h∨, rank, fund_dim, h_fund at level 1)
# ═══════════════════════════════════════════════════════════════════════
#
# For a simple Lie algebra g at level k=1:
#   c(g) = k·dim(g) / (k + h∨(g)) = dim(g) / (1 + h∨(g))
#   h(fund) = C₂(fund) / (k + h∨(g))  where C₂ is the quadratic Casimir
#
# All values from standard Lie algebra tables (Humphreys / Di Francesco).

def _lie_data():
    """Return dict of simple Lie algebras up to rank 8 with their data.

    Centre orders (|Z(G)|):
        A_n: Z_{n+1},  B_n: Z_2,  C_n: Z_2,
        D_n: Z_4 (n odd) or Z_2×Z_2 (n even),
        G₂: trivial,  F₄: trivial,  E₆: Z_3,  E₇: Z_2,  E₈: trivial.
    """
    algebras = {}

    # A_n = SU(n+1): dim = n(n+2), h∨ = n+1, fund_dim = n+1,
    #   C₂(fund) = n(n+2)/(2(n+1)), centre = Z_{n+1}
    for n in range(1, 9):
        name = f"A{n}"
        dim_g = n * (n + 2)
        hv = n + 1
        fund = n + 1
        C2_fund = Fraction(n * (n + 2), 2 * (n + 1))
        algebras[name] = {
            'dim': dim_g, 'hv': hv, 'rank': n,
            'fund_dim': fund, 'C2_fund': C2_fund,
            'centre_order': n + 1,
            'family': 'A'
        }

    # B_n = SO(2n+1): dim = n(2n+1), h∨ = 2n−1, fund_dim = 2n+1,
    #   C₂(fund) = n  (vector rep), centre = Z_2
    for n in range(2, 9):
        name = f"B{n}"
        dim_g = n * (2 * n + 1)
        hv = 2 * n - 1
        fund = 2 * n + 1
        C2_fund = Fraction(n, 1)
        algebras[name] = {
            'dim': dim_g, 'hv': hv, 'rank': n,
            'fund_dim': fund, 'C2_fund': C2_fund,
            'centre_order': 2,
            'family': 'B'
        }

    # C_n = Sp(2n): dim = n(2n+1), h∨ = n+1, fund_dim = 2n,
    #   C₂(fund) = (2n+1)/4, centre = Z_2
    for n in range(2, 9):
        name = f"C{n}"
        dim_g = n * (2 * n + 1)
        hv = n + 1
        fund = 2 * n
        C2_fund = Fraction(2 * n + 1, 4)
        algebras[name] = {
            'dim': dim_g, 'hv': hv, 'rank': n,
            'fund_dim': fund, 'C2_fund': C2_fund,
            'centre_order': 2,
            'family': 'C'
        }

    # D_n = SO(2n): dim = n(2n−1), h∨ = 2n−2, fund_dim = 2n,
    #   C₂(fund) = (2n−1)/2, centre = Z_4 (n odd) or Z_2×Z_2 (n even)
    for n in range(3, 9):
        name = f"D{n}"
        dim_g = n * (2 * n - 1)
        hv = 2 * (n - 1)
        fund = 2 * n
        C2_fund = Fraction(2 * n - 1, 2)
        algebras[name] = {
            'dim': dim_g, 'hv': hv, 'rank': n,
            'fund_dim': fund, 'C2_fund': C2_fund,
            'centre_order': 4,  # Z₄ (n odd) or Z₂×Z₂ (n even): order 4 both
            'family': 'D'
        }

    # Exceptionals
    algebras['G2'] = {
        'dim': 14, 'hv': 4, 'rank': 2,
        'fund_dim': 7, 'C2_fund': Fraction(2, 1),
        'centre_order': 1,
        'family': 'E'
    }
    algebras['F4'] = {
        'dim': 52, 'hv': 9, 'rank': 4,
        'fund_dim': 26, 'C2_fund': Fraction(6, 1),
        'centre_order': 1,
        'family': 'E'
    }
    algebras['E6'] = {
        'dim': 78, 'hv': 12, 'rank': 6,
        'fund_dim': 27, 'C2_fund': Fraction(26, 3),
        'centre_order': 3,
        'family': 'E'
    }
    algebras['E7'] = {
        'dim': 133, 'hv': 18, 'rank': 7,
        'fund_dim': 56, 'C2_fund': Fraction(57, 4),
        'centre_order': 2,
        'family': 'E'
    }
    algebras['E8'] = {
        'dim': 248, 'hv': 30, 'rank': 8,
        'fund_dim': 248, 'C2_fund': Fraction(30, 1),
        'centre_order': 1,
        'family': 'E'
    }

    return algebras


ALGEBRAS = _lie_data()


def central_charge(name):
    """Central charge c = dim/(1 + h∨) at level k=1."""
    d = ALGEBRAS[name]
    return Fraction(d['dim'], 1 + d['hv'])


def h_fund(name):
    """Conformal weight of the fundamental at level k=1: h = C₂/(k+h∨)."""
    d = ALGEBRAS[name]
    return d['C2_fund'] / (1 + d['hv'])


def derive():
    """
    Exhaustive coherence-gate test over the declared invariant table.

    For each rank-8 algebra P, test all table pairs (G, F) with
    rank(G) + rank(F) ≤ rank(P) for all six conditions.  Pair membership
    as an actual subalgebra is outside this numerical scan.
    """
    print("\n" + "=" * 78)
    print("  E₈ CONFORMAL EMBEDDING UNIQUENESS (coherence gate)")
    print("=" * 78)

    rank8 = {name: data for name, data in ALGEBRAS.items()
             if data['rank'] == 8}
    # skip isomorphic duplicates so the pair count is over DISTINCT
    # algebras: D3 ≅ A3, C2 ≅ B2 (both remain in the table for data
    # completeness; only one representative each enters the scan)
    _ISO_DUPES = {'D3', 'C2'}
    sub_candidates = {name: data for name, data in ALGEBRAS.items()
                      if data['rank'] <= 8 and name not in _ISO_DUPES}

    results = {}
    survivors = []
    near_misses = []     # pairs failing EXACTLY one gate

    for p_name, p_data in sorted(rank8.items(),
                                  key=lambda x: x[1]['dim']):
        c_P = central_charge(p_name)
        print(f"\n  Parent: {p_name}  dim={p_data['dim']}  "
              f"h∨={p_data['hv']}  c={float(c_P):.4f}")

        found = False
        # Try all ordered pairs (G, F) with G ≠ F
        for g_name, g_data in sorted(sub_candidates.items(),
                                      key=lambda x: x[1]['dim']):
            if g_data['rank'] > p_data['rank']:
                continue
            c_G = central_charge(g_name)
            if c_G >= c_P:
                continue

            for f_name, f_data in sorted(sub_candidates.items(),
                                          key=lambda x: x[1]['dim']):
                if f_name == g_name:
                    continue
                if g_data['rank'] + f_data['rank'] > p_data['rank']:
                    continue
                c_F = central_charge(f_name)

                # # Gate (a): c_coset = 0 #
                c_coset = c_P - c_G - c_F
                if c_coset != 0:
                    continue

                # Found a conformal sub-pair with c_coset = 0
                h_G = h_fund(g_name)
                h_F = h_fund(f_name)
                h_br = h_G + h_F

                n_G = g_data['fund_dim']
                n_F = f_data['fund_dim']
                bridge_dim = n_G * n_F
                dim_total = g_data['dim'] + f_data['dim'] + bridge_dim

                # # Gate (e): parent self-adjoint #
                self_adj = (p_data['fund_dim'] == p_data['dim'])

                # # Gate (f): both factors have trivial centre #
                trivial_centres = (g_data['centre_order'] == 1
                                   and f_data['centre_order'] == 1)

                gates = {
                    '(a) c_coset=0':       c_coset == 0,
                    '(b) h_bridge=1':      h_br == 1,
                    '(c) dim_match':       dim_total == p_data['dim'],
                    '(d) bridge>0':        bridge_dim > 0,
                    '(e) self-adjoint':    self_adj,
                    '(f) trivial centres': trivial_centres,
                }
                all_pass = all(gates.values())

                status = "★ ALL GATES PASS" if all_pass else "✗ FAILS"
                print(f"    {g_name} × {f_name}:  c_coset={float(c_coset)}  "
                      f"h_bridge={float(h_br):.4f}  "
                      f"bridge={n_G}×{n_F}={bridge_dim}  "
                      f"dim_check={dim_total}/{p_data['dim']}  "
                      f"{status}")

                if not all_pass:
                    fails = [k for k, v in gates.items() if not v]
                    print(f"      fails: {', '.join(fails)}")

                key = (p_name, g_name, f_name)
                results[key] = gates
                if all_pass:
                    survivors.append(key)
                elif sum(1 for v in gates.values() if not v) == 1:
                    failed_gate = next(k for k, v in gates.items()
                                       if not v)
                    near_misses.append((key, failed_gate))
                found = True

        if not found:
            print("    no conformal sub-pair with c_coset = 0 found")

    # # Verdict (computed from the scan, not asserted as literals) #
    print(f"\n  {'─' * 70}")
    print(f"  SURVIVORS (all six gates): {len(survivors)}")
    for s in survivors:
        print(f"    {s[0]} ⊃ {s[1]} × {s[2]}")
    print(f"  CONTENDERS (fail exactly one gate): {len(near_misses)}")
    for (p, g, f), failed in near_misses:
        print(f"    {p} ⊃ {g} × {f}   fails only {failed}")

    # unordered survivor sets: the scan sees both orderings of a pair
    unordered = {(p, frozenset((g, f))) for p, g, f in survivors}

    assert len(unordered) == 1, (
        f"expected exactly 1 embedding up to ordering, got "
        f"{sorted(unordered)}")
    parent, pair = next(iter(unordered))
    survivor = (parent,) + tuple(sorted(pair, key=lambda n:
                                        ALGEBRAS[n]['dim']))
    # the one genuine contender: E7 × A1 passes (a)-(e), fails only
    # the trivial-centre gate (f).  Surfaced and frozen so the decisive
    # framework closure requirement remains visible.
    nm_keys = {(p, frozenset((g, f))): failed
               for (p, g, f), failed in near_misses}
    assert nm_keys.get(('E8', frozenset(('E7', 'A1')))) \
        == '(f) trivial centres', (
        "E8 > E7 x A1 must appear as the gate-(f) contender")
    print(f"  (both orderings of the same unique embedding)")
    print(f"  UNIQUENESS ESTABLISHED: {parent}(1) ⊃ "
          f"{survivor[1]}(1) × {survivor[2]}(1) through all six gates; "
          f"contender table above")

    return {
        'n_rank8_tested': len(rank8),
        'n_pairs_tested': len(results),
        'n_survivors': len(unordered),
        'survivor': survivor,
        'near_misses': [(p, g, f, failed)
                        for (p, g, f), failed in near_misses],
    }


if __name__ == "__main__":
    derive()
