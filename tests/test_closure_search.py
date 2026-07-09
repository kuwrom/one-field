"""
test_closure_search.py: freeze the forward search (D1 run as a
generator, blind gates first, data contact once).

What is frozen:
  1. the stage table (candidates → blind survivors → hits @1% → @0.1%),
     so the look-elsewhere denominator cannot drift silently;
  2. canon's membership in the hit set at every stage;
  3. the two derived edges (m_s norm conversion, m_d on-menu
     identity) with their derivations, and the down-gen-1 stage's
     unique conditioned hit.  If the committed menus are ever
     widened, this test fails and the freeze is re-committed:
     widening the grammar is a committed act, exactly as
     edge_grammar.py stipulates.
"""

import pytest

import closure_search


@pytest.fixture(scope="module")
def search():
    lines = []
    out = closure_search.run(report=lines.append)
    return out, "\n".join(lines)


def test_stage_table_frozen(search):
    out, text = search
    assert out["stages"] == {
        "leptons":                    (14, 1, 1, 1),
        "up_gen1":                    (28, 28, 1, 0),
        "up_gen2":                    (56, 56, 18, 3),
        "up_gen3":                    (21, 21, 7, 4),
        "down_b":                     (9, 8, 1, 0),
        "down_s (menu-only)":         (9, 9, 0, 0),
        "down_s (with derived edge)": (9, 9, 4, 0),
    }
    assert out["joint_candidates"] == 37340352
    assert out["joint_blind"] == 2370816
    assert out["joint_hits"] == 504
    assert out["joint_tight"] == 0     # tighter than u/s/b data precision
    assert "structural wiring fixed: E8(1) > G2(1) x F4(1)" in text
    assert "504 hit(s) in 2370816 eligible mass grammar realizations" in text


def test_canon_is_a_hit_at_every_stage(search):
    """Canon is a hit at every stage, asserted on hit identities."""
    out, _ = search
    h = out["stage_hits"]
    assert any(abs(t - closure_search.CANON["theta"]) < 1e-12
               for t in h["leptons"])
    for gen, ident in enumerate(closure_search.CANON["up"], start=1):
        assert ident in h[f"up_gen{gen}"]
    assert closure_search.CANON["koide_b"] in h["down_b"]
    assert closure_search.CANON["koide_s"] in \
        h["down_s (with derived edge)"]
    # canon m_s requires the derived bridge; it is not a menu-only hit
    assert closure_search.CANON["koide_s"] not in h["down_s (menu-only)"]


def test_theta_is_unique_blind_survivor(search):
    """The strongest single result of the forward search: of the 14
    committed angles, exactly ONE survives the blind positivity gate
    at all - and it is h₁₀ = 2/9, and it hits both lepton ratios."""
    out, text = search
    # the stage table is frozen in test_stage_table_frozen; this
    # asserts the unique survivor's identity
    assert len(out["stage_hits"]["leptons"]) == 1
    assert abs(out["stage_hits"]["leptons"][0] - 2.0 / 9.0) < 1e-12
    assert "θ = 0.2222222222" in text


def test_cross_derived_edges(search):
    """The two edges beyond the word/correction menus are derivations
    with stated witnesses (m_s: Schur-unique norm conversion; m_d:
    the exact on-menu identity), and the report only down-gen-1 stage
    shows canon (aa + h10) as the unique ladder-conditioned hit."""
    out, _ = search
    assert tuple(name for name, _, _ in out["derived"]) == (
        "m_s bridge", "m_d vent")
    assert all(deriv for _, _, deriv in out["derived"])
    assert ("aa", "h10") in out["down_gen1"]["uncond"]
    assert out["down_gen1"]["cond"] == ["aa"]


def test_word_counts_unique_at_1pct(search):
    """At the 1% band the data selects the word
    COUNT uniquely at each up-type stage; the sub-percent corrections
    are what the hits differ in.  U2: all 18 hits carry count 12.
    U3: all 7 hits carry count 97."""
    import edge_grammar as eg
    out, _ = search
    w3, term = eg.word_menu(3), eg.terminus_menu()
    c2 = {w3[w] for w, _ in out["stage_hits"]["up_gen2"]}
    c3 = {term[w] for w, _ in out["stage_hits"]["up_gen3"]}
    assert c2 == {12}
    assert c3 == {97}
