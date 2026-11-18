"""
test_emergence_ladder.py: freeze the record of couplings per rung.

What is frozen:
  1. placement: every web node sits on exactly one rung, every channel
     maps to a rung, and every source a rung is read from exists;
  2. the common mode: Sigma_vent(UV) is fed by exactly the four rungs
     that carry degrees of freedom, its rows sum to Sigma, and without
     the bridge row the sum is -10/3 (the sign gate as entered in the ledger);
  3. one path for the deepest backreaction: every depth 4 echo couples
     to the depth 3 e<->q loop;
  4. sparse rungs: u, d, c, t are read from their lepton by ratio, and the rungs
     without a web node are exactly the algebra, the substrate, hadrons,
     atoms and the cosmological constant.
"""

import pytest

import emergence_ladder


@pytest.fixture(scope="module")
def ladder():
    return emergence_ladder.run(report=lambda *a: None)


def test_every_node_channel_and_source_is_placed(ladder):
    assert ladder["unplaced_nodes"] == []
    assert ladder["unplaced_channels"] == []
    assert ladder["doubly_placed"] == []
    assert ladder["phantom_nodes"] == []
    assert ladder["missing_sources"] == []


def test_common_mode_is_fed_by_every_rung_with_degrees_of_freedom(ladder):
    assert ladder["feeders"]["Sigma_vent(UV)"] == {
        "bridge", "forces", "leptons", "quarks"}
    assert ladder["sigma_sum"] == pytest.approx(ladder["sigma"], abs=1e-12)
    assert ladder["sigma"] == pytest.approx(27.0, abs=1e-4)
    assert ladder["sigma_without_bridge"] == pytest.approx(-10 / 3, abs=1e-12)


def test_depth_four_echoes_read_the_charged_loop(ladder):
    deep = [e for e in ladder["edges"] if e["depth"] == 4]
    assert len(deep) == 7 and all(e["parent"] == "eq_loop" for e in deep)
    readers = [e for e in deep if e["node"] != "common_mode"]
    assert len(readers) == 5 and all(e["path"].startswith("read:eq_loop") for e in readers)
    assert sorted(e["path"] for e in deep if e["node"] == "common_mode") == ["Z3(q=1)", "Z3(q=2)"]


def test_quarks_read_their_lepton_and_the_web_stops_where_stated(ladder):
    reads = {e["node"]: e["path"] for e in ladder["edges"]
             if e["kind"] == "ratio" and e["node"] in ("m_u", "m_d", "m_c", "m_t")}
    assert reads == {"m_u": "reads m_e", "m_d": "reads m_e",
                     "m_c": "reads m_mu", "m_t": "reads m_tau"}
    bare = {r["key"] for r in ladder["rungs"] if not r["nodes"]}
    assert bare == {"substrate", "algebra", "hadrons", "atoms", "lambda"}


def test_gravity_rung_lists_its_gates_and_readers(ladder):
    """The lane's six gates sit on the gravity rung, with verdict, tier
    and witnesses taken from the engine. The nodes of other rungs that
    depend on the common mode are the masses the anchor transmutes."""
    gates = ladder["gates"]
    assert len(gates) == 6 and all(v[0] == "PASS" and v[1] == "derived" and v[2] >= 4 for v in gates.values())
    assert ladder["gravity_readers"] == ["m_e", "m_mu", "m_tau", "v_EW_GeV"]
