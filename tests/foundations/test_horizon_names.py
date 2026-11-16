"""Freezes tests/probes/horizon_names.py: the bridge as a multiplicity of the E8(1) unit."""


def test_the_bridge_is_a_multiplicity_of_the_unit():
    from probes import horizon_names
    r = horizon_names.run(report=lambda *a: None)
    assert r["det"] == 1 and r["even"] and r["weight_one"] == 248
    assert r["split"] == (14, 52, 182) and r["zero_divisor"]
