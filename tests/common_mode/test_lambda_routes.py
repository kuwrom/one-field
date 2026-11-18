"""Freezes the Lambda routes A, B and C.  Route B
(lambda_unit.py): the metric's unit cannot absorb V0, so V0 = -rho_vac
is a tuning, and off it the vacuum is de Sitter with no graviton mass.
Route A (composite_volume.py): the carrier's volume element is nonzero
only where its metric is Euclidean and vanishes in expectation in every
vacuum with a symmetric gradient correlator."""


def _silent(*args, **kwargs):
    pass


def test_lambda_route_b_the_vacuum_constant_is_a_tuning():
    from probes import lambda_unit as lu

    out = lu.run(report=_silent)
    assert all(abs(v - 1) < 1e-12 for v in out["b1"]["ratio"].values())
    assert out["b2"][-1.0] == 0.0 and out["b2"][-0.999] * out["b2"][-1.001] < 0
    assert out["b3"]["second_order"] < 1e-12


def test_lambda_route_a_the_carrier_volume_element_is_excluded():
    from probes import composite_volume as cv

    out = cv.run(report=_silent)
    assert abs(out["a1"]["weight_one"] - 1) < 1e-10
    assert out["a3"]["min_eig_when_phi"] > 0
    assert out["a3"]["phi_when_rank3"] < 1e-10
    assert out["a4"]["lorentz"] < 1e-12 and out["a4"]["real_symmetric"] < 1e-10
    assert out["a4"]["frame_fixing"] > 1e-3


def test_lambda_route_c_the_grain_is_a_unit():
    """Route C (composite_q.py): the one candidate that passes the form
    tests is the regulator's grain.  It is the empty field's only scale.
    A q-reading of it would leave the dimension two moment."""
    from probes import composite_q as cq

    out = cq.run(report=_silent)
    assert all(d < 1e-12 for d in out["c2"].values())
    for r in out["c3"].values():
        assert abs(r["numeric"] - r["closed"]) < 1e-6 * r["quartic"] + 1e-13


def test_lambda_d3_record_count_overshoots_the_window():
    """everpresent_grain.py: a record count at the anchor's grain gives a
    fluctuation amplitude above the 0.01-0.02 window in every convention."""
    from probes import everpresent_grain as eg

    out = eg.run(report=_silent)
    assert all(b > 0.02 for b in out["reduced"].values())
    assert all(b > 0.02 for b in out["unreduced"].values())


def test_lambda_lane_the_trace_free_form_is_the_one_that_closes():
    """lambda_lane.py: under the axiom (a constant in the potential weighs
    no history) the vacuum term's density has three options: zero, a
    closed form built from the field, or fixed.  The first has no lump.
    The second has no Lorentzian vacuum.  The third is the traceless
    form, which keeps the full form's saddle structure."""
    from probes import lambda_lane as ll

    out = ll.run(report=_silent)
    assert out["l1"]["full_second"] < 1e-5 and out["l1"]["trace_free"] < 1e-6
    assert out["l2"]["min_eig"] > 0 and out["l2"]["lorentz"] < 1e-9
    assert out["l5"]["passing"] == 4 and out["l5"]["kin"] < 1e-6
    assert out["l2b"]["lorentz"] < 1e-8 * max(1.0, out["l2b"]["framed"])
    assert out["l2b"]["framed"] > 1e-3
    assert out["l6"]["trace_free"] < 1e-6 and out["l6"]["full"] < 1e-6


def test_lambda_scope_reduces_to_the_lanes_premises():
    """lambda_scope.py: a null Lagrangian is a divergence built from the
    field plus its value on the empty field. The divergence averages to zero in a
    translation invariant vacuum and not on a ruler profile."""
    from probes import lambda_scope as ls

    out = ls.run(report=_silent)
    assert out["s2"]["div_vs_phi"] < 1e-6 and out["s2"]["phi_at_zero"] == 0.0
    assert out["s3"]["periodic_average"] < 1e-12 and abs(out["s3"]["ruler_average"]) > 1.0
