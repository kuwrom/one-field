"""Freezes tests/probes/unnamed_metric.py: with no names, a metric free
symmetric tensor is built from scalars alone; 2-forms need an oriented
triple (Urbantke), 3-forms give nothing."""


def _silent(*args, **kwargs):
    pass


def test_unnamed_metric_is_scalar_built():
    from probes import unnamed_metric as um

    out = um.run(report=_silent)
    for k, (s, a) in out["two_forms"].items():
        if k == "Urbantke eps_ijk":
            assert s > 1e-2 and a < 1e-12
        else:
            assert s < 1e-12
    assert out["collapse"] < 1e-12
    assert out["three_forms"] < 1e-12
    assert out["scalar_rank"] == 4
