"""Freezes tests/probes/kernel_sign.py: the one list's spin-2 form factor in ball arithmetic."""

import pytest

pytest.importorskip("flint")


def _silent(*args, **kwargs):
    pass


def test_kernel_sign_closed_form_and_limit():
    """The closed form equals the node kernel to 1e-6 for y <= 10 and its
    two integral forms agree to 1e-8. At y = 5e4 the node kernel is
    negative where the form factor is 8.1241. Certified positive on
    [1e-3, 0.1] and [1e3, 1e12]. The limit is 65/8, positive for more
    than 117 scalars, and the label's xi on the amount puts a zero
    between y = 6.5 and 6.7."""
    from probes import kernel_sign as ks

    r = ks.run(report=_silent, full=False)
    assert r["node_agreement"] < 1e-6 and r["form_agreement"] < 1e-8
    assert r["node_at_5e4"] < 0 and r["phi_at_5e4"] == pytest.approx(8.1241, abs=1e-4)
    assert r["least"] > 0 and r["limit"] == 65 / 8 and r["scalars_needed"] == 117
    assert r["phi_110"] < 0 and r["label_zero"] == (6.5, 6.7)


@pytest.mark.slow
def test_kernel_sign_is_positive_on_the_whole_axis():
    """Every y in [1e-3, 1e12] lies in an interval whose rigorous
    enclosure of the form factor is positive."""
    from probes import kernel_sign as ks

    r = ks.run(report=_silent, full=True)
    assert r["certified"] == ((1e-3, 1e3), (1e3, 1e12)) and r["least"] > 0.003
