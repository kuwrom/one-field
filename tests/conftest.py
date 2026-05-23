"""Shared fixtures for the scorecard test suite.

The full derivation chain is run once per test session, silently.  Every
test then asserts against the cached results, so the test suite is fast
(<1s total) and the layer-by-layer console output never pollutes pytest.
"""

import contextlib
import io

import pytest

from E8 import (
    algebra,
    alpha_s,
    baryogenesis,
    ckm,
    gravity,
    higgs,
    leptons,
    neutrinos,
    octonions,
    pmns,
    quarks,
    scale,
    wzw,
)


@pytest.fixture(scope="session")
def res():
    """Run the full E8 derivation chain once, silently."""
    with contextlib.redirect_stdout(io.StringIO()):
        alg = algebra.derive()
        scl = scale.derive(alg)
        lep = leptons.derive(scl)
        qrk = quarks.derive(alg, scl, lep)
        wz = wzw.derive()
        oc = octonions.derive(wz)
        ck = ckm.derive(wz)
        pm = pmns.derive(wz, lep, ck)
        nu = neutrinos.derive(alg)
        als = alpha_s.derive(alg, scl, qrk)
        hg = higgs.derive(scl, qrk, als, lep, alg)
        gr = gravity.derive(alg, scl)
        br = baryogenesis.derive(alg, pm)
    return dict(
        alg=alg, scl=scl, lep=lep, qrk=qrk, wz=wz, oc=oc, ck=ck,
        pm=pm, nu=nu, als=als, hg=hg, gr=gr, br=br,
    )
