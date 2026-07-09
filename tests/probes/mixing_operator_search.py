"""Exact mixing-sector path-space audit, before mixing data are consulted.

The existing D6* closure words are polynomials in commuting fusion
operators, so their scalar path sums cannot produce a relative
mass-frame diagonalizer.  The exact D^(6)* cells change that conclusion
on the open path space: their braid generators do not commute, satisfy
Yang--Baxter, and at depth four contain two three-dimensional carriers.
What remains is the dynamical assignment of those carriers to sectors.

The Albert carrier does split explicitly as 26 = 2 + 8 + 8 + 8, but the
current bridge acts identically on the irreducible 26 and the depth-4
transport reads the closed path g^-1 g = 1.  Neither operation selects sector frames.

This probe records those facts only.  It does not scan candidate mixing
maps or compare with CKM/PMNS data.
"""

import os
import sys
from fractions import Fraction

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
for _path in (os.path.dirname(_here),
              os.path.join(os.path.dirname(os.path.dirname(_here)),
                           "interference"),
              os.path.join(_here, "interference")):
    if _path not in sys.path:
        sys.path.insert(0, _path)
import _bootstrap  # noqa: E402,F401

from words import N_FUND  # noqa: E402
from mixing import _compute_d6_boundary  # noqa: E402

try:  # noqa: E402
    from probes.wiring_theorem_probe import albert_unital_response
except ImportError:  # direct execution from tests/probes
    from wiring_theorem_probe import albert_unital_response


def _spin8_branching_projectors():
    """Orthogonal coordinate projectors for J3(O)=1+2+8+8+8."""
    unit = np.zeros((27, 27))
    unit[:3, :3] = np.ones((3, 3)) / 3.0
    diagonal = np.zeros((27, 27))
    diagonal[:3, :3] = np.eye(3) - unit[:3, :3]
    arms = []
    for start in (3, 11, 19):
        projector = np.zeros((27, 27))
        projector[start:start + 8, start:start + 8] = np.eye(8)
        arms.append(projector)
    projectors = (unit, diagonal, *arms)
    assert np.allclose(sum(projectors), np.eye(27))
    assert all(np.allclose(p @ p, p) for p in projectors)
    assert all(np.allclose(p @ q, 0.0)
               for i, p in enumerate(projectors)
               for j, q in enumerate(projectors) if i != j)
    ranks = tuple(int(round(np.trace(p))) for p in projectors)
    assert ranks == (1, 2, 8, 8, 8)
    return ranks


def _spectra_do_not_fix_relative_frames():
    """Exhibit two operator pairs with the same spectra but different mixing."""
    up = np.diag([1.0, 2.0, 4.0])
    down = np.diag([3.0, 5.0, 9.0])
    angle = 0.37
    c, s = np.cos(angle), np.sin(angle)
    rotation = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    rotated_down = rotation @ down @ rotation.T

    assert np.allclose(np.linalg.eigvalsh(down),
                       np.linalg.eigvalsh(rotated_down))
    identity_relative_frame = np.eye(3)
    rotated_relative_frame = rotation
    assert not np.allclose(identity_relative_frame, rotated_relative_frame)
    return {
        "same_up_spectrum": tuple(np.linalg.eigvalsh(up)),
        "same_down_spectrum": tuple(np.linalg.eigvalsh(down)),
        "relative_frames_differ": True,
    }


def run(report=print):
    report("MIXING SECTOR OPEN-PATH AUDIT: no mixing data")
    report("=" * 58)

    n10 = np.asarray(N_FUND, dtype=float)
    n01 = n10.T
    n11 = n10 @ n01 - np.eye(6)
    generators = {"n10": n10, "n01": n01, "n11": n11}
    commutators = {
        f"[{left},{right}]": float(np.linalg.norm(a @ b - b @ a))
        for left, a in generators.items()
        for right, b in generators.items() if left < right
    }
    assert set(commutators.values()) == {0.0}
    report("  D6* fusion generators commute: closure-word operators share")
    report("  an eigenbasis after open paths have been summed.")

    d6 = _compute_d6_boundary()
    assert d6["type_i_residual"] < 1e-12
    assert d6["connection_unitarity_residual"] < 1e-12
    assert d6["braid_residual"] < 1e-12
    assert d6["braid_commutator"] > 1e-12
    assert sorted(d6["depth4_carrier_dimensions"]) == [2, 3, 3]
    report("  Exact D6* cells restore the open information: the braid")
    report("  connection is unitary/Yang--Baxter and noncommuting.")
    report("  Depth four splits as 3+3+2, exposing two generation-sized")
    report("  carriers before the scalar path sum is taken.")

    ranks = _spin8_branching_projectors()
    report("  Albert branching projectors: J3(O)=1+2+8+8+8;")
    report("  therefore P26=2+8+8+8, particle labels follow from the mass hierarchy (masses.py).")

    albert = albert_unital_response()
    assert albert["junction"] == 1
    assert albert["response"] == Fraction(27, 2)
    report("  The unital closed transport g^-1 g=1 reads only")
    report("  its trace response 27/2; the relative endpoint map cancels.")

    spectral = _spectra_do_not_fix_relative_frames()
    report("  Scalar mass spectra do not determine eigenvectors: the same")
    report("  up/down spectra admit both aligned and rotated relative frames.")

    report("  The mixing sector now requires a dynamical sector assignment for the two")
    report("  three-dimensional carriers, before tracing.")
    return {
        "d6_commutators": commutators,
        "d6_connection_magnitude": d6["connection_magnitude"],
        "d6_braid_residual": d6["braid_residual"],
        "d6_braid_commutator": d6["braid_commutator"],
        "depth4_carrier_dimensions": d6["depth4_carrier_dimensions"],
        "j3_branching_projector_ranks": ranks,
        "p26_branching_projector_ranks": ranks[1:],
        "albert_junction": albert["junction"],
        "albert_trace_response": albert["response"],
        "spectral_nonuniqueness": spectral,
        "selected_open_transport": None,
    }


if __name__ == "__main__":
    run()
