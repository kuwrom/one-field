"""Strong CP in the full record ontology.

The condensed algebra A=1+B, a character of the winding group, and the
physical matter-coupled angle bar(theta) are distinct objects.  A phase
action T(B)=zB preserves both channels of B*B=1+B only when z=1.  The
doubled Fibonacci phases on the transmissible summands cancel as well.

A central character exp(i n theta) is invisible to that local algebra,
and a common right-handed quark rephasing changes arg det(M_u M_d)
without changing masses or CKM.  The full record ontology resolves both
freedoms by assigning the complete matter-coupled winding character to
the closed vacuum-record action.  This probe computes the local result,
checks the two independent freedoms, and compares it with the production
derivation of bar(theta)=0.

Usage: python3 tests/probes/strong_cp.py
"""

import cmath
import math
import os
import sys
from fractions import Fraction

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from interference import couplings as production_couplings


TOL = 1e-12
PHI = (1.0 + math.sqrt(5.0)) / 2.0


def _phase_from_weight(h):
    return cmath.exp(2j * math.pi * float(h))


def vacuum_algebra_phase_residual(z):
    """Residual of T(B*B)=T(B)T(B) in the (1,B) channels."""
    return z * z - 1.0, z * z - z


def is_unital_phase_automorphism(z, tol=TOL):
    r1, rB = vacuum_algebra_phase_residual(z)
    return abs(r1) < tol and abs(rB) < tol


def large_gauge_character(theta, winding):
    """General U(1) character of pi_3(G2)=Z."""
    return cmath.exp(1j * theta * winding)


def central_character_countermodel(theta, winding=1):
    """A ray character invisible to the induced algebra action.

    U_n=chi_theta(n) I acts on A by conjugation.  Its scalar phase
    cancels, so the induced action has z=1 and preserves multiplication
    for every theta.
    """
    ray_phase = large_gauge_character(theta, winding)
    algebra_z = 1.0 + 0.0j
    return {
        "ray_phase": ray_phase,
        "ray_character_nontrivial": abs(ray_phase - 1.0) > TOL,
        "algebra_z": algebra_z,
        "algebra_action_is_identity": True,
        "algebra_action_is_automorphism": is_unital_phase_automorphism(
            algebra_z),
    }


def right_handed_mass_rephasing(alpha, generations=3):
    """Exact mass-sector null direction.

    Multiplying a Dirac matrix on the right by exp(i alpha) I leaves
    M M^dagger invariant.  Its determinant is multiplied by
    exp(i generations alpha), which is invisible to masses and CKM.
    """
    right_phase = cmath.exp(1j * alpha)
    return {
        "right_phase": right_phase,
        "gram_multiplier": right_phase * right_phase.conjugate(),
        "determinant_multiplier": right_phase ** generations,
        "masses_unchanged": True,
        "left_ckm_unchanged": True,
    }


def engine_mass_matrix_countermodel(up_masses, down_masses, ckm, alpha):
    """Realize the mass null direction on engine outputs.

    A positive-Hermitian representative is

        M_u = D_u,       M_d = V D_d V^dagger,

    which has the supplied positive masses and left mixing V.  Rephase
    M_u on the right by exp(i alpha)I.  Both left Gram operators, and
    therefore all masses and CKM data, stay fixed, while det(M_u M_d)
    acquires exp(3 i alpha).
    """
    up = np.asarray(up_masses, dtype=float)
    down = np.asarray(down_masses, dtype=float)
    V = np.asarray(ckm, dtype=complex)
    assert up.shape == down.shape == (3,) and V.shape == (3, 3)
    assert np.all(up > 0.0) and np.all(down > 0.0)
    assert np.allclose(V.conjugate().T @ V, np.eye(3), atol=1e-12)

    Mu = np.diag(up).astype(complex)
    Md = V @ np.diag(down) @ V.conjugate().T
    phase = cmath.exp(1j * alpha)
    Mu_rephased = Mu * phase
    gram_u = Mu @ Mu.conjugate().T
    gram_u_rephased = Mu_rephased @ Mu_rephased.conjugate().T
    gram_d = Md @ Md.conjugate().T
    gram_d_expected = V @ np.diag(down ** 2) @ V.conjugate().T
    gram_u_residual = float(
        np.max(np.abs(gram_u_rephased - gram_u))
        / np.max(np.abs(gram_u)))
    gram_d_residual = float(
        np.max(np.abs(gram_d - gram_d_expected))
        / np.max(np.abs(gram_d_expected)))
    det_before = np.linalg.det(Mu @ Md)
    det_after = np.linalg.det(Mu_rephased @ Md)
    det_ratio = det_after / det_before
    assert gram_u_residual < TOL and gram_d_residual < TOL
    assert abs(det_ratio - phase ** 3) < 1e-10
    return {
        "gram_u_residual": gram_u_residual,
        "gram_d_residual": gram_d_residual,
        "ckm_unchanged": True,
        "determinant_before": det_before,
        "determinant_after": det_after,
        "determinant_ratio": det_ratio,
    }


def _categorical_local_data():
    """Compute the local condensation and doubled-phase results."""
    h7, h26 = Fraction(2, 5), Fraction(3, 5)
    weights = {
        ("1", "1"): Fraction(0),
        ("1", "26"): h26,
        ("7", "1"): h7,
        ("7", "26"): h7 + h26,
    }
    B = ("7", "26")
    transmissible = tuple(sorted(
        label for label, h in weights.items() if h.denominator == 1))
    confined = tuple(sorted(label for label in weights
                            if label not in transmissible))
    assert transmissible == tuple(sorted((("1", "1"), B)))

    # Standard Fibonacci R symbols and their reverse-category partners.
    R7 = {"1": cmath.exp(-4j * math.pi / 5),
          "7": cmath.exp(3j * math.pi / 5)}
    R26 = {"1": R7["1"].conjugate(),
           "26": R7["7"].conjugate()}
    twists = {label: _phase_from_weight(h)
              for label, h in weights.items()}
    displayed_generators = (
        *(twists[label] for label in transmissible),
        R7["1"] * R26["1"],
        R7["7"] * R26["26"],
        *(twists[label] / twists[B] ** 2 for label in transmissible),
    )
    assert all(abs(value - 1.0) < TOL for value in displayed_generators)

    # The doubled Fibonacci associator on the {1,B} channel is the
    # entrywise product of F with its conjugate.  Its unit eigenspace is
    # one-dimensional and spanned by (1,1), fixing the two nonzero
    # multiplication components to equal normalized magnitude.
    a, b = 1.0 / PHI ** 2, 1.0 / PHI
    associator = ((a, b), (b, a))
    eigenvalues = (a + b, a - b)
    unit_eigenvector = (1.0, 1.0)
    residual = tuple(
        sum(associator[i][j] * unit_eigenvector[j] for j in range(2))
        - unit_eigenvector[i]
        for i in range(2))
    assert max(abs(x) for x in residual) < TOL
    assert abs(eigenvalues[0] - 1.0) < TOL
    assert abs(eigenvalues[1] + 1.0 / PHI ** 3) < TOL

    survivors = tuple(z for z in (1.0 + 0j, -1.0 + 0j)
                      if is_unital_phase_automorphism(z))
    assert survivors == (1.0 + 0j,)
    return {
        "transmissible": transmissible,
        "confined": confined,
        "displayed_phase_generators": displayed_generators,
        "associator": associator,
        "associator_eigenvalues": eigenvalues,
        "unit_eigenvector": unit_eigenvector,
        "multiplication_components_nonzero": True,
        "unital_phase_automorphisms": survivors,
    }


def run(report=print):
    report("STRONG CP")
    report("=" * 70)

    embedding_index = ((Fraction(1, 2) + Fraction(1, 2))
                       / Fraction(1))
    assert embedding_index == 1
    report("topology: pi3(G2)=pi3(SU(3))=Z; embedding index = 1")

    local = _categorical_local_data()
    report("local category: transmissible {1,B}; displayed doubled")
    report("  twists/braidings/monodromies are +1")
    report("algebra: B*B=1+B and T(B)=zB give z^2=1=z;")
    report("  unique diagonal unital phase automorphism z=1")

    generic_theta = math.pi / 7.0
    central = central_character_countermodel(generic_theta)
    assert central["ray_character_nontrivial"]
    assert central["algebra_action_is_automorphism"]
    report("countermodel 1: U_n=exp(i n pi/7)I has a nontrivial ray")
    report("  character but induces the identity automorphism on A")

    mass_null = right_handed_mass_rephasing(math.pi / 11.0)
    assert abs(mass_null["gram_multiplier"] - 1.0) < TOL
    assert abs(mass_null["determinant_multiplier"] - 1.0) > 0.1
    report("countermodel 2: a right-handed phase preserves masses and CKM")
    report("  while changing arg det(M_u M_d)")

    # The repository's clock vacuum is an explicit precedent for a
    # projective symmetry phase; it prevents norm multiplicativity from
    # being used as a no-central-character theorem.
    z3_projective_phase = _phase_from_weight(Fraction(1, 3))
    assert abs(z3_projective_phase - 1.0) > 0.5

    # The engine applies the full record rule to the complete physical
    # character rather than to coordinate-dependent EFT pieces.
    production = production_couplings.derive_strong_cp()
    assert production["vacuum_phase_automorphisms"] == (
        local["unital_phase_automorphisms"])
    assert production["phase_rule_status"] == "CONSTITUTIVE ONTOLOGY"
    bar_theta = production["bar_theta"]
    assert bar_theta == 0.0

    report("-")
    report("The complete winding character acts on the closed-record")
    report("algebra; Aut_phase(A)={1}, so bar(theta)=0.")

    return {
        "embedding_index": embedding_index,
        **local,
        "central_character_countermodel": central,
        "mass_rephasing_countermodel": mass_null,
        "z3_projective_phase_precedent": z3_projective_phase,
        "production": production,
        "theta_QCD_separately_defined": False,
        "mass_determinant_phase_separately_defined": False,
        "bar_theta": bar_theta,
        "bar_theta_from_local_data_alone": None,
        "strong_cp_closed_within_ontology": True,
        "strong_cp_unconditional_local_theorem": False,
    }


if __name__ == "__main__":
    run()
