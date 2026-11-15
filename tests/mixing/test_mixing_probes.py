"""Mechanism probes of the mixing branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The wiring and the dictionary ##



def test_orientation_bit(res):
    """One orientation datum: the CP-odd sign flips coherently with
    the Fano orientation."""
    from probes import orientation_bit
    r = orientation_bit.run(report=_silent)
    assert r["flip_coherent"]
    assert set(r["cp_odd_forward"]) == {"eta_bar", "J_CKM", "J_lep", "eta_B"}
    for name, value in r["cp_odd_forward"].items():
        assert r["cp_odd_reversed"][name] == pytest.approx(-value)
    assert r["cp_odd_forward"]["eta_bar"] == pytest.approx(res["mx"]["eta"])
    assert r["cp_odd_forward"]["J_CKM"] == pytest.approx(res["mx"]["J_ckm"])
    assert r["cp_odd_forward"]["J_lep"] == pytest.approx(res["mx"]["J_lep"])
    assert r["cp_odd_forward"]["eta_B"] == pytest.approx(res["g"]["eta_B"])
    assert r["ckm_magnitudes_invariant"]
    assert r["pmns_angles_invariant"]


def test_sign_exhaustion_theorem(res):
    """D1: the PMNS sign rows against the data.  Under the one
    orientation bit axiom, the NuFIT 3 sigma range and eta_B > 0 leave
    one row.  It is the same row the blind gates of pmns_sign_bits.py
    leave.  U_TBM is real, so epsilon enters the lepton sector only
    through delta."""
    from probes import sign_exhaustion_theorem
    r = sign_exhaustion_theorem.run(report=_silent)
    assert r["u_tbm_real"]
    assert r["j_lep_flips"]
    assert r["angles_invariant"]
    assert r["closing_count"] == 1
    assert r["closing_row"]["delta_sign"] == -1
    assert r["closing_row"]["eta_sign"] == -1
    assert r["eta_B_positive"]


def test_strong_cp_closure_and_null_directions(res):
    """Freeze the full record closure and the countermodels
    that return when its constitutive phase rule is removed."""
    from interference import couplings as production_couplings
    from probes import strong_cp
    r = strong_cp.run(report=_silent)
    assert r["embedding_index"] == 1
    assert all(abs(g - 1.0) < 1e-12
               for g in r["displayed_phase_generators"])
    assert r["associator_eigenvalues"][0] == pytest.approx(1.0)
    assert abs(r["associator_eigenvalues"][1]) < 1.0
    assert r["unit_eigenvector"] == (1.0, 1.0)
    assert r["multiplication_components_nonzero"]
    assert r["unital_phase_automorphisms"] == (1.0 + 0j,)
    assert r["central_character_countermodel"][
        "ray_character_nontrivial"]
    assert r["central_character_countermodel"][
        "algebra_action_is_automorphism"]
    assert r["mass_rephasing_countermodel"]["masses_unchanged"]
    assert r["mass_rephasing_countermodel"]["left_ckm_unchanged"]
    assert abs(r["mass_rephasing_countermodel"][
        "determinant_multiplier"] - 1.0) > 0.1
    assert r["strong_cp_closed_within_ontology"]
    assert not r["strong_cp_unconditional_local_theorem"]
    assert r["bar_theta_from_local_data_alone"] is None
    assert not r["theta_QCD_separately_defined"]
    assert not r["mass_determinant_phase_separately_defined"]
    assert r["bar_theta"] == 0.0
    assert r["production"] == res["c"]["strong_cp"]
    assert res["c"]["bar_theta"] == 0.0
    assert production_couplings.projective_phase_status(
        has_relative_carrier=True) == "CARRIED"
    assert production_couplings.projective_phase_status(
        has_relative_carrier=False) == "QUOTIENT-IDENTICAL"

    # Instantiate the mass null direction on this repository's actual
    # outputs, not merely on an abstract Dirac matrix.
    m = res["m"]
    actual = strong_cp.engine_mass_matrix_countermodel(
        (m["m_u"], m["m_c"], m["m_t"]),
        (m["m_d"], m["m_s"], m["m_b"]),
        res["mx"]["V"], math.pi / 11.0)
    assert actual["gram_u_residual"] < 1e-12
    assert actual["gram_d_residual"] < 1e-8
    assert actual["ckm_unchanged"]
    assert abs(actual["determinant_ratio"] - 1.0) > 0.1


def test_mixing_open_path_carriers():
    """D6* has 3D carriers. Sector assignment is not web computable."""
    from probes import mixing_operator_search
    r = mixing_operator_search.run(report=_silent)
    assert set(r["d6_commutators"].values()) == {0.0}
    assert r["d6_connection_magnitude"] == pytest.approx(
        1/math.sqrt(2), abs=1e-15)
    assert r["d6_braid_residual"] < 1e-12
    assert r["d6_braid_commutator"] > 1e-12
    assert sorted(r["depth4_carrier_dimensions"]) == [2, 3, 3]
    assert r["j3_branching_projector_ranks"] == (1, 2, 8, 8, 8)
    assert r["p26_branching_projector_ranks"] == (2, 8, 8, 8)
    assert r["albert_junction"] == 1
    assert r["albert_trace_response"] == Fraction(27, 2)
    assert r["spectral_nonuniqueness"]["relative_frames_differ"]
    assert r["selected_open_transport"] is None


def test_mixing_depth4():
    """Depth-4 frames are a real finite menu. Depth 4 moves no mixing."""
    from probes import mixing_depth4
    r = mixing_depth4.run(report=_silent)
    a, b, c = r["A"], r["B"], r["C"]
    assert (a["blocks"], a["carriers"]) == (12, 24)
    assert a["root_dims"] == (3, 3, 2)
    assert a["menu"] == (Fraction(0), Fraction(1, 12), Fraction(1, 6),
                         Fraction(1, 4), Fraction(1, 3), Fraction(1, 2),
                         Fraction(2, 3), Fraction(3, 4), Fraction(1))
    assert a["nonzero_floor"] == Fraction(1, 12)
    assert a["snap_residual"] < 1e-12 and a["j_max"] < 1e-12
    assert set(b["deep_nodes"]) <= {"inv_alpha", "S_quark",
                                    "m_e", "m_mu", "m_tau", "common_mode"}
    assert b["v_shift"] < 1e-12
    assert 0 < b["ratio_shift"] < 1e-4
    assert max(abs(v) for v in c["pmns_dpull"].values()) < 1e-3
    assert c["lambda_dpull_menu"] < 0.3
    assert abs(c["tan_h10_need"]) > 5 * b["ceiling_menu"]
    assert c["ckm_best_deviation"] > c["ckm_identity_deviation"]


def test_mixing_z3_dichotomy():
    """The closure theorem gives V=1 or exact democracy."""
    import contextlib
    import io

    import masses
    import root as _root
    from probes import mixing_z3_dichotomy

    with contextlib.redirect_stdout(io.StringIO()):
        R = _root.derive()
        M = masses.derive(R)
    reading = {k: M[f"m_{k}"] for k in
               ("e", "mu", "tau", "u", "d", "s", "c", "b", "t")}
    reading["s_raw"] = M["m_s"] / M["bridge"]

    r = mixing_z3_dichotomy.run(report=_silent, engine_masses=reading)

    # 1. two Z3-equivariant sector transports commute: V = 1, exactly
    assert r["circulant_offdiagonal"] < 1e-9
    assert r["sector_commutator"] < 1e-9
    # 2. the Delta(27) escape is exact, not a leading order effect
    assert r["unbiased_deviation"] < 1e-12
    assert r["democratic_jarlskog"] == pytest.approx(
        1.0 / (6.0 * math.sqrt(3.0)), abs=1e-12)
    # 3. the derived amplitude sqrt(d10) matches the crossing triplets
    amps = r["triplet_amplitudes"]
    for name in ("leptons (e,mu,tau)", "crossing (c,b,t)",
                 "crossing (s,c,b)"):
        assert amps[name]["deviation"] < 0.01
    for name in ("up-type (u,c,t)", "down-type (d,s,b)"):
        assert amps[name]["deviation"] > 0.05
    # 4. the down channel is the simple current fixed point: the modular
    #    amplitude to the up channel is zero, and the whole menu is 0/1/2
    assert r["modular_menu"]["fixed_point_charged_overlap"] < 1e-12
    assert r["modular_menu"]["up_to_down"] < 1e-12
    assert set(round(v) for v in r["modular_menu"]["menu"]) == {0, 1, 2}
    assert r["modular_menu"]["up_to_antifundamental"] == pytest.approx(
        1.0, abs=1e-12)
    # 5. eta_bar = pi/d11^2 is ONE datum seen three ways.  Through
    #    T_ll = exp(2 pi i (h_l - c/24)), the R-matrix eigenvalue,
    #    arg S_ff and the connection identity are all functions of the
    #    same conformal weights.  Their agreement is a consistency check
    #    of the modular data, NOT three independent confirmations.
    phase = r["braiding_phase"]
    for key in ("r_matrix_symmetric", "arg_S_ff", "connection_identity"):
        assert phase[key] == pytest.approx(math.pi / 9.0, abs=1e-12)
    # 6. the CKM coordinate map has a four dimensional parameter rank.
    #    That does not make the correlated scorecard independent data.
    assert r["ckm_parameter_rank"] == 4


def test_pmns_form_sensitivity_is_scoped():
    """Freeze the declared local scan and the diagnostics that prevent
    it from being mislabeled as a complete PMNS form space grammar."""
    from probes import pmns_form_sensitivity
    r = pmns_form_sensitivity.run(report=_silent)
    assert r["legacy_n_rows"] == 144
    assert r["legacy_unique_fingerprints"] == 144
    assert r["legacy_n_within_1sigma"] == 2
    assert r["legacy_n_within_half_sigma"] == 1
    assert r["distinct_plane_n_rows"] == 96
    assert r["distinct_plane_hits"] == {1.0: 2, 0.5: 1}
    assert r["extended_n_rows"] == 288
    assert r["extended_hits"] == {1.0: 3, 0.5: 2}
    forks = r["fork_comparison"]
    assert forks["h10"]["max_abs_pull"] < 0.1
    assert 0.9 < forks["lambda"]["max_abs_pull"] < 1.0
    assert forks["atan_lambda"]["max_abs_pull"] < 0.3
    assert forks["h10"]["diagonal_score"] < \
        forks["atan_lambda"]["diagonal_score"] < \
        forks["lambda"]["diagonal_score"]
    assert r["theta_sign_promotions"] == []
    assert r["permutation_promotions"] == []
    assert r["delta_conjugation_invariant"]
    assert 0.15 < r["solar_cosine"] < 0.16
    lo, hi = r["principal_delta_band_deg"]
    assert 100.0 < lo < r["delta_deg"] < hi < 107.0


def test_pmns_sign_is_fixed_by_two_blind_gates(res):
    """pmns_sign_bits.py: the invariant PMNS phase depends only on the
    product of two bits.  One bit is the lepton pairing.  The other is
    the quark rotation that carries the CP phase.  Two gates that use no
    data rule out all but one of the four rows.  P3: only the Z_C
    pairing builds three generation carriers with a democratic block.
    P4: the charge dictionary puts the phase of S[(1,0),(1,0)] on the
    up-type rotation.  The NuFIT range picks the same row
    (sign_exhaustion_theorem.py), so the data check the sign and do not
    choose it."""
    from probes import pmns_sign_bits as pb
    from probes import sign_exhaustion_theorem as se

    out = pb.run(report=_silent)
    r = out["rows"]
    assert abs(r["U_u|Z_C"] - 275.7354) < 1e-3 and abs(r["U_d|diagonal"] - r["U_u|Z_C"]) < 1e-9
    assert abs(r["U_d|Z_C"] - 84.2646) < 1e-3
    g = pb.gates()
    assert g["blocks"] == {"Z_C": (3, True), "diagonal": (0, False)}
    assert g["rotation"] == "U_u" and out["survivors"] == ["U_u|Z_C"]
    assert r["U_u|Z_C"] == pytest.approx(res["mx"]["delta_invariant_deg"], abs=1e-6)
    data = se.run(report=_silent)
    assert data["closing_count"] == 1 and data["closing_row"]["delta_sign"] == out["delta_sign"]


def test_lambda_menu_closes_on_one_reading(res):
    """lambda_admissibility.py: the twelve entry menu lambda = h10 (1 + x/K^n)
    closes on one reading without using a measured value.  The order
    lemma, the SU(3)_3 weights, the simple current's trivial boundary
    action and f x fbar = 1 + 8 leave three entries.  The overlap
    V = U_u^dagger U_d is a Hermitian pair, which leaves the paired
    insertion h11/K^2, lambda = 73/324.
    The data check it: of the twelve entries it is the only one within 1 sigma,
    and the two single insertions the lemmas leave miss by more than 3 sigma."""
    from probes import lambda_admissibility as la
    from interference import mixing

    out = la.run(report=_silent)
    assert out["after"] == [8, 6, 4, 3, 1]
    assert sorted(out["admissible"]) == [("h10", 3, "single"), ("h11", 2, "paired"), ("h11", 3, "single")]
    assert out["survivor"] == ("h11", 2) and out["lam"] == Fraction(73, 324)
    assert res["mx"]["lam"] == pytest.approx(73 / 324, rel=1e-12)
    data = mixing.lambda_echo_menu()
    assert data["hit_keys"] == [out["survivor"]] and data["committed"] == out["survivor"]
    assert all(abs(data["menu"][(nm, m)][1]) > 3 for nm, m, kind in out["admissible"] if kind == "single")
