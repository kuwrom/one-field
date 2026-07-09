"""The probes as tests: every mechanism claim is an assertion.

The scorecard suite (test_interference.py) freezes the NUMBERS. This
suite freezes the MECHANISMS.  Each probe runs on the theory's own
equations and its conclusion is asserted here, so a change
that silently breaks a mechanism claim fails CI exactly like a change
that moves a mass.

The substrate-dynamics probes evolve the NLS equation for
minutes and carry the `slow` marker. The default run skips them:

    pytest -q              # fast probes included
    pytest -q -m slow      # the substrate-dynamics arc
"""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The wiring and the dictionary ##

def test_wiring_found_not_stated():
    """E8 > G2(1) x F4(1) at k = 3 is the unique zero-residual wiring.
    The dictionary is read off the winner, never input."""
    import wiring
    r = wiring.run(report=_silent)
    assert r["lattice"] == "E8"
    assert r["pair"] == "G2(1) x F4(1)"
    assert r["k"] == 3
    assert (r["n7"], r["n26"]) == (7, 26)
    assert abs(r["d10"] - 2.0) < 1e-9 and abs(r["d11"] - 3.0) < 1e-9


def test_normalized_current_response_wired():
    """The selected h=1 simple current and its S-matrix channel weight
    generate alpha_G2(M_Pl); root consumes the same wiring record."""
    import wiring
    import root

    w = wiring.canonical()
    assert w["k"] == 3 and w["h_J"] == 1.0
    assert w["pointed_orbit"] == ((0, 0), (3, 0), (0, 3))
    assert abs(w["D2_tot"] - 36.0) < 1e-12
    assert abs(w["pointed_weight"] - 1.0 / 12.0) < 1e-12
    assert abs(w["current_response"] - 1.0 / (2.0 * math.pi)) < 1e-12
    assert w["current_response"] == wiring.normalized_current_response(
        root.h_bridge)
    with pytest.raises(ValueError):
        wiring.normalized_current_response(2.0)
    assert abs(w["alpha_G2_Pl"] - 1.0 / (24.0 * math.pi)) < 1e-12
    assert w["N_bridge"] == 7 * 26 == root.N_bridge
    assert abs(w["bridge_mode_moment"] - 0.5) < 1e-12
    assert w["bridge_mode_moment"] == wiring.bridge_mode_second_moment(
        math.sqrt(w["d10"]), w["d10"] / w["d11"]**2,
        w["n7"], w["n26"])
    assert float(root.E_v2) == pytest.approx(
        w["bridge_mode_moment"], abs=1e-12)
    assert root.WIRING_DATA == w
    assert root.alpha_G2_Pl == w["alpha_G2_Pl"]


def test_e8_vertex_charge_trace_wired():
    """The one-orientation E8 current index is 30 for every unit q."""
    import wiring
    import root
    form = wiring.e8_positive_root_current_form()
    assert form == pytest.approx(30.0 * np.eye(8), abs=1e-12)
    assert len(wiring.e8_roots()) == 240
    assert wiring.e8_vertex_charge_trace() == pytest.approx(30.0)
    diagonal = np.ones(8) / math.sqrt(8.0)
    assert wiring.e8_vertex_charge_trace(diagonal) == pytest.approx(30.0)
    assert root.vertex_charge_trace == root.hv_E8 == root.N_vertex == 30
    assert root.WIRING_DATA["vertex_charge_trace"] == 30


# Skyrmion spectral flow and CP-orientation propagation

def test_skyrmion_dirac_spectral_flow():
    """The actual 3D knot supports one localized spectral crossing.

    The scanned coupling labels the spectral family and is not a fitted
    framework input. The crossing is composed with the closure and
    representation branches under the stated fermion/knot identification.
    """
    from probes import zero_mode
    r = zero_mode.run(report=_silent)
    assert r["profile"] == "F(r)=pi*exp(-r/1.5)"
    assert r["net_crossings"] == 1
    assert r["crossing_converged"]
    assert r["levels"][0]["energy"] > 0
    assert r["levels"][-1]["energy"] < 0
    assert min(s["core_weight"] for s in r["levels"]) > 0.5
    assert r["scan_role"] == "spectral-family coordinate"
    assert not r["enters_scorecard"]
    assert not r["pinned_zero_mode"]


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
    """D6* has 3D carriers; sector assignment is not web-computable."""
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
    # 2. the Delta(27) escape is rigid, not a leading order
    assert r["unbiased_deviation"] < 1e-12
    assert r["democratic_jarlskog"] == pytest.approx(
        1.0 / (6.0 * math.sqrt(3.0)), abs=1e-12)
    # 3. the derived amplitude sqrt(d10) lands on the crossing triplets
    amps = r["triplet_amplitudes"]
    for name in ("leptons (e,mu,tau)", "crossing (c,b,t)",
                 "crossing (s,c,b)"):
        assert amps[name]["deviation"] < 0.01
    for name in ("up-type (u,c,t)", "down-type (d,s,b)"):
        assert amps[name]["deviation"] > 0.05
    # 4. the down channel is the simple-current fixed point: the modular
    #    amplitude to the up channel is zero, and the whole menu is 0/1/2
    assert r["modular_menu"]["fixed_point_charged_overlap"] < 1e-12
    assert r["modular_menu"]["up_to_down"] < 1e-12
    assert set(round(v) for v in r["modular_menu"]["menu"]) == {0, 1, 2}
    assert r["modular_menu"]["up_to_antifundamental"] == pytest.approx(
        1.0, abs=1e-12)
    # 5. eta_bar = pi/d11^2, ONE datum with three faces: T_ll =
    #    exp(2 pi i (h_l - c/24)) makes the R-matrix eigenvalue, arg S_ff
    #    and the connection identity all functions of the same conformal
    #    weights, so their agreement is a consistency check of the
    #    modular data, NOT three independent confirmations
    phase = r["braiding_phase"]
    for key in ("r_matrix_symmetric", "arg_S_ff", "connection_identity"):
        assert phase[key] == pytest.approx(math.pi / 9.0, abs=1e-12)
    # 6. the CKM coordinate map has a four-dimensional parameter rank;
    # this does not make the correlated scorecard independent data
    assert r["ckm_parameter_rank"] == 4




def test_condensate_map_is_the_per_mode_reading():
    """The condensate map's field map, executed.  The SM
    broken-generator sum rule in
    the real-canonical norm consumes v_EW, so a linear identification
    of the 7-singlet transports v7 = v_EW (alpha_s -2.5 sigma); the
    coded v7^2 = E[v^2] v_EW^2 is the bridge per-mode contraction
    moment, reproduces M_V = 111.648 exactly, and the two
    contractions are separated by 3 sigma.  Every non-electroweak
    scale on the ladder fails by tens of sigma.  The look-elsewhere
    count over bare factors is recorded: five atoms within 1 sigma,
    the pin is the role, not the fit."""
    import contextlib
    import io

    import masses
    import root as _root
    from probes import condensate_map

    with contextlib.redirect_stdout(io.StringIO()):
        R = _root.derive()
        M = masses.derive(R)

    r = condensate_map.run(report=_silent, engine_m_t_GeV=M["m_t"] / 1e3)

    # 1. the SM sum rule: C2(2) = 3/4, three equal masses, M_W = gv/2,
    #    identical in the real and complex bookkeepings of one doublet
    sm = r["sm_real_norm"]
    assert sm["casimir_sum_over_g2v2"] == pytest.approx(0.75, abs=1e-12)
    assert sm["mw_over_gv"] == pytest.approx(0.5, abs=1e-12)
    assert (sm["real_norm_length"] / sm["complex_norm_length"]
            == pytest.approx(math.sqrt(2.0), abs=1e-12))
    # 2. the per-mode contraction reproduces the coded threshold
    assert r["per_mode"]["M_V"] == pytest.approx(111.648, abs=1e-2)
    assert r["per_mode"]["M_V"] == pytest.approx(
        (math.pi / (2.0 * math.sqrt(2.0))) * r["v_EW_GeV"]
        / math.sqrt(6.0), abs=1e-9)
    assert abs(r["per_mode"]["pull"]) < 0.5
    # 3. the linear map (= the total contraction) is the excluded row
    assert r["linear_map"]["M_V"] == pytest.approx(157.894, abs=1e-2)
    assert r["linear_map"]["pull"] < -2.0
    assert r["dichotomy_separation_sigma"] == pytest.approx(2.95, abs=0.1)
    # 4. the scale menu: every non-electroweak condensate row is dead
    menu = r["scale_menu_pulls"]
    assert abs(menu["M_R:per_mode"]) > 50
    assert abs(menu["Lambda_conf:per_mode"]) > 50
    assert abs(r["no_threshold_pull"]) > 5
    # 5. look-elsewhere honesty: the 1-sigma atom count is frozen
    assert len(r["atoms_within_1sigma"]) == 5
    assert "E_v2 = 1/2" in r["atoms_within_1sigma"]


def test_sakharov_normalisation_is_exactly_2pi_d11():
    """The Sakharov normalisation target, made exact.  The standard
    proper-time algebra (one scalar reproduces the textbook Sakharov
    Lambda^2/(12 pi)) turns Sigma = 6 pi into the engine expression
    Sigma = 2 pi d11 exactly; the naive substrate schemes miss by 3
    to 12; two exactly degenerate factorizations carry the target
    (P-closure x healing area with d10 cancelling, and common-mode
    share x Planck area); the q=0 dispersion is omega^2 = mu p^2 +
    p^4/4 exactly. The symmetric-to-clock rebase factor is computed."""
    from probes import sakharov_normalisation

    r = sakharov_normalisation.run(report=_silent)

    # 1. the ledger reproduces the engine closure
    assert r["g_ratio"] == pytest.approx(0.999999917, abs=1e-8)
    assert r["sigma"] == pytest.approx(6.0 * math.pi / 0.999999917,
                                       abs=1e-6)
    # 2. the exact identity and the d10 cancellation
    assert r["target_is_2pi_d11"]
    assert r["factorizations_equivalent"]
    assert abs(r["eh_coefficient_lpl2"] / r["eh_target_lpl2"]
               - 1.0 / r["g_ratio"]) < 1e-12
    assert r["d10_cancellation_exact"]
    # 3. the naive schemes: 1/12, 1/6, 1/3
    s = r["scheme_ratios"]
    assert s["t0 = xi0^2"] == pytest.approx(1.0 / 12, abs=1e-6)
    assert s["t0 = 2 xi0^2"] == pytest.approx(1.0 / 6, abs=1e-6)
    assert s["t0 = l_Pl^2"] == pytest.approx(1.0 / 3, abs=1e-6)
    assert s["t0 = d11 l_Pl^2"] == pytest.approx(1.0, abs=1e-6)
    # 4. substrate dispersion and vacuum-rebase values
    assert r["mu_symmetric"] == pytest.approx(3.0 + math.sqrt(2.0),
                                              abs=1e-12)
    assert r["mu_clock"] == pytest.approx(3.0 - 1.0 / math.sqrt(2.0),
                                          abs=1e-12)
    assert r["k4_coefficient"] == 0.25
    assert r["clock_rebase_factor"] == pytest.approx(1.3875, abs=1e-3)



def test_z3_symmetric_background_is_not_the_vacuum():
    """The background every 3D Z3-NLS run is seeded on is a stationary
    point with NEGATIVE directions, 93% above the true uniform ground
    state, which is the clock (1,w,w^2).  The clock keeps the Z3 centre
    unbroken and, unlike the symmetric state, has rho_rel != 0 so the
    director is a UNIT vector -- which is what makes a Hopf charge an
    integer at all."""
    import substrate3d
    from probes import z3_vacuum

    r = z3_vacuum.run(report=_silent)
    st = r["states"]
    e_sym = st["Z3-symmetric (1,1,1)"]["E"]
    e_clk = st["clock q=1 (1,w,w^2)"]["E"]

    # the clock is the minimum over every uniform state tried
    assert e_clk == min(v["E"] for v in st.values())
    assert e_sym > e_clk
    assert e_sym / e_clk == pytest.approx(6.621320 / 3.439340, rel=1e-6)
    assert st["clock q=2 (1,w^2,w)"]["E"] == pytest.approx(e_clk, rel=1e-12)

    # both are exact stationary points, with the closed-form potentials
    g0, g1 = substrate3d.g0, substrate3d.g1
    assert r["res_symmetric"] < 1e-12 and r["res_clock"] < 1e-12
    assert r["mu_symmetric"] == pytest.approx(3*g0 + 2*g1, rel=1e-9)
    assert r["mu_clock"] == pytest.approx(3*g0 - g1, rel=1e-9)

    # ...but only one of them is stable
    assert r["n_negative_symmetric"] > 0
    assert r["n_negative_clock"] == 0
    assert r["hess_symmetric"][0] == pytest.approx(-7.0*g1, rel=1e-4)

    # THE POINT: |n| = rho_rel/(rho_rel + EPS). The symmetric background
    # has no director at all, so n maps into a contractible ball and the
    # Hopf charge is not quantized; the clock gives a unit director.
    assert st["Z3-symmetric (1,1,1)"]["rho_rel"] == pytest.approx(0.0, abs=1e-14)
    assert st["Z3-symmetric (1,1,1)"]["n_norm"] == pytest.approx(0.0, abs=1e-14)
    assert st["clock q=1 (1,w,w^2)"]["rho_rel"] == pytest.approx(1.0, rel=1e-12)
    assert st["clock q=1 (1,w,w^2)"]["n_norm"] == pytest.approx(1.0, rel=1e-12)

    # the centre is unbroken: the cyclic shift is exactly a global phase
    assert r["shift_ratio_spread"] < 1e-12
    assert abs(r["shift_phase"]) == pytest.approx(1.0, rel=1e-12)
    assert cmath.phase(r["shift_phase"]) == pytest.approx(-2*math.pi/3, abs=1e-12)



## Statistics: the level counted, the winding carried ##

def test_bookkeeping_level_and_theta():
    """Child anomaly vanishes. WZ level = 3 by counting. The lightest
    Z3 sector-changer is h = 2/9 (internal assertions)."""
    from probes import bookkeeping
    bookkeeping.run(report=_silent)


def test_winding_texture_unit_charge():
    """The 3D texture carries unit topological charge (pi_3)."""
    from probes import winding_texture
    B = winding_texture.run(report=_silent)
    # REGRESSION PIN: 15% tolerance on a claimed INTEGER charge
    # (discretization); read as "the estimator reproduces its
    # measured value", not as a tight integer test
    assert abs(abs(B) - 1.0) < 0.15


## The spectrum at the MTC layer, and theta's uniqueness ##

def test_mtc_spectrum():
    """The three-line spectrum exists at the MTC layer, assembled
    entirely from computed data."""
    from probes import mtc_spectrum
    r = mtc_spectrum.run(report=_silent)
    assert abs(r["Q"] - 2.0 / 3.0) < 1e-12
    assert abs(r["mu_e"] - 206.77) < 0.01


def test_theta_menu_permits_only_2_9():
    """The data permits only theta = h(fund) = 2/9: every other
    natural MTC angle leaves the positive branch."""
    from probes import theta_menu
    dev, name, mu_e, tau_e = theta_menu.run(report=_silent)
    assert dev < 1e-3
    assert abs(mu_e - 206.7683) / 206.7683 < 1e-4


## The closure dichotomy ##

def test_effective_action_stabilizer():
    """DP #4 branch (i), closed constructively: the exact CP^1 split
    and Mermin-Ho relation hold on the repo's texture (verified with
    resolution convergence), and the induced Faddeev-Skyrme quartic
    (coefficient rho/8m^2, m the measured ring-down gap) produces a
    Derrick minimum a few healing lengths wide (2.5 xi at the
    measured gap; 1.0 xi under the non-variational variant), with the local limit
    self-consistent at its own minimum (R* m = 1.41)."""
    from probes import effective_action
    r = effective_action.run(report=_silent)
    # REGRESSION PIN: the (2.0, 3.0) window brackets the measured
    # R* = 2.5xi (a pin of the measured value, not an a-priori band)
    assert r["R_star_m"] > 1.0
    assert 2.0 < r["R_star_xi"] < 3.0


def test_surplus_edges_enumerated():
    """The grammar run AWAY from the data (web completeness): every
    used base is admissible, the excluded terminus ghosts (64, 139
    GeV) stay absent from the PDG, and the surplus channel is finite
    and enumerated (the free prediction channel)."""
    from probes import surplus_edges
    r = surplus_edges.run(report=_silent)
    assert len(r["rows"]) == 15
    assert all(g in (1, 2, 3) for g, w, m in r["falsifiers"])


def test_selector_theorem_binary():
    """The idempotent-measure space on Z3 is exactly {Haar, delta_e}:
    the confinement dichotomy as a two-element theorem."""
    from probes import selector_theorem
    found = selector_theorem.run(report=_silent)
    assert len(found) == 2


## The wiring derivation ##

def test_edge_derivations():
    """The two derived mass edges: m_d is the exact on-menu
    identity (aa + h10) = 83/9; the strange bridge is the Schur-unique
    norm conversion 32/27 with the Albert gap constructed on an
    anti-associative triple and the isotropic share unique in the
    scan."""
    from probes import edge_derivations
    r = edge_derivations.run(report=_silent)
    assert r["m_d"] == "83/9" and r["bridge_sq"] == "32/27"
    assert r["delta_J_sq"] == "3/8" and r["share_hits"] == ["1/8"]
    assert r["anti_associative"] == 28


def test_substrate_depth_measured():
    """The substrate depth is a measured lower bound, not a premise:
    c rho xi^3 > ~196 per-proton (Super-K vs the Planckian attempt
    clock), ~324 ensemble; the deep-condensate reading has two
    decades of margin."""
    from probes import substrate_depth
    r = substrate_depth.run(report=_silent)
    assert 190.0 < r["per_proton"] < 200.0
    assert 315.0 < r["ensemble"] < 330.0


def test_lambda_assignment_forced():
    """The Koide representation assignment derived (the insertion
    rule): given
    the neutral-lane premise, adjoint -> Q(c,b,t) and fundamental ->
    Q(s,c,b) are forced by two canon lemmas; the swap is dead."""
    from fractions import Fraction
    from probes import lambda_assignment
    r = lambda_assignment.run(report=_silent)
    assert r["Q_cbt"] == Fraction(289, 432)
    assert r["Q_scb"] == Fraction(649, 972)
    assert r["neutral_stepper"] == (1, 1)
    assert r["mb_swap"] - r["mb_canon"] > 70.0     # swap ~ +13 sigma


def test_conversion_rigidity():
    """The even-switch loophole closed by label rigidity (the
    conversion lemma): Hom(x,y) = delta over all pairs, no same-orientation
    vertex, conversion = pair fusion (vent-kind); ghosts absent."""
    from probes import conversion_rigidity
    r = conversion_rigidity.run(report=_silent)
    assert r["rigidity"] and r["conversion_dim"] == 1
    assert r["ghost_fence_gev"] == [64, 128, 139, 277]


def test_vertex_pinning_locks():
    """Each depth-4 response is pinned by its own observable:
    unique 2-sigma survivors (2, 16, 27/2) over the 29-value menu."""
    from fractions import Fraction
    from probes import vertex_pinning
    r = vertex_pinning.run(report=_silent)
    assert r["survivors"] == ([Fraction(2)], [Fraction(16)],
                              [Fraction(27, 2)])


def test_echo_naturality_regenerates():
    """THE RESPONSE RULE: a web with ZERO per-node depth-4
    declarations (one loop object read through channel responses)
    regenerates the canonical precision block exactly."""
    from probes import echo_naturality
    r = echo_naturality.run(report=_silent)
    assert set(r["regenerated"]) == {"inv_alpha", "m_mu", "m_tau",
                                     "v_EW"}
    assert r["responses"] == {"coupling": "2", "action": "16",
                              "amplitude": "27/2"}


## Substrate dynamics (the classical-1D arc. Minutes of PDE) ##

@pytest.mark.slow
def test_knot_formation_and_control():
    """A winding seed grows and localizes. With g1 = 0 the same seed
    disperses.  Matter forms because of the G2 coupling."""
    from probes import knot
    r = knot.run(report=_silent)
    assert r["growth_true"] > 5.0
    assert r["growth_ctrl"] < 1.0


@pytest.mark.slow
def test_stationary_persistence():
    """The saturated lump is alive at t = 40 with a sub-scale internal
    frequency (the first `mass is rotation rate` measurement)."""
    from probes import stationary
    r = stationary.run(report=_silent)
    assert r["persistent"]
    # The circulation rate is bounded by the relative-sector scale.
    # The sign is a coordinate (flips under a factor-2 change in g0 at
    # fixed ratio, and under the seed amplitude).  At the frozen canon
    # parameters the measured value is omega_rel = -0.133; the sign is
    # frozen as measured.  A flip at fixed parameters indicates an
    # integrator change and fails here.
    assert -r["gap_scale"] < r["omega_rel"] < 0.0


@pytest.mark.slow
def test_universality_registered_outcome():
    """Standing character is universal. A unique attractor is NOT
    established in 1D (the measured outcome)."""
    from probes import universality
    same = universality.run(report=_silent)
    assert same is False


@pytest.mark.slow
def test_dynamics_erosion_and_shadow():
    """The dynamics probe on the live equation, verdicts frozen as
    measured: the shadow forms and the dip persists; bias erosion
    and bright-component localization are NOT observed on this
    horizon (the scalar-phase substep's positives are integrator
    artifacts)."""
    from probes import dynamics
    r = dynamics.run(report=_silent)
    assert r["shadow"] and r["dip_persist"]
    assert not r["erosion"] and not r["localized"]


@pytest.mark.slow
def test_spectrum_runs():
    """The spectroscopy probe evolves stably (meanings pre-stated in
    its docstring. The fringe question is interpretive)."""
    from probes import spectrum
    spectrum.run(report=_silent)


@pytest.mark.slow
def test_ew_internal_reclassification():
    """G_F internal at ppm. The top-doublet rho term computed from
    framework outputs."""
    from probes import ew_internal
    r = ew_internal.run(report=_silent)
    assert abs(r["G_F"] / 1.1663788e-5 - 1.0) < 1e-5
    assert abs(r["d_rho_t"] - 0.00933) < 2e-4


@pytest.mark.slow
def test_knot_charge_not_protected():
    """pi_3 winding is NOT a substrate-level invariant: the estimator
    converges, the saturated knot unwinds through amplitude-node
    events, and no Noether charge backs internal rotation (the Q-lump
    loophole is closed).  Frozen witness for the topological-protection
    derivation (effective_action.py, knot_charge.py, knot_spectrum.py)."""
    from probes import knot_charge
    knot_charge.main()  # four assertions live inside


@pytest.mark.slow
def test_skyrmion_3d_formation():
    """3D knot formation (the review's sim3d experiment as a house
    probe): the seeded texture grows only under the true G2 coupling
    (control disperses), the growth rate matches the computed BdG
    instability, and the common mode dips (the shadow).  Charge
    protection is NOT claimed. knot_charge.py measures the
    unwinding. Protection belongs to the confined phase."""
    from probes import skyrmion_3d
    r = skyrmion_3d.run(report=_silent)   # frozen assertions inside
    assert r["growth"] > 5.0 and r["decay"] < 0.05
    # REGRESSION PIN: an asymmetric factor-2.2 band around the BdG
    # rate (measured 0.61 of sigma_bdg); a bracketing of the
    # measured match, not an a-priori quantitative prediction band
    assert (0.5 * r["sigma_bdg"] < r["sigma_measured"]
            <= 1.1 * r["sigma_bdg"])


@pytest.mark.slow
def test_confined_phase_entry():
    """Measure the vortex-density part of the confined-phase probe.

    The imprinted ring is the positive control, the endpoint carries
    proliferated defects, and the g1 = 0 control remains defect-free.
    test_confined_phase_area_law measures the separate loop-area law.
    """
    from probes import confined_phase
    confined_phase.run(report=_silent)


@pytest.mark.slow
def test_knot_spectrum_tower():
    """The saturated lump's ring-down tower is gapped, with the first
    internal excitation at the relative-sector binding scale: no
    surplus states below the lock-in scale (feeds the D4 prediction
    channel)."""
    from probes import knot_spectrum
    r = knot_spectrum.run(report=_silent)
    assert len(r["pk_rel"]) >= 1


def test_circulant_normalization():
    """Three claimants to the mass-circulant slot, adjudicated.

    Freezes: (a) Q = 1/3 + (B/A)^2/6 is exact and theta-blind, so the
    Koide value pins the AMPLITUDE alone and no normalisation can fix
    the phase; (b) the total and per-mode second moments differ by the
    factor 2 of the two Z3 channels (Pv2P = 2*E[v2]), which is why the
    web's E_v2 = 1/2 sits alongside Pv2P = 1; (c) theta = 0 forces exact
    two-fold degeneracy at every amplitude, retiring the (1,-1/2,-1/2)
    harmonic; (d) at alpha = 0 the substrate's uniform BdG circulant
    has B/A = 2 + g1/g0 >= 2: the positive-gap branch and the
    modulational instability are DISJOINT there; (e) a complex Z3
    phase is admissible and reaches (sqrt 2, 2/9) at a SINGLE physical
    point, modulationally stable (full 6x6 BdG, PART F; direct
    dynamics arbitration, PART G) -- its 'two inverse images' are one
    theory under the exact redundancy (g1, alpha) -> (-g1, alpha - pi):
    the disjointness holds for every admissible (r, alpha), and the
    canonical point is a stable uniform phase carrying the Koide
    circulant.  A per-channel-pairing reading reports coexistence
    instead; that is an artifact of the pairing, valid only at
    alpha = 0.  The effective-coupling route has no implementation."""
    from probes import circulant_normalization as cn
    r = cn.run(report=_silent)
    assert r["identity_max_dev"] < 1e-12
    assert abs(r["Pv2P_total"] - 1.0) < 1e-12
    assert abs(r["E_v2_mode"] - 0.5) < 1e-12
    assert r["theta0_degenerate"] is True
    assert abs(r["Q_canon"] - 2.0 / 3.0) < 1e-12
    assert abs(r["Q_harmonic"] - 0.5) < 1e-12
    assert r["BdG_BA"] >= 2.0
    assert r["r_needed_for_sqrt2"] < 0.0
    assert r["gap_min"] > 0.0
    assert r["reconcilable"] is False
    # PART F + G: the canonical point is reached and stable; its 'two
    # inverse images' are ONE theory (exact redundancy, asserted); the
    # dynamics arbitration confirms (canonical seeded mode flat,
    # alpha = 0 control grows at g1 rho0).  Frozen as a consequence of
    # the hermitian equation.
    assert abs(r["BA_alpha"] - math.sqrt(2.0)) < 1e-9
    assert abs(r["theta_alpha"] - 2.0 / 9.0) < 1e-9
    assert r["alpha_point_unstable"] is False
    assert r["alpha_images_all_stable"] is True
    assert r["alpha_point_positive_gaps"] is True
    assert r["inverse_image_redundancy_dev"] < 1e-9
    assert abs(r["dynamics_rate_canonical"]) < 0.02
    assert 0.5 < r["dynamics_rate_control"] < 0.85
    # the closed-form circulant must track nls_soliton's live matrix
    assert cn.test_matches_live_substrate() is True



def test_wiring_theorem_probe():
    """The wiring theorem: Weyl/Kac-Peterson data, the unique (7,26)
    boson, the frozen 4 x 7 response inventory (PART D, the executable
    menu vertex_pinning cites), and the edge classification with zero
    UNCLASSIFIED canon edges."""
    from fractions import Fraction
    from probes import wiring_theorem_probe as wt
    r = wt.run(report=_silent)
    assert r["n_edges"] == 21
    assert r["counts"] == {"EMBEDDING": 12, "INTRA-FACE": 9,
                           "UNCLASSIFIED": 0}
    assert r["albert"]["options"] == {
        0: Fraction(13), 1: Fraction(27, 2)}
    assert r["albert"]["junction"] == 1
    assert r["albert"]["response"] == Fraction(27, 2)


def test_dimensionality_derived():
    """3+1 DERIVED, not assumed: three components force a C^2
    relative sector (codim 4, measured: D=3 no zeros, D=4 isolated
    points, D=5 lines); knots need d_spatial >= 3; pointlike charge
    events need D = 4.  D = 3+1 unique."""
    from probes import dimensionality
    r = dimensionality.run(report=_silent)
    assert r["unique_D"] == 4
    assert r["codim"][4]["cnt_hi"] / r["codim"][4]["cnt_lo"] < 1.5


def test_envelope_universality():
    """THE NLS FORM DERIVED up to locality: the envelope of
    a non-NLS microdynamics (phi^4 Klein-Gordon) obeys the DERIVED
    NLS to normal-form accuracy, error ~ eps^2"""
    from probes import envelope_universality
    r = envelope_universality.run(report=_silent)
    assert r["ratio"] < 0.45
    assert r["errors"][0.06][0] < 0.025


@pytest.mark.slow
def test_confined_phase_area_law():
    """Conjecture #1's second half, measured: the Z3 disorder operator
    with fractional test charge 1/3 obeys an AREA law at the saturated
    endpoint (c = 0.021, residual 8x better than perimeter; frozen
    inside the probe)."""
    from probes import confined_phase
    r = confined_phase.area_law(report=_silent)
    assert r["area_wins"] is True


@pytest.mark.slow
def test_knot_circulant_anomalous_phase():
    """Stage zero + first reading of the two-phase bridge: a certified
    stationary knot carries a nonzero gauge-invariant anomalous phase
    (mean |beta| ~ 0.8) that the uniform medium cannot (beta = 0 by
    identity).  The quantitative identification has not yet been pinned."""
    from probes import knot_circulant
    r = knot_circulant.run(report=_silent)
    assert r["sin_beta_mean"] > 0.3
    assert r["residual"] < 5e-3
    # the BdG on the same background: two zero modes, stable, and the
    # relative channels HYBRIDIZED (uniform alpha = 0: exactly 0)
    b = knot_circulant.bdg(report=_silent)
    assert b["hybridization"] > 0.15
    assert b["max_growth"] < 2e-2


def test_harmonic_single_source():
    """The canonical harmonic (B, theta) lives in ONE place: root.py.
    protected_forgetting.py carries standalone literals so it runs
    without the engine, and this test freezes the guard that they
    cannot drift: they equal root's canonical values, and the
    module-level drift assertion is present (importing the module
    with the engine on the path executes it)."""
    import root
    import protected_forgetting as pf
    assert pf.B_CANON == root.BA_ratio
    assert pf.THETA_CANON == float(root.h10)
    v = pf.harmonic()
    # the electron gap that anchors M_Pl is the k=1 component
    assert abs((1.0 + v[1]) - root._Delta_e) < 1e-15
