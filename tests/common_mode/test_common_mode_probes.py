"""Mechanism probes of the common mode branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The wiring and the dictionary ##



def test_bridge_decoupling(res):
    """Count and decoupling: 182 bridge modes from E8 branching. They decouple from
    the running but gravitate through the scalar heat kernel weight."""
    from probes import bridge_decoupling
    r = bridge_decoupling.run(report=_silent)
    assert r["dim_E8"] == 248
    assert r["N_bridge"] == 182
    assert r["c_coset_zero"]
    assert r["h_bridge_one"]
    assert r["bridge_share_zero"]
    assert r["bridge_a1_positive"]
    from interference import boundary
    assert boundary.a1("vector", 1.0, 0.0) == pytest.approx(-1 / 3)
    assert r["bridge_contribution_positive"]
    assert r["sigma_without_negative"]
    assert r["sigma_with_positive"]


def test_sakharov_ratio_is_a_regulator_dependent_diagnostic():
    """The Sakharov ratio 6 pi/Sigma is a regulator dependent
    diagnostic (sakharov_normalisation.py FACT 7). Under the same rows
    the vector contact coefficient +1, 0 and -1 gives 0.690, 0.416 and
    0.298, all attractive. No stiffness node reads the carrier, so no G
    is read from the ledger. The cutoff schemes sit in the ratios
    1 : 2 : 4 : 12, the q=0 dispersion is omega^2 = mu p^2 + p^4/4
    exactly, and the symmetric to clock rebase factor is computed."""
    from probes import sakharov_normalisation

    r = sakharov_normalisation.run(report=_silent)

    # 1. the ledger reproduces the engine reading
    assert r["g_ratio"] == pytest.approx(0.689617900, abs=1e-8)
    assert r["sigma"] == pytest.approx(6.0 * math.pi / 0.689617900,
                                       abs=1e-6)
    # 2. the ratio moves with the regulator and no G is read from it
    v = r["regulator_variants"]
    assert v[1.0] == pytest.approx(0.690, abs=1e-3)
    assert v[0.0] == pytest.approx(0.416, abs=1e-3)
    assert v[-1.0] == pytest.approx(0.298, abs=1e-3)
    assert all(x > 0 for x in v.values())
    assert r["gravity_reads_carrier"] is False
    assert abs(r["eh_coefficient_lpl2"] / r["eh_target_lpl2"]
               - 1.0 / r["g_ratio"]) < 1e-12
    # 3. the cutoff schemes sit in the exact ratios 1 : 2 : 4 : 12
    s = r["scheme_ratios"]
    g_r = r["g_ratio"]
    assert s["t0 = xi0^2"] == pytest.approx(g_r / 12, abs=1e-9)
    assert s["t0 = 2 xi0^2"] == pytest.approx(g_r / 6, abs=1e-9)
    assert s["t0 = l_Pl^2"] == pytest.approx(g_r / 3, abs=1e-9)
    assert s["t0 = d11 l_Pl^2"] == pytest.approx(g_r, abs=1e-9)
    # 4. substrate dispersion and vacuum rebase values
    assert r["mu_symmetric"] == pytest.approx(3.0 + math.sqrt(2.0),
                                              abs=1e-12)
    assert r["mu_clock"] == pytest.approx(3.0 - 1.0 / math.sqrt(2.0),
                                          abs=1e-12)
    assert r["k4_coefficient"] == 0.25
    assert r["clock_rebase_factor"] == pytest.approx(1.3875, abs=1e-3)


def test_one_cone_relative_doublet_is_klein_gordon_on_its_own_cone():
    """Internal speeds of a uniform clock medium: shifting the relative branches by the
    charge term 3r/4 gives Omega^2 = E0^2 + c_rel^2 p^2 + p^4/4 with
    c_rel = sqrt(5r/4) = 0.62 c_ph. The band mass sqrt(21)/5 is shared,
    the gaps split by 3r/2, and m* c_ph^2/Delta is 7.51 and 1.57."""
    from probes import one_cone

    r = one_cone.run(report=_silent)
    c = r["constants"]
    assert r["kg_error"] < 1e-12
    assert c["c_rel"] == pytest.approx(math.sqrt(1.25 * r["r"]), rel=1e-12)
    assert c["m_star"] == pytest.approx(math.sqrt(21.0) / 5, rel=1e-12)
    assert c["c_rel"] / c["c_ph"] == pytest.approx(0.620877, abs=1e-6)
    assert r["ratio_lower"] == pytest.approx(7.511626242, abs=1e-6)
    assert r["ratio_upper"] == pytest.approx(1.567767682, abs=1e-6)
    assert r["gap_split"] == pytest.approx(1.5 * r["r"], rel=1e-12)


def test_orbital_lock_has_no_invariant_term_and_no_isotropic_state():
    """Gate T at first order, the lock.  The Z3 x SO(3) invariant
    tensors are all spatially delta.  A twisted C3 lock is anisotropic
    at rank two.  A twisted T lock hides the anisotropy, but its rank-2
    content is delta.
    The rest frame counterflow tensor has rank <= 2. The g1 density is
    bounded below by -rho_tot^2/6, attained by the clock, so every
    locked state sits above it."""
    from probes import orbital_lock

    r = orbital_lock.run(report=_silent)
    sch = r["schur"]
    assert sch["Z3 x SO(3)"]["neutral_aniso"] < 1e-12
    assert sch["twisted C3"]["neutral_aniso"] > 0.5
    assert sch["twisted T"]["neutral_aniso"] < 1e-12
    assert r["boost"]["rest_min_over_max"] < 1e-12
    assert r["bound"]["min_ratio"] >= -1.0 / 6 - 1e-9
    assert r["bound"]["clock_ratio"] == pytest.approx(-1.0 / 6, abs=1e-12)
    for c in r["counterflow"]:
        assert c["excess"] == pytest.approx(c["predicted"], rel=1e-9)
        assert c["relaxed_excess"] > 0


def test_induced_horava_tensor_speed_is_on_neither_cone():
    """Gauge the Galilean symmetry and induce the Newton-Cartan
    metric's action from the clock vacuum quasiparticles.  The TT mode
    is healthy (A > 0, B > 0) with speed c_T = 2.19 c_ph = 3.53 c_rel.
    The phonon channel alone gives B/A = 5 mu, and lambda = 1/8."""
    from probes import induced_horava

    r = induced_horava.run(report=_silent)
    assert r["A"] > 0 and r["B"] > 0
    assert r["cT_over_cph"] == pytest.approx(2.1910, abs=2e-3)
    assert r["cT_over_crel"] == pytest.approx(3.5289, abs=3e-3)
    assert r["phonon_ratio"] == pytest.approx(5.0, abs=1e-6)
    assert r["lam"] == pytest.approx(0.125, abs=1e-9)
    assert r["A_spectral"] == pytest.approx(r["A"], rel=1e-6)
    assert abs(r["cT_over_cph"] - 1) > 1 and abs(r["cT_over_crel"] - 1) > 1


def test_relativistic_carrier_nr_limit_is_exact_and_one_cone():
    """Three complex Klein-Gordon fields with the U(1) x Z3 quartic have
    the Z3-NLS as their exact nonrelativistic limit. The invariant
    quartic space is 6-dimensional and a real carrier reaches only its
    4-dimensional symmetric part. The phonon has c_s^2 = c^2 k/(1+3k).
    The NR residual of the doublet scales as kappa."""
    from probes import relativistic_carrier as rc

    r = rc.run(report=_silent)
    assert r["nr_map_error"] < 1e-10
    assert r["invariants"] == 6 and r["real_image"] == 4
    k = r["per_nu"][0.01]["kappa"]
    assert r["per_nu"][0.01]["c_s_over_c"] == pytest.approx(
        math.sqrt(k / (1 + 3 * k)), rel=1e-6)
    assert r["nr_scaling"] == pytest.approx(10.0, abs=0.5)


def test_carrier_sakharov_tensor_speed_needs_scale_separation():
    """The actionless metric at the relativistic carrier level.  A
    covariant regulator gives c_T = c and lambda = 1 without a
    condensate.  A condensate pulls c_T below c as
    (m/Lambda)^2 sqrt(mu_N/m), which is order one at a Planckian healing
    length.  A spatial cutoff gives 5/7 even in vacuum (a regulator
    artefact)."""
    from probes import carrier_sakharov as cs

    r = cs.run(report=_silent)
    assert r["lambda_li"] == pytest.approx(1.0, abs=1e-6)
    for p in r["m_power"].values():
        assert p == pytest.approx(2.0, abs=0.01)
    assert r["mu_power"] == pytest.approx(0.5, abs=0.01)
    assert r["nosep"][0.1][0] < -0.3 and r["nosep"][0.01][0] < -0.3
    assert r["xi_min"] > 1e6
    assert r["spatial_ratio"] == pytest.approx(5.0 / 7.0, abs=2e-3)


def test_arena_rulers_one_clock_zero_rulers():
    """The continuous symmetry algebra of V is one dimensional (the
    common phase). The clock Hessian at fixed rho_tot has exactly one
    zero. The relative block is (5 -/+ sqrt13) g1 rho/2, paired by the
    twisted Z3 structure J (J^2 = -1, [H, J] = 0), with no zero. The
    U(1)^3 extension gives at most two rulers. Three rulers are
    incompatible with the locked circulant and with the rank bound."""
    from probes import arena_rulers

    r = arena_rulers.run(report=_silent)
    assert r["symmetry_dim"] == 1
    assert abs(r["sphere_hessian"][0]) < 1e-5
    assert min(r["relative_hessian"]) > 0.4
    g1 = 1.0 / math.sqrt(2.0)
    assert r["relative_hessian"][0] == pytest.approx((5 - math.sqrt(13)) * g1 / 2, abs=1e-5)
    assert r["relative_hessian"][-1] == pytest.approx((5 + math.sqrt(13)) * g1 / 2, abs=1e-5)
    assert r["j_error"] < 1e-12 and r["commutator"] < 1e-5
    assert r["diagonal_symmetry_dim"] == 3 and r["diagonal_relative_null"] == 2
    assert r["rulers"] == 0 and r["three_rulers_possible"] is False


def test_composite_metric_tensor_term_is_the_vacuum_term(res):
    """Written without a background metric the carrier is the theory of
    an actionless metric g = h/U, h the gradient bilinear, and its
    energy momentum tensor vanishes identically. The traceless form is
    L = 2|det h|^(1/4) - U. At one loop under a covariant regulator the
    q = 0 tensor term is -rho/4 in the full form and zero identically in
    the traceless form (I2 = 3 I1).  The q^2 term is the heat kernel
    weight a1 = 1/6 per real component.  Classical h is positive
    semidefinite: the clock gives rank one, and a texture gives spatial
    rank three."""
    from probes import composite_metric

    r = composite_metric.run(report=_silent)
    el, lp, st, led = r["elimination"], r["loop"], r["states"], r["ledger"]
    assert max(el["stress"], el["lagrangian"], el["trace_free"], el["value"]) < 1e-10
    assert r["expansion_ratio"] == pytest.approx(4.0, abs=0.05)
    assert r["expansion_ratio_tf"] == pytest.approx(4.0, abs=0.05)
    assert lp["I2"] / lp["I1"] == pytest.approx(3.0, abs=1e-6)
    assert lp["w0_full"] == pytest.approx(-lp["rho_closed"] / 4, rel=1e-6)
    assert abs(lp["w0_trace_free"]) < 1e-6 * abs(lp["w0_full"])
    assert 4 * lp["kappa"] == pytest.approx(lp["stiffness"], rel=1e-6)
    assert lp["kappa"] > 0
    assert st["clock_eigs"][2] < 1e-12 < st["clock_eigs"][3]
    assert st["random_min"] > -1e-10 and st["texture_rank3"] > 0.9
    assert st["q_identity"] < 1e-10
    assert led["str1"] == 131
    assert led["sigma"] == pytest.approx(res["g"]["Sigma_total_UV"], rel=1e-12)
    assert led["sigma_matter"] == pytest.approx(-3.0, abs=1e-12)
    assert led["lambda_G"] == pytest.approx(-0.2754, abs=1e-3)


def test_determinant_vectors_closed_form_is_maxwell_on_the_composite_metric():
    """A gauge field has a closed metric free form: the determinant of
    h + i sqrt(U) F over U.  The four dimensional determinant identity
    holds with the Pfaffian term.  The second order in F is the Maxwell
    density on the composite metric, with U cancelling.  Exact
    elimination of the auxiliary metric agrees at order F^2.  At order
    F^4 it differs by the contact term
    -(1/32)(F_h^2)^2 + (1/2) Pf^2/det h.  The density is gauge invariant
    and a scalar density.  On h = c eta the quadratic form in A is the
    flat Maxwell form for every c."""
    from probes import determinant_vectors

    r = determinant_vectors.run(report=_silent)
    assert r["identity"] < 1e-12
    assert r["maxwell"] < 1e-5 and r["u_free"] < 1e-5
    assert abs(r["order2"][2]) < 0.3 * abs(r["order2"][0])
    assert r["contact"][0] == pytest.approx(-1 / 32, abs=2e-4)
    assert r["contact"][1] == pytest.approx(0.5, abs=5e-3)
    assert r["contact_misfit"] < 1e-3 and r["residual"] < 1e-10
    assert r["gauge"] < 1e-12 and r["covariance"] < 1e-12
    assert r["flat"] < 1e-6


def test_signature_of_the_saddle_is_the_carriers():
    """The determinant form has constant saddles of every sign pattern.
    The carrier's Hamiltonian is bounded below on exactly the four
    (1,3) saddles that put the evolution coordinate on the odd sign.
    The nonhyperbolic readings have growth rates unbounded in |k|.  The
    opposite sign Lorentzian saddle is well posed, with growth bounded
    by m, but has no ground state.  In that orientation the positive
    norm Fock vacuum fills the bilinear with the trace part
    m^2 <|phi|^2>/2 times eta, plus a frame part that the covariant
    regulator does not have.  The kinetic form's cone is the cone of
    g = h/U."""
    from probes import signature_of_the_saddle

    r = signature_of_the_saddle.run(report=_silent)
    assert max(r["first_order"]) < 1e-6
    assert len(r["passing"]) == 4 and r["unbounded"] == 60
    assert max(r["rates"]["(+,-,-,-) along the +"]) == 0.0
    assert min(r["rates"]["(+,+,+,+)"][2], r["rates"]["(+,+,-,-) along a +"][2]) > 90.0
    assert all(x <= 1.0 for x in r["rates"]["(-,+,+,+) along the -"])
    assert r["rates"]["(-,+,+,+) along the -"][0] > 0.99
    for v in r["vacuum"].values():
        assert v["C"] > 0 and v["D"] > v["C"] and abs(v["virial"]) < 1e-8
    I1, contact, field = r["euclidean"]
    assert I1 == pytest.approx(contact - field, rel=1e-12)
    c2, c2b, on_cone, on_cone_b = r["cone"]
    assert c2 == pytest.approx(4.0, abs=1e-6) and c2b == pytest.approx(4.0, abs=1e-6)
    assert max(on_cone, on_cone_b) < 1e-6


def test_composite_pole_kernel_is_transverse_and_einstein():
    """Full one loop tensor kernel of the carrier's components under a
    covariant regulator (Pauli-Villars and the proper time cutoff).
    The q = 0 term is the vacuum term, and the flat metric is the
    fluctuation term at the same point.  The kernel annihilates every
    gauge direction at every momentum and has two form factors.  At low
    momentum 4 F2/q^2 is the a1 = 1/6 stiffness with F0/F2 = -2: the
    massless pole, Newton's law, gamma = 1.  With the proper time cutoff
    1/G is FACT 1's (N/6)/(2 pi t0)."""
    from probes import composite_pole

    r = composite_pole.run(report=_silent)
    for a in r["audits"].values():
        assert a["vacuum"] < 1e-9 and a["gap"] == pytest.approx(1.0, abs=1e-9)
        assert all(row["leak"] < 1e-6 and row["residual"] < 1e-6 for row in a["rows"])
        assert a["rows"][0]["stiff_ratio"] == pytest.approx(1.0, abs=1e-4)
        assert a["rows"][0]["f0_over_f2"] == pytest.approx(-2.0, abs=1e-3)
        assert a["newton"] == pytest.approx(1.0, abs=1e-4)
    assert r["I1_proper_time"] > 0 > r["rho_proper_time"]
    assert r["G_inverse"] == pytest.approx(r["fact1"], rel=1e-12)
    assert r["below"][10.0] == (True, True) and r["below"][3.0] == (True, False)
    assert abs(r["departure"][100.0]) < 2e-3
    assert r["local"] < 1e-4 and r["off_ratio"] > 1e-2
    for ratios in r["cut"].values():           # the cutoff adds no absorptive part
        assert all(abs(c - 1) < 1e-4 for c in ratios)
    assert r["below_cut"] < 1e-8
    assert r["exchange_trace_free"] == pytest.approx(1.0, abs=1e-4)


def test_ledger_kernel_rows_are_einstein_with_the_ledger_weights():
    """The second order heat trace reproduces the direct scalar kernel in
    all components (minimal and with xi R phi^2) and the charged scalar
    polarization. Each row is the linearized Einstein-Hilbert operator
    at low momentum with the ledger's weight per state: scalar 1/6 - xi
    (the amount minimally coupled, xi = 0),
    Weyl 1/12, gauge vector -1/3, Proca -1/6. The ledger's sum has
    stiffness Sigma. On the Euclidean axis its spin-2 and spin-0 form
    factors have no zero from q^2 t0 = 1e-3 to 1e3. In real time
    the kernel of a massive component is Einstein like below the cutoff
    scale and its real part changes sign between s t0 = 10 and 12."""
    from probes import ledger_kernel

    r = ledger_kernel.run(report=_silent)
    assert max(max(v) for v in r["check"].values()) < 1e-10 and r["omega_error"] < 1e-10
    want = {"scalar": 1 / 6, "Dirac": 1 / 12, "gauge vector": -1 / 3, "Proca": -1 / 6}
    for label, value in want.items():
        assert r["weights"][label][0] == pytest.approx(value, abs=2e-4)
        assert r["weights"][label][1] == pytest.approx(-2.0, abs=5e-3)
    assert r["sigma_rows"] == pytest.approx(82 / 3, abs=1e-3)
    assert r["sigma_kernel"] == pytest.approx(r["sigma_rows"], abs=2e-2)
    assert r["zeros_spin2"] == [] and r["zeros_spin0"] == []
    assert r["real_time"][10.0].real > 0 > r["real_time"][12.0].real
    assert min(r["next_order"].values()) > 19 and max(r["next_order"].values()) < 25


def test_regulator_free_ratios_higgs_is_the_dimension_two_ratio(res):
    """The induced mass term over the Einstein stiffness is 6 per real
    component under every regulator weight, while the vacuum term over
    the stiffness squared depends on the weight's shape. Under one
    regulator the Higgs mass term over 1/(16 pi G) is C_H/Sigma. With
    the engine's couplings the Standard Model rows give the Veltman
    combination -0.214 at M_Pl/sqrt(d11), zero only near 1.7e23 GeV, while
    the data need C_H/Sigma = -2.6e-33. The engine's quartic boundary
    lambda(M_Pl) = -0.01245 sits near zero with a small beta."""
    from probes import regulator_free_ratios as rf

    r = rf.run(report=_silent, res=res)
    for v in r["mass_over_stiffness"].values():
        assert v == pytest.approx(6.0, abs=1e-4)
    shape = r["vacuum_over_stiffness2"]
    assert shape["proper time t0 = 0.0001"] == pytest.approx(1.0, abs=5e-3)
    assert shape["Pauli-Villars L = 30 m"] == pytest.approx(2.107, abs=0.05)
    assert r["str1"] == 131
    assert r["transport_check"] < 1e-8
    assert r["CV_reg"] == pytest.approx(-0.214, abs=2e-3)
    assert r["CV_Pl"] == pytest.approx(-0.202, abs=2e-3)
    assert 1e23 < r["mu_zero"] < 4e23
    assert abs(r["data_ratio"]) < 1e-32 and r["induced_ratio"] < -1e-2
    assert r["yt_needed"] < r["yt_reg"] - 0.04
    assert r["bridge_need"] == pytest.approx(0.428, abs=5e-3)
    assert r["lambda_Pl"] == pytest.approx(-0.01245, abs=1e-5)
    assert abs(r["beta_lambda_Pl"]) < 1e-3
    portal = r["bridge_portal"]
    assert min(portal.values()) > 8 * r["bridge_need"]
    assert portal["m_e"] == pytest.approx(3.76, abs=0.05)
    assert portal["the regulator scale"] == pytest.approx(80.8, abs=0.5)
    assert r["M1_needed"] > 5e3
    assert r["v_over_MPl"] == pytest.approx(r["transmuted"], rel=1e-12)
    assert r["alpha_s_MZ_chain"] == pytest.approx(0.11869, abs=2e-5)
    assert r["kappa_g_tracking"] == pytest.approx(-9 / 2, abs=1e-6)
    assert r["kappa_g_both"] > 0 > r["kappa_g_top"]


def test_log_class_stiffness_holds_and_extra_modes_sit_at_the_cutoff(res):
    """The logarithmic class of the induced action is regulator free:
    over the ledger's rows c = 91/20 and d = 21/8. Switching every
    Standard Model mass off moves G by 1.2e-33. The scalaron lies
    within a factor 3.2 of the reduced Planck mass, 5.9e8 or more short
    of Starobinsky inflation, and m0^2/|m2^2| = c/(3d)."""
    from probes import log_class

    r = log_class.run(report=_silent)
    assert r["c"] == pytest.approx(91 / 20, abs=1e-12)
    assert r["d"] == pytest.approx(2.625, abs=1e-12)
    assert 1e-33 < r["accumulation"] < 2e-33
    assert r["m0_over_M"][100.0] == pytest.approx(0.317, abs=2e-3)
    assert r["starobinsky_shortfall"] > 5.9e8
    assert r["m0_over_m2_sq"] == pytest.approx(26 / 45, abs=1e-4)


@pytest.mark.slow
def test_higgs_boundary_vent_menu(res):
    """The bridge boundary, run with the Standard Model equations and
    matched to the pole mass, gives M_h = 125.08 GeV.  The vent menu
    holds two readings within 1 sigma, (1 - h10) and 1/4 (125.24 GeV).
    So the elementary route takes M_h as data, and the engine, which reads
    the boson as the condensate's radial mode, excludes it
    (regulator_free_ratios.py R6)."""
    from probes import higgs_mass

    r = higgs_mass.run(report=_silent)
    assert r["lam_mt"] == pytest.approx(0.125940, abs=2e-6)
    assert r["tree"] == pytest.approx(123.57, abs=0.02)
    assert r["M_h"] == pytest.approx(125.08, abs=0.02)
    assert abs(r["pull"]) < 1.0
    assert r["scan"][0.0] == pytest.approx(129.13, abs=0.05)
    assert r["M_h_msbar"] == pytest.approx(125.21, abs=0.03)
    assert r["landings"] == ["1/4", "h10"]
    assert r["vent_mh"]["1/4"] == pytest.approx(125.24, abs=0.03)


def test_self_consistent_graviton_row(res):
    """The on shell graviton row is -8/3, an additive row on the
    engine's full ledger.  The Coleman-Weinberg sextic is ruled out
    (200x too small).  The Sakharov ratio remains a diagnostic."""
    from probes import self_consistent_gravity

    r = self_consistent_gravity.run(report=_silent)
    assert r["sigma_ledger"] == pytest.approx(82 / 3, abs=1e-3)
    assert r["sigma_ledger"] == pytest.approx(res["g"]["Sigma_total_UV"], rel=1e-12)
    assert r["literature"]["on_shell"] == pytest.approx(-8 / 3, abs=1e-9)
    assert r["kappa_CW"] < 0.01


def test_closed_event_scale_record_closure_needs_a_planckian_carrier():
    """The smallest length the lump interior resolves is at least half
    its carrier's Compton length (c_s^2/c^2 = 2|mu0|/mc^2 <= 1/2 x 2).
    So the record closure xi = 1.930 l_Pl needs m >= 0.26 M_Pl, and the
    identification xi0 = l_Pl/2 is m = M_Pl at the stability edge.  A
    smallest closed event at the electron mass has a keV carrier whose
    interior records 8e-52 of 1/(4G).  The matter and regulator
    readings of the lump thus differ by 3e24 in carrier mass, and no
    gate tests it.  At depth 3 the minimal lump spans 654 of its own
    Compton lengths at any carrier mass, and r_s/r_half =
    689 (m/M_Pl)^2: no carrier mass makes it a particle."""
    from probes import closed_event_scale

    r = closed_event_scale.run(report=_silent)
    assert r["cs2_identity_error"] < 1e-12
    assert r["t0"] == pytest.approx(4.350, abs=2e-3)
    assert r["xi_over_lambda_c"][0.5] == pytest.approx(0.5, abs=1e-12)
    ps = r["readings"]["per state record"]
    assert ps["m_min"] == pytest.approx(0.2590, abs=5e-3)
    assert ps["E2"] == pytest.approx(0.473, abs=5e-3)
    assert r["readings"]["one record per quantum"]["m_min"] == pytest.approx(0.110, abs=2e-3)
    ident = r["readings"]["identification l_Pl/d10"]
    assert ident["m_min"] == pytest.approx(1.0, abs=1e-12)
    assert ident["E2"] == pytest.approx(2.12e-3, rel=0.02)
    assert r["per_state_m_min_o4"] == pytest.approx(0.366, abs=5e-3)
    e = r["electron"]
    assert e["binding_fraction"] == pytest.approx(0.964, abs=2e-3)
    assert e["m_carrier_MeV"] == pytest.approx(1.08e-3, rel=0.01)
    assert 5e24 < e["xi_int_lpl"] < 6e24
    assert e["record_over_quarter"] < 1e-50
    assert r["dichotomy"] > 1e24
    assert r["own_compton"] == pytest.approx(654, abs=2)
    assert r["electron_bound"] == pytest.approx(2.85e-7, rel=0.01)
    assert r["electron_lump_m"] == pytest.approx(2.53e-10, rel=0.01)
    assert r["horizon_coefficient"] == pytest.approx(689, abs=2)
    assert r["rs_over_rhalf"]["per state record"] == pytest.approx(46.8, abs=1)
    assert r["m_weak_field"] == pytest.approx(0.038, abs=1e-3)
