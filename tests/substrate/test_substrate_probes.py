"""Mechanism probes of the substrate branch: every mechanism claim is an assertion."""

import cmath
import math
from fractions import Fraction

import numpy as np
import pytest

def _silent(*args, **kwargs):
    pass


## The wiring and the dictionary ##



# Skyrmion spectral flow and CP orientation propagation

def test_skyrmion_dirac_spectral_flow():
    """The repository's static hedgehog profile supports one localized
    spectral crossing.

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
    assert r["scan_role"] == "spectral family coordinate"
    assert not r["enters_scorecard"]
    assert not r["pinned_zero_mode"]


def test_z3_symmetric_background_is_not_the_vacuum():
    """The background every 3D Z3-NLS run is seeded on is a stationary
    point with NEGATIVE directions, 93% above the true uniform ground
    state, which is the clock (1,w,w^2).  The clock keeps the Z3 centre
    unbroken.  Unlike the symmetric state it has rho_rel != 0, so the
    director is a UNIT vector.  That is what makes a Hopf charge an
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

    # both are exact stationary points, with the closed form potentials
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
    # Hopf charge is not quantized. The clock gives a unit director.
    assert st["Z3-symmetric (1,1,1)"]["rho_rel"] == pytest.approx(0.0, abs=1e-14)
    assert st["Z3-symmetric (1,1,1)"]["n_norm"] == pytest.approx(0.0, abs=1e-14)
    assert st["clock q=1 (1,w,w^2)"]["rho_rel"] == pytest.approx(1.0, rel=1e-12)
    assert st["clock q=1 (1,w,w^2)"]["n_norm"] == pytest.approx(1.0, rel=1e-12)

    # the centre is unbroken: the cyclic shift is exactly a global phase
    assert r["shift_ratio_spread"] < 1e-12
    assert abs(r["shift_phase"]) == pytest.approx(1.0, rel=1e-12)
    assert cmath.phase(r["shift_phase"]) == pytest.approx(-2*math.pi/3, abs=1e-12)


def test_winding_texture_unit_charge():
    """The 3D texture carries unit topological charge (pi_3)."""
    from probes import winding_texture
    B = winding_texture.run(report=_silent)
    # REGRESSION PIN: 15% tolerance on a claimed INTEGER charge
    # (discretization).  It checks that the estimator reproduces its
    # measured value.  It is not a tight integer test.
    assert abs(abs(B) - 1.0) < 0.15


def test_substrate_depth_measured():
    """The substrate depth is a measured lower bound, not a premise:
    c rho xi^3 > ~196 per proton (Super-K vs the Planckian attempt
    clock), ~324 ensemble. The deep condensate reading has two
    decades of margin."""
    from probes import substrate_depth
    r = substrate_depth.run(report=_silent)
    assert 190.0 < r["per_proton"] < 200.0
    assert 315.0 < r["ensemble"] < 330.0


## Substrate dynamics (classical 1D and 3D runs, minutes of PDE time) ##

@pytest.mark.slow
def test_spectrum_runs():
    """The spectroscopy probe evolves stably, and the prestated
    outcome (a) is absent: classically no three lines sit in the
    circulant fringe ratios 1 : 14.4 : 59 (mtc_spectrum.py reads the
    fringe at the MTC layer)."""
    from probes import spectrum
    r = spectrum.run(report=_silent)
    assert r["lines"]["a1"] and r["lines"]["a2"]
    assert r["fringe_found"] is False


def test_effective_action_stabilizer():
    """Setting: the gas sign on a fixed uniform background, which the
    knot_charge, knot_in_lump and substrate_depth probes build on.  The
    exact CP^1 split and the Mermin-Ho relation hold on the repo's
    texture (checked by resolution convergence).  The induced
    Faddeev-Skyrme quartic (coefficient rho/8m^2, m the measured ring
    down gap) gives a Derrick minimum a few healing lengths wide: 2.5 xi
    at the measured gap, 1.0 xi under the nonvariational variant.  The
    local limit is self-consistent at its own minimum (R* m = 1.41)."""
    from probes import effective_action
    r = effective_action.run(report=_silent)
    # R* kappa = sqrt(I4/I2) = 1.41 holds for every conversion of the
    # measured frequency gap to a screening wavenumber. The width sits
    # between 1.9 and 5.3 healing lengths across the two conversions
    # that bracket the dispersion inside the knot, which is not measured.
    assert r["R_star_m"] > 1.0
    assert 1.5 < r["R_star_xi_lo"] < r["R_star_xi_hi"] < 6.0


@pytest.mark.slow
def test_knot_charge_not_protected():
    """pi_3 winding is NOT a substrate level invariant.  The estimator
    converges, the saturated knot unwinds through amplitude node
    events, and no Noether charge protects internal rotation (the
    Q-lump loophole is closed).  This freezes the result for the
    topological protection derivation (effective_action.py,
    knot_charge.py, knot_spectrum.py)."""
    from probes import knot_charge
    knot_charge.main()  # four assertions live inside


@pytest.mark.slow
def test_skyrmion_3d_formation():
    """3D knot formation (the review's sim3d experiment, rebuilt as a
    repo probe).  The seeded texture grows only under the true G2
    coupling (the control disperses).  The growth rate matches the
    computed BdG instability, and the common mode density dips at the
    core.  Charge protection is NOT claimed: knot_charge.py measures
    the unwinding, and protection belongs to the confined phase."""
    from probes import skyrmion_3d
    r = skyrmion_3d.run(report=_silent)   # frozen assertions inside
    assert r["growth"] > 5.0 and r["decay"] < 0.05
    # REGRESSION PIN: an asymmetric factor 2.2 band around the BdG
    # rate (measured 0.61 of sigma_bdg).  The band brackets the
    # measured match.  It is not an a priori prediction band.
    assert (0.5 * r["sigma_bdg"] < r["sigma_measured"]
            <= 1.1 * r["sigma_bdg"])


@pytest.mark.slow
def test_confined_phase_entry():
    """Measure the vortex density part of the confined phase probe.

    The imprinted ring is the positive control, the endpoint carries
    proliferated defects, and the g1 = 0 control remains defect free.
    test_confined_phase_area_law measures the separate loop area law.
    """
    from probes import confined_phase
    confined_phase.run(report=_silent)


@pytest.mark.slow
def test_knot_spectrum_tower():
    """The saturated lump's ring down tower is gapped, with the first
    internal excitation at the relative sector binding scale.  There
    are no extra states below the lock in scale.  This result feeds the
    D4 prediction channel."""
    from probes import knot_spectrum
    r = knot_spectrum.run(report=_silent)
    assert len(r["pk_rel"]) >= 1


@pytest.mark.slow
def test_confined_phase_area_law():
    """The second half of conjecture #1, measured.  The Z3 disorder
    operator with fractional test charge 1/3 obeys an AREA law at the
    saturated endpoint (c = 0.021, residual 8x better than a perimeter
    law).  The values are frozen inside the probe."""
    from probes import confined_phase
    r = confined_phase.area_law(report=_silent)
    assert r["area_wins"] is True


def test_substrate_coupling_is_the_inverse_depth_over_8pi():
    """The medium mediated force between two quanta, in units of hbar c_s,
    is 1/(8 pi rho xi^3) for any equation of state; alpha_G2 = 1/(24 pi)
    then gives rho xi^3 = 3."""
    from probes import substrate_coupling

    r = substrate_coupling.run(report=_silent)
    assert r["max_deviation"] < 1e-6
    assert r["depth"] == pytest.approx(3.0, abs=1e-12)
    assert r["depletion"] == pytest.approx(1.40e-3, abs=2e-5)


def test_lump_in_empty_field_clock_inside_weight_is_energy():
    """In the empty field the self-bound sign holds Z3 lumps with
    the clock inside: N_min = 43.87, bound for N > 55.67, the clock rate
    below m c^2, dQ/dOmega < 0 on the branch. The integrated stress
    vanishes (von Laue), so weight = energy with binding included and
    gamma = 1. The Eotvos parameter between the smallest and largest
    lump is below 1e-10 against 0.29 for weightless binding. In 1+1
    dimensions with the full three component force a boosted lump moves
    at v with E = gamma E0, P = gamma E0 v, at 0.9 c above its interior
    sound speed, 0.757 c at the probe's unit choice m c^2 = 3.49 |mu0|,
    and noise leaves it a clock."""
    from probes import lump_in_empty_field

    r = lump_in_empty_field.run(report=_silent)
    assert r["force_error"] < 1e-12
    assert r["n_min"] == pytest.approx(43.8675, abs=2e-3)
    assert r["n_bound"] == pytest.approx(55.673, abs=2e-2)
    for d, l in zip(r["lumps"], r["relativistic"]):
        assert abs(d["s"]) < 1e-10 * d["K"]
        assert abs(l["weight"] - 1) < 1e-10 and abs(l["gamma_ppn"] - 1) < 1e-10
        assert l["omega"] < 1
    assert r["relativistic"][-1]["binding"] == pytest.approx(-0.24964, abs=1e-4)
    assert r["eotvos"] < 1e-10 and r["eotvos_unbound"] == pytest.approx(0.288, abs=2e-3)
    assert r["c_sound"] == pytest.approx(0.75671, abs=1e-5)
    for v, b in r["boosts"].items():
        assert b["v_fit"] == pytest.approx(v, abs=2e-4)
        assert b["E"] == pytest.approx(b["gamma"] * r["rest_1d"]["E0"], rel=1e-6)
    assert r["boosts"][0.9]["P"] == pytest.approx(
        r["boosts"][0.9]["gamma"] * r["rest_1d"]["E0"] * 0.9, rel=1e-6)
    assert max(r["noisy"]["off_clock"]) < 2 * r["noisy"]["off_clock"][0]


def test_lump_stability_common_term_binds_g1_makes_the_clock():
    """Second variation of E - mu N on the lump in all six real
    components: one negative eigenvalue (l = 0), the phase zero at
    l = 0, the translation zero at l = 1, every relative direction
    positive. At g1 = 0 the lump still binds and its l = 0 kernel is
    five dimensional (phase and four rotations). The quartic is lowest
    on the clock. The carrier potential is nonnegative on every
    orientation iff m c^2 >= 2|mu0|. Gradient flow at fixed N from
    random data ends on the clock lump, q = 1 or its conjugate q = 2."""
    from probes import lump_stability

    r = lump_stability.run(report=_silent)
    assert r["hessian_error"] < 1e-5
    assert r["counts"] == {0: (1, 1), 1: (0, 1), 2: (0, 0)}
    assert min(r["relative_min"].values()) > 0.02
    assert r["relative_min"][0] == pytest.approx(0.0723, abs=2e-3)
    assert r["kernel_free"] == 5
    assert r["floor"] >= -0.61786 and r["floor_free"] == pytest.approx(-0.5, abs=1e-9)
    assert r["edge"] == pytest.approx(0.57261, abs=1e-5)
    assert {fl["mode"] for fl in r["flows"]} == {1, 2}
    for fl in r["flows"]:
        assert fl["share"] > 1 - 1e-8
        assert fl["E"] == pytest.approx(r["E_1d"], rel=2e-4)


def test_saturation_origin_is_the_depth_not_the_carrier():
    """The carrier's dispersion adds an attractive R^2 term to mu.  The
    static one loop sextic is attractive.  The zero point energy of the
    relative pair would set an equilibrium five decades below the
    measured depth floor. With the sextic the depth is
    sqrt(kappa)/(2 sqrt 3 a^2). Any total density saturating term gives
    the same stable count. The cubic norm sextic opens the relative
    sector at every kappa."""
    from probes import saturation_origin as so

    r = so.run(report=_silent)
    assert r["kinematic_R2"] < 0 and r["kinematic_R2"] == pytest.approx(-2 * so.A**2, rel=1e-2)
    assert np.all(r["hessian_eigs"] < 0) and r["sum_cubes"] < 0
    assert r["branch_error"] < 1e-6 and r["common_imag"] > 1.0
    assert r["lhy_check"] == pytest.approx(1.0, abs=1e-6)
    assert r["depth_lhy"] < 1e-2 and r["correction_at_floor"] < 1e-4
    assert r["depth_formula"] == pytest.approx(1.0, abs=1e-9)
    for f in r["forms"].values():
        assert abs(f["laue"]) < 1e-8
        assert f["l0"] == (1, 1) and f["l1"] == (0, 1) and f["relative_min"] > 0.01
    assert r["norm_l0"][0] >= 2 and r["norm_l1"][0] >= 2
    assert r["uniform_norm"] < -2 < 0 < r["uniform_total"]
    assert r["uniform_norm_kappa_0.25"] == pytest.approx(r["uniform_norm_kappa_4.0"], rel=1e-6)
    assert 0.0 < r["x_critical"] < 0.1


@pytest.mark.slow
def test_lump_internal_modes_soft_mode_does_not_radiate():
    """The 1D clock lump sits far from the Manakov point. Its relative
    internal modes lie below the continuum edge at 0.0393, 0.1726,
    0.2376 with first harmonics in the continuum 7, 2, 2 and no growing
    mode. Seeded small, the K = 2 mode loses a third of its lock deficit
    by t = 4000 and the K = 7 mode keeps it; a large soft excess drains
    and slows."""
    from probes import lump_internal_modes as im

    r = im.run(report=_silent)
    assert min(r["off_manakov"].values()) > 0.5
    assert r["growth"] < 1e-6 and r["clock_plane_max"] < 1e-6
    assert np.allclose(r["omegas"], [0.0393, 0.1726, 0.2376], atol=5e-4)
    assert np.allclose(r["omegas"], r["omegas_coarse"], atol=1e-4)
    assert r["harmonics"] == [7, 2, 2]
    soft, hard, large = r["soft small"], r["hard small"], r["soft large"]
    assert abs(soft["deficit"][-1] / soft["deficit"][0] - 1) < 1e-3
    assert hard["deficit"][-1] / hard["deficit"][0] < 0.75
    assert large["excess"][-1] / large["excess"][0] < 0.5
    drops = -np.diff(large["excess"]) / np.diff([125, 375, 750, 1500, 2500, 3500])
    assert drops[-1] < drops[0]
    assert np.allclose(r["bound_3d"], [0.0406, 0.1362], atol=1e-3)


def test_spinning_lump_exists_above_threshold():
    """Stationary spinning clock lumps psi_a = F(rho, z) e^{i m phi}
    w^a / sqrt 3: m = 1 is bound from N = 500, with mu/mu_co = 0.59,
    0.70, 0.78 at N = 500, 1000, 2000, the second and third on the
    stable side of the cubic quintic threshold 0.69; m = 2 is bound
    from N = 500."""
    from probes import spinning_lump as sl

    r = sl.run(report=_silent)
    tori = r["tori"]
    assert tori[(1, 200)]["E"] > 0 > tori[(1, 500)]["E"]
    assert tori[(2, 200)]["E"] > 0 > tori[(2, 500)]["E"]
    assert tori[(1, 500)]["mu_ratio"] < 0.69 < tori[(1, 1000)]["mu_ratio"] < tori[(1, 2000)]["mu_ratio"]
    assert tori[(1, 2000)]["ring"] > tori[(1, 1000)]["ring"] > tori[(1, 500)]["ring"] > 3.0


@pytest.mark.slow
def test_spinning_lump_persists_or_splits_in_three_components():
    """Full three component real time at 48^3 with noise in every
    component: the m = 1 ring at N = 500 splits by t = 130, at N = 2000
    it holds its shape and angular momentum with the clock locked."""
    from probes import spinning_lump as sl

    r = sl.run(report=_silent, three_d=True)
    small, large = r[("3d", 500)][-1], r[("3d", 2000)][-1]
    assert max(small["harm"]) > 0.3
    assert max(large["harm"]) < 0.05
    assert large["Lz"] / large["N"] == pytest.approx(1.0, abs=2e-3)
    assert large["lock"] > 0.999


@pytest.mark.slow
def test_knot_in_lump_unwinds_and_drains():
    """A unit hedgehog written into the relative spinor of the 3D clock
    lump on the self-bound sign flows at fixed N to the plain clock lump:
    the winding estimate goes to zero, the relative charge drains, and
    the energy ends at the clock lump's."""
    from probes import knot_in_lump as kl

    r = kl.run(report=_silent)
    for (n, amp, size), case in r.items():
        rows = case["rows"]
        assert abs(rows[0]["winding"]) > 0.8 and abs(rows[-1]["winding"]) < 0.05
        assert rows[-1]["relative"] < 1e-3 * rows[0]["relative"]
        assert rows[-1]["E"] == pytest.approx(case["E_clock"], abs=5e-3)
        assert all(b["E"] <= a["E"] + 1e-9 for a, b in zip(rows, rows[1:]))


@pytest.mark.slow
def test_lump_formation_on_the_self_bound_sign():
    """Real time formation in the empty field, one dimension, full
    three component dynamics on g0 < 0 < g1: a random cloud leaves a
    bound lump at the centre. The G2 term raises the phase lock from the
    cloud's value and it does not fall back, while at g1 = 0 the lump
    binds with its internal state frozen. The lock stays below one and
    the lump keeps a positive excess over the clock lump at its charge.
    On the gas sign the same cloud disperses."""
    from probes import lump_formation

    r = lump_formation.run(report=_silent)
    assert r["force_error"] < 1e-12 and r["closed_dN"] < 1e-6 and r["closed_dE"] < 1e-5
    for seed in (1, 2):
        m = r["main"][seed]
        assert m["late"]["E_w"] < 0 and m["late"]["off_centre"] < 2.0
        assert m["early"]["lock"] > m["start"] + 0.5
        assert m["late"]["lock"] > m["early"]["lock"] - 0.05
        assert m["late"]["lock"] < 1.0 and 0.0 < m["excess"] < 0.1
    first, second = r["free"]
    assert second["E_w"] < 0 and abs(second["lock"] - first["lock"]) < 1e-2
    assert r["gas"]["peak"] < 0.02 and r["gas"]["E_w"] > 0


def test_self_bound_vacuum_gas_sign_holds_nothing_self_bound_sign_binds():
    """On the gas sign (+,+) the clock has p = eps = a R^2 with
    a = g0/2 - g1/6 > 0, so it binds nothing without external pressure. At
    |g1/g0| = 1/sqrt2 only g0 < 0 < g1 keeps the clock as the phase
    ground state with a < 0. Adding a total density sextic gives a
    zero pressure point with eps - mu n = 0, a real Bogoliubov
    spectrum, the relative gaps unchanged, and c_rel/c_ph = 0.691."""
    from probes import self_bound_vacuum

    r = self_bound_vacuum.run(report=_silent)
    assert r["a_gas"] > 0
    assert all(abs(p - eps) < 1e-6 for _, eps, _, p in r["gas"])
    assert abs(r["p0"]) < 1e-8 and abs(r["interior_term"]) < 1e-8
    assert r["ground"] and r["stable"]
    assert np.max(abs(r["gaps"] - r["gaps_expected"])) < 1e-4
    assert r["ratio_liquid"] == pytest.approx(0.6905, abs=1e-3)
    assert r["ratio_gas"] == pytest.approx(0.6209, abs=1e-3)
