"""parameter_ledger.py: what the framework holds fixed, by kind, and the facts that tie it.

THE QUESTION. The substrate's "pixel" (healing length), "impedance"
(dimensionless coupling) and smallest quantum, and the ledger's scales:
which are units, which are derived or closure fixed, which are
identifications, and which numbers are genuinely free?

THE INVENTORY (status words are the repository's).
  Units: hbar, m, c, |g0|, and the anchor m_e ("SOLE DIMENSIONAL INPUT").
  Closure fixed: the four integers (wiring scan), the sign g0 < 0 (the
    self-bound gate), the liquid density (zero pressure), the integer
    windings (phase closure: the smallest quantum is one closed winding).
  Free continuous numbers: the depth rho xi^3, equivalently the sextic
    kappa ("the coefficient is the depth", saturation_origin.py O4),
    bounded only by data. And nu = mu/mc^2, a scanned regime number on
    which no gate depends at leading order.
  Identifications, substrate: g1/g0 = 1/sqrt 2, the Z3 phase alpha = 0,
    xi0 = C_EFF l_Pl = l_Pl/d10.
  Identifications, ledger (discrete): the boundary M_b read as the full
    Planck mass (with G anchored, t0 = Sigma/(2 pi) l_Pl^2 is the same
    identification, root.py), the trace 16, the h-dual metric, B/A's channel map,
    theta's slot, the per mode E[v^2], t0 = d11 l_Pl^2 (the diagnostic of
    sakharov_normalisation.py FACT 7), the scheme
    constant at mu*, the one loop truncation,
    the (1 - h10) vent, the CKM/PMNS maps, the bridge read as amount
    past the horizon, the declared premise T1 of
    notes/open_derivation_tasks.md (boundary.py bridge="amount": no share in b0,
    weight 1, light at the regulator;
    tests/probes/horizon_amount.py A2, A3). G_F tests C_EFF, the depth-2
    coefficient and N together (g_path_trials.py G4-G6), and each is an
    object the web already holds (tests/probes/frame_coefficients.py).
  Identifications read as objects the web holds, each alternative
    failing against data: C_EFF = 1/2, the depth-2 coefficient,
    w4 = 27/2 (tests/probes/frame_coefficients.py).
  Read from computations with a failing case: the vector contact term
    (the Weitzenbock term) and the amount's zero share, tested by
    1/alpha(0) (tests/probes/gravity_premises.py).
  Imports: H0, Omega_Lambda, m_H, Delta r_W, rho hat, cited QCD constants.

FACTS (the numbers below check each one).
  P1  Pixel, impedance and quanta per pixel are one number. In 3D,
      xi^2 = hbar^2/(m g rho), so the dimensionless coupling is
      g m rho^(1/3)/hbar^2 = (rho xi^3)^(-2/3): the impedance is the depth.
  P2  The probes run below the depth that alpha_G2 implies. At g0 = -1,
      g1/g0 = 1/sqrt 2, kappa = 1 the clock quartic is a = -0.618 and the
      depth sqrt(kappa)/(2 sqrt 3 a^2) = 0.756 in O4's convention. The
      committed sign's depth 3 of P9 is quoted in the range convention,
      2^(3/2) smaller than O4's, so it needs kappa = 125.9. Kappa = 15.7
      is depth 3 in O4's convention (channel_response.py X3).
      Classically kappa scales away. Quantum mechanically it is the depth.
  P3  Gravity does not see the depth. The Bogoliubov vacuum's
      entanglement per area in healing units is the same at mu = 1 and
      mu = 4 (record_area.py): the record area depends on the pixel size
      alone. The pixel is bounded below by half the carrier's Compton
      length, so a Planckian record implies a Planckian carrier mass
      (closed_event_scale.py S1, S2).
  P4  C_EFF carries two roles that the closure separates. As the
      electron's amplitude prefactor it is pinned at 1/2 by G_F. As
      xi0/l_Pl the record closure would place it at 1.918 (record_area.py
      R4). Gate C excludes that value for a vacuum regulator (R5: a
      medium filling space needs xi > 2.7e6 l_Pl). So the closure
      describes a record inside a lump and does not fix G's
      normalization. The same holds on the committed clock medium
      (clock_record.py K2-K4), where no reading of the record gives 0.5.
      The second role is an identification that rests on the first.
  P5  Two healing conventions are in use. The first is xi^2 = 1/mu
      (effective_action.py, sakharov_normalisation.py FACT 5 and
      FACT 6(c)). The second is xi^2 = 1/(4 mu), the screening range
      (2 c_s)^-2 on the gas sign (record_area.py, substrate_coupling.py).
      saturation_origin.py's 1/(2|a|R0) equals 1/(4|mu0|) on the clock,
      which is 2 range^2 there (channel_response.py X3). The
      identification xi0 = l_Pl/d10 is taken in the first convention by
      sakharov_normalisation.py and in the second by record_area.py.
      They differ by 2 in length and 8 in xi^3, the factor
      substrate_depth.py calls its one O(1) sensitivity.
  P6  The node event cost c is not computed in the repository
      (substrate_depth.py: "the effective theory has not computed"), so
      the depth floor is on the product c rho xi^3.

  P7  The inversion alpha_G2 = 1/(24 pi) -> depth is a menu until the
      gauge sector is derived as substrate excitations: the readings of
      "the substrate's coupling" span depth 0.32 to 655 (repo healing
      convention). They fall in two classes. Read as a loop correction
      (alpha or alpha/2pi equal to the quantum depletion or the LHY
      correction), the depth is 0.3 to 6.4 and the node event floor
      needs c of 30 to 600. Read as the tree level coupling (1/depth, or
      the medium's impedance g n^(1/3) = alpha, convention free), the
      depth is 75 to 655 and the floor needs c of 0.3 to 2.6. With c an
      O(1) action coefficient the floor admits the tree level class
      alone: g n^(1/3) = 1/(24 pi), depth 81.8 in the screening
      convention (654.7 for xi^2 = 1/mu), the floor met for c > 2.4.
  P8  The floor does not apply to the committed sign. It treats a
      proton as a pi_3 texture of a populated relative sector, which is
      the gas sign picture. The self-bound clock lump's vacuum manifold
      is a circle times two orientations, with trivial pi_3. A seeded
      knot drains (knot_in_lump.py K1, K2). What persists is held by
      charge and angular momentum conservation (spinning_lump.py), with
      no exp(-c rho xi^3) barrier. On the committed sign the depth is
      bounded only by classicality: the quantum depletion
      4.21e-3/(rho xi^3) is well below one (1.3 percent at the smallest
      reading). Both classes of P7 survive. Only the reading of alpha
      as a tree level coupling, not a loop correction, separates them.
      The depth is free in the committed theory.

  P9  A coupling defined by physics, not by normalization: the force the
      medium carries between two of its quanta, in units of hbar c_s, is
      1/(8 pi rho xi^3) for any equation of state (substrate_coupling.py).
      Read as alpha_G2 = 1/(24 pi), the depth is 3. The P7 split by
      class was a split by normalization: the physically defined
      exchange coupling is tree level and gives depth 3. At depth 3
      the gas sign node event floor would need c > 65, so it fails for
      any O(1) c. This is P8 seen from the other side: that floor
      belongs to the gas sign knot channel, not to the committed sign.
      Two qualifications (channel_response.py). First, the law is a
      one component statement. On the three component clock, a quantum
      coupled by d mu/dR gives depth 3, and one coupled by its own
      Hartree potential gives 1.418 times that, 4.25. This is a binary
      choice, not a continuum. Second, the depth is in the range convention
      xi = 1/(2 c_s). On the clock mu0 = c_s^2/2, so O4's convention
      xi^2 = 1/(4|mu0|) quotes it 2^(3/2) larger: range depth 3 is O4
      depth 8.49 and kappa = 126, not 15.7.

RESULT. Once units and closures are counted, the substrate holds one
free physical number, the depth, with no data floor on the committed
sign, plus a regime number nu on which no gate depends. Its other
content is two identifications (g1/g0, alpha). The link to the ledger is
one more identification, xi0/l_Pl. It shares its value with a G_F
pinned prefactor, but the record closure places it at 1.9 to 2.9.
Gravity's normalization is the boundary identification (root.py).
xi0/l_Pl describes a lump interior, and the smallest closed event ties
it to a Planckian carrier mass. Matter stability and the saturation
depend only on the depth. The carrier mass and nu enter no result,
because no carrier mass makes a lump a point particle
(closed_event_scale.py).

Usage: python3 tests/audit/parameter_ledger.py
"""

import math
import os
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_here), "probes"))
import _bootstrap  # noqa: E402,F401
import record_area as ra  # noqa: E402
import channel_response as cr  # noqa: E402


def run(report=print):
    out = {}
    # P1: impedance = depth^(-2/3), checked on random (g, rho) with hbar = m = 1
    rng = np.random.default_rng(0)
    g, rho = rng.uniform(0.1, 5, 5), rng.uniform(0.1, 5, 5)
    xi = 1 / np.sqrt(g * rho)
    out["P1"] = float(np.max(abs(g * rho**(1 / 3) - (rho * xi**3)**(-2 / 3))))
    # P2: the probes' depth
    g1 = 1 / math.sqrt(2)
    a = -(0.5 + g1 / 6)
    depth = 1 / (2 * math.sqrt(3) * a * a)
    assert abs(cr.medium(g1)["a"] - a) < 1e-15
    kappa_range3 = (3 * 2**1.5 * 2 * math.sqrt(3) * a * a)**2       # depth 3, range convention (P9)
    kappa_o4_3 = (3 * 2 * math.sqrt(3) * a * a)**2                  # depth 3, O4's convention
    out.update(a=a, depth_probe=depth, kappa_for_range_depth_3=kappa_range3, kappa_for_o4_depth_3=kappa_o4_3)
    # P3: record area in healing units at two chemical potentials
    def kappa_ent(mu):
        kp = np.linspace(1e-3, 12, 60) * math.sqrt(mu)
        s = np.array([ra.chain_entropy(x, 0.1 / math.sqrt(mu), 800, mu=mu) for x in kp])
        return float(np.trapezoid(kp * s, kp) / (2 * math.pi)) / (4 * mu)
    out["kappa_mu1"], out["kappa_mu4"] = kappa_ent(1.0), kappa_ent(4.0)
    # P4: the two roles of C_EFF
    r = ra.run(report=lambda *x: None)
    out["xi0_closure"], out["C_EFF"] = r["xi_closure"][1.0], 0.5
    # P5: the healing conventions on the clock: 1/(2|a|R0) against 1/(4|mu0|)
    R0 = 3 * abs(a) / 2
    mu0 = a * R0 / 2
    out["conv_match"] = abs(1 / (2 * abs(a) * R0) - 1 / (4 * abs(mu0)))
    al = 1 / (24 * math.pi)
    K = (16 * math.pi)**-1.5
    menu = {"depletion = alpha": K * 8 / (3 * math.sqrt(math.pi)) / al,
            "depletion = alpha/2pi": K * 8 / (3 * math.sqrt(math.pi)) / (al / (2 * math.pi)),
            "LHY = alpha/2pi": K * 128 / (15 * math.sqrt(math.pi)) / (al / (2 * math.pi)),
            "1/depth = alpha": 1 / al,
            "impedance g n^(1/3) = alpha": (4 * al)**-1.5,
            "impedance, xi^2 = 1/mu": al**-1.5}
    loop = [menu[k] for k in ("depletion = alpha", "depletion = alpha/2pi", "LHY = alpha/2pi")]
    tree = [menu[k] for k in ("1/depth = alpha", "impedance g n^(1/3) = alpha", "impedance, xi^2 = 1/mu")]
    out.update(menu=menu, loop_max=max(loop), tree_min=min(tree))
    depletion = {k: 4.21e-3 / v for k, v in menu.items()}
    out["depletion_max"] = max(depletion.values())
    assert out["depletion_max"] < 0.02
    assert max(loop) < 7 and min(tree) > 75 and abs(menu["impedance g n^(1/3) = alpha"] - 81.84) < 0.05
    assert out["P1"] < 1e-12
    assert abs(depth - 0.756) < 1e-3 and abs(kappa_range3 - 125.9) < 0.1 and abs(kappa_o4_3 - 15.7) < 0.1
    assert abs(out["kappa_mu1"] / out["kappa_mu4"] - 1) < 1e-6
    assert abs(out["xi0_closure"] / out["C_EFF"] - 3.84) < 0.05
    assert out["conv_match"] < 1e-12
    report(f"P1 impedance = depth^(-2/3): max deviation {out['P1']:.1e}")
    report(f"P2 probe depth {depth:.3f} at kappa = 1 (a = {a:.3f}). Depth 3 (range convention) needs "
           f"kappa = {kappa_range3:.1f}, depth 3 in O4's convention kappa = {kappa_o4_3:.1f}")
    report(f"P3 record area in healing units: {out['kappa_mu1']:.5f} at mu = 1, {out['kappa_mu4']:.5f} at mu = 4")
    report(f"P4 C_EFF = 1/2 (G_F pinned prefactor) against xi0/l_Pl = {out['xi0_closure']:.3f} (record closure)")
    report("P5 1/(2|a|R0) = 1/(4|mu0|) on the clock. The other convention, 1/mu, is 2x longer")
    report("P7 alpha -> depth menu: " + ", ".join(f"{k} {v:.2f}" for k, v in menu.items())
           + f". Loop class <= {max(loop):.1f}, tree class >= {min(tree):.1f}")
    report("RESULT: one free physical number in the substrate (the depth), a regime number nu "
           "no gate depends on, two substrate identifications, and xi0/l_Pl as the ledger link, "
           "which ties the regulator to a Planckian carrier mass that no gate depends on.")
    return out


if __name__ == "__main__":
    run()
