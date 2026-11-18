"""The smallest closed event against the regulator: what its count fixes.

THE QUESTION. The gravity lane's coefficient is 1/G = Sigma/(2 pi t0)
under a proper time cutoff t0 (sakharov_normalisation.py FACT 1). G is
the anchor inversion m_e -> M_Pl (root.py). So t0 = Sigma/(2 pi) l_Pl^2.
This is the boundary identification of root.py in proper time form:
4.350 l_Pl^2 on the counted ledger this probe reads (FACT 7's
reference), 4.165 on the engine's one list. A count of the states and
area of the smallest closed event could fix t0. clock_record.py K4
makes that count. The minimal lump of lump_in_empty_field.py holds
N_min = 164.1 x depth quanta on an area of 95.4 xi^2, and gives three
readings of xi/l_Pl (1.930, 4.54, 0.5). What does each reading fix?
Can the same lump be the particle of the matter reading and the pixel
of the regulator reading?

PREMISES. P1 the carrier of lump_in_empty_field.py, three complex
Klein-Gordon fields of mass m, the self-bound sign g0 = -1,
g1 = 1/sqrt 2, kappa = 1, and its lift to the relativistic energy
(lift). P2 record_area.py R4's premise that every ledger state is an
excitation of the substrate, so the record's pixel is the interior's
screening range xi = hbar/(2 m c_s) (the range convention of
clock_record.py). P3 G anchored, and Sigma the counted rows
(engine_sigma, FACT 7's reference), all massless at the regulator.
The carrier's own components vent E2(m^2 t0) of a massless row
(composite_pole.py K5).
P4 the records and readings of clock_record.py.

STATEMENTS (the numbers below check each one).

  S1 (the pixel is at least half a Compton length). On the clock at
     zero pressure the interior's sound speed obeys c_s^2/c^2 =
     2|mu0|/(m c^2) =: 2 nu0 exactly (lump_in_empty_field.py line
     261 at general m, c: the denominator m c^2 + 2 mu0 + R0 dmu/dR
     collapses to m c^2). lump_stability.py S3 keeps the empty field
     the global minimum iff m c^2 >= 2|mu0|, so nu0 <= 1/2, c_s <= c,
     and xi = hbar/(2 m c_s) >= lambda_C/2, lambda_C = hbar/(m c),
     with equality at the stability edge. The "0.757 c" of the lump
     probe is the value at its unit choice m c^2 = 3.49|mu0|. As the
     binding vanishes (nu0 -> 0) the pixel grows only as a power. So a
     ratio e^S0 between pixel and Compton length needs
     nu0 = 1/(8 e^(2 S0)) = 3.3e-40. Nearness to the self-binding
     threshold gives no exponential. (Vladimirov and Diakonov,
     arXiv:1208.1254, obtain their scale separation from a phase
     transition.)
  S2 (a Planck pixel needs a Planck carrier). Each reading xi/l_Pl
     bounds the carrier mass: m/M_Pl >= 1/(2 xi/l_Pl), larger by
     1/sqrt(2 nu0) below the edge. Per state record 1.930: m >= 0.259
     M_Pl. One record per quantum 4.54: m >= 0.110 M_Pl. The
     identification xi0 = l_Pl/d10 = 0.5 is exactly xi = hbar/(2 m c_s)
     at m = M_Pl, c_s = c. The carrier's own components then vent
     E2(m^2 t0) = 0.47, 0.82 and 2.1e-3 of a massless row. At the
     identification the carrier quanta sit above the regulator. The
     smallest closed event has mass N_min x 0.964 x m, from 52 to 474
     M_Pl at depth 3. In the O4 convention (sqrt 2 longer pixel) the
     per state bound is 0.366 M_Pl: Planckian in either convention.
  S3 (the electron is not the smallest closed event). If the minimal
     lump had mass m_e, the carrier would be m = m_e/(N_min x 0.964) =
     1.08e-3 MeV, 8.8e-26 M_Pl, its interior pixel at least 5.7e24
     l_Pl, and the interior's record per area 8e-52 of 1/(4 l_Pl^2).
  S4 (no reader, so no count fixes t0). No gate reads m or nu0
     (parameter_ledger.py). The matter reading of the lump (m a unit)
     and the regulator reading (m Planckian) need carrier masses 3e24
     apart. So the count of the smallest closed event fixes t0 only by
     fixing m. That trades one number for one number: an
     identification.
     The composite metric literature has the same
     structure, M_P = m [N Gamma(1 - D/2)/(6 (4 pi)^{D/2})]^{1/(D-2)}
     with the regulator a place holder (Carone, Claringbold and Vaman
     arXiv:1710.09367 Eq. 13; Carone, Erlich and Vaman
     arXiv:1610.08521).
  S5 (the lump is no point particle). At the edge the minimal lump's
     half radius is r_half = 2.755 xi (clock_record.py K4) >= 1.378
     lambda_C, and its energy is E = N_min x 0.964 m c^2. So it spans
     r_half E/(hbar c) >= 654 of its own Compton lengths at depth 3
     (218 x depth), for any m. The electron is pointlike below
     2.8e-7 of its own Compton length: r_e < 1.1e-19 m at LEP2
     (arXiv:1302.3415 Eq. 3.9). An electron mass minimal lump has a
     half density radius of 2.5e-10 m.
  S6 (a Planckian lump is inside its horizon). With G anchored,
     r_s/r_half = 689 (m/M_Pl)^2 at the edge, and more below it. The
     readings of S2 give 47, 8.3 and 689. A weak field outside a
     minimal lump needs m < 0.038 M_Pl. That is below every reading
     that ties the record to the regulator.

RESULT. The lump interior's pixel is bounded below by half the
carrier's Compton length. So every reading of the record ties the
regulator to a Planckian carrier mass that nothing reads. An electron
mass lump records nothing of 1/(4G), so the smallest closed event
does not set t0. t0 is the boundary identification of
root.py in proper time form, and G stays with the anchor. No carrier
mass makes the minimal lump a particle (S5, S6), so the carrier mass
is not a leftover free number. The lump is read in carrier units. It
is evidence of closure, not the ledger's point particle.

Usage: python3 tests/probes/closed_event_scale.py     (about twenty seconds)
"""

import math
import os
import sys

from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import clock_record as CR  # noqa: E402
import lump_in_empty_field as L  # noqa: E402
from root import C_EFF, M_E_ANCHOR_MEV, M_Pl_MeV, S_lepton  # noqa: E402
from sakharov_normalisation import engine_sigma  # noqa: E402

NU0_EDGE = 0.5                       # lump_stability S3: m c^2 >= 2|mu0|
DEPTH = 3.0                          # range convention depth, channel_response.py
HBARC_MEV_M = 197.3269804e-15        # hbar c in MeV m
R_E_LEP2_M = 1.1e-19                 # electron radius bound, arXiv:1302.3415 Eq. 3.9


def cs2_over_c2(c, m=1.0):
    """lump_in_empty_field line 261 at general (m, c): R0 dmu/dR over m c^2 + 2 mu0 + R0 dmu/dR."""
    rdmu = -L.A * L.R0               # = -2 mu0 on the clock at zero pressure
    return rdmu / (m * c * c + 2 * L.MU0 + rdmu)


def xi_over_lambda_c(nu0):
    """Range convention pixel hbar/(2 m c_s) over the Compton length hbar/(m c)."""
    return 1 / (2 * math.sqrt(2 * nu0))


def run(report=print):
    out = {}
    sigma = engine_sigma()[0]
    t0 = sigma / (2 * math.pi)                                   # l_Pl^2, G anchored
    out["sigma"], out["t0"] = sigma, t0
    assert abs(t0 - 4.350) < 2e-3

    # S1: the sound speed identity and the pixel bound.
    err = max(abs(cs2_over_c2(c) - 2 * (-L.MU0) / (c * c)) for c in (1.0, math.sqrt(-2 * L.MU0), 2.0))
    out["cs2_identity_error"] = err
    assert err < 1e-12 and abs(math.sqrt(cs2_over_c2(1.0)) - 0.75671) < 1e-5
    nu0_unit = -L.MU0                                            # the probe's m = c = 1 choice
    nus = (NU0_EDGE, nu0_unit, 0.1, 0.01)
    out["xi_over_lambda_c"] = {nu: xi_over_lambda_c(nu) for nu in nus}
    assert out["xi_over_lambda_c"][NU0_EDGE] == 0.5 and all(v >= 0.5 for v in out["xi_over_lambda_c"].values())
    report(f"S1 c_s^2/c^2 = 2|mu0|/mc^2 (identity error {err:.1e}); 0.757 c is the unit choice "
           f"mc^2 = {1 / nu0_unit:.2f}|mu0|. Pixel/Compton = 1/(2 sqrt(2 nu0)) = "
           + ", ".join(f"{v:.3f} at nu0 = {nu:.3f}" for nu, v in out["xi_over_lambda_c"].items()) + ".")

    nu_hier = 1 / (8 * math.exp(2 * S_lepton))                   # xi/lambda_C = e^S0
    assert abs(xi_over_lambda_c(nu_hier) / math.exp(S_lepton) - 1) < 1e-9 and 3e-40 < nu_hier < 3.5e-40
    out["nu0_for_e_S0"] = nu_hier
    report(f"   a pixel/Compton ratio e^S0 needs nu0 = 1/(8 e^(2 S0)) = {nu_hier:.1e}: a power, no exponential.")

    # S2: each reading of xi/l_Pl bounds the carrier mass.
    cr = CR.run(report=lambda *a, **k: None)
    k0, krel = cr["kappa_0"], cr["kappa_rel"]
    n_min = cr["quanta_per_depth"] * DEPTH
    d = L.lump(0.13624 * L.MU0)
    c_edge = math.sqrt(-2 * L.MU0)                                # m c^2 = 2|mu0| at m = 1
    frac = L.lift(d, c_edge)["E"] / (d["N"] * c_edge**2)         # E_min/(N m c^2) at the edge
    assert 0.95 < frac < 0.98
    readings = {"per state record": cr["readings"]["per state Bogoliubov"],
                "one record per quantum": cr["readings"]["one record per quantum, depth 3.0"],
                "identification l_Pl/d10": C_EFF}
    out["readings"] = {}
    for name, xi in readings.items():
        m_min = 1 / (2 * xi)                                      # M_Pl, at the edge
        row = dict(xi_lpl=xi, m_min=m_min, m2t0=m_min**2 * t0, E2=float(expn(2, m_min**2 * t0)),
                   E_min_over_MPl=n_min * frac * m_min,
                   m_at_nu0={nu: m_min / math.sqrt(2 * nu) for nu in nus})
        out["readings"][name] = row
    r_ps, r_q, r_id = (out["readings"][k] for k in readings)
    assert abs(r_ps["m_min"] - 0.2590) < 5e-3 and abs(r_q["m_min"] - 0.110) < 2e-3
    assert r_id["m_min"] == 1.0 and r_id["E2"] < 3e-3 and 0.4 < r_ps["E2"] < 0.5
    out["per_state_m_min_o4"] = r_ps["m_min"] * math.sqrt(2)
    report(f"S2 t0 = {t0:.3f} l_Pl^2, N_min = {n_min:.1f} at depth {DEPTH:.0f}, E_min/(N m c^2) = {frac:.3f} at the edge;")
    for name, r in out["readings"].items():
        report(f"   {name}: xi = {r['xi_lpl']:.3f} l_Pl -> m >= {r['m_min']:.3f} M_Pl, m^2 t0 = {r['m2t0']:.3f}, "
               f"E2 = {r['E2']:.2e}, E_min = {r['E_min_over_MPl']:.0f} M_Pl. At nu0 = 0.01, m = {r['m_at_nu0'][0.01]:.2f} M_Pl.")
    report(f"   O4 convention: per state bound {out['per_state_m_min_o4']:.3f} M_Pl.")

    # S3: the electron as the smallest closed event.
    m_carrier = M_E_ANCHOR_MEV / (n_min * frac)                   # MeV
    xi_e = (M_Pl_MeV / m_carrier) / 2                             # l_Pl, at the edge
    out["electron"] = dict(binding_fraction=frac, m_carrier_MeV=m_carrier, m_over_MPl=m_carrier / M_Pl_MeV,
                           xi_int_lpl=xi_e, record_over_quarter=4 * (k0 + krel) / xi_e**2)
    assert 1e24 < xi_e < 1e25 and out["electron"]["record_over_quarter"] < 1e-50
    report(f"S3 electron as the minimal lump: carrier m = {m_carrier:.2e} MeV = {m_carrier / M_Pl_MeV:.1e} M_Pl, "
           f"pixel >= {xi_e:.1e} l_Pl, interior record {out['electron']['record_over_quarter']:.1e} of 1/(4 l_Pl^2).")

    # S4: the dichotomy.
    out["dichotomy"] = r_ps["m_min"] * M_Pl_MeV / m_carrier
    assert out["dichotomy"] > 1e24
    report(f"S4 the two readings of the lump need carrier masses {out['dichotomy']:.1e} apart. No gate reads m or nu0.")

    # S5: the minimal lump in its own Compton lengths, at any carrier mass.
    r_half_xi = cr["r_half_xi"]
    own = r_half_xi * 0.5 * n_min * frac                          # r_half E/(hbar c) at the edge
    lep2 = R_E_LEP2_M / (HBARC_MEV_M / M_E_ANCHOR_MEV)
    size_e = r_half_xi * 0.5 * HBARC_MEV_M / m_carrier            # the electron mass minimal lump, m
    out.update(own_compton=own, own_compton_per_depth=own / DEPTH, electron_bound=lep2, electron_lump_m=size_e)
    assert abs(own - 654) < 2 and lep2 < 3e-7 and abs(size_e - 2.53e-10) < 2e-12
    report(f"S5 the minimal lump spans >= {own:.0f} of its own Compton lengths ({own / DEPTH:.0f} x depth) at any m; "
           f"the electron is below {lep2:.1e} of its own (LEP2). An electron mass lump has r_half = {size_e:.2e} m.")

    # S6: a Planckian lump against its own Schwarzschild radius.
    horizon = 4 * n_min * frac / r_half_xi                        # r_s/r_half per (m/M_Pl)^2 at the edge
    rs = {name: horizon * row["m_min"]**2 for name, row in out["readings"].items()}
    out.update(horizon_coefficient=horizon, rs_over_rhalf=rs, m_weak_field=horizon**-0.5)
    assert abs(horizon - 689) < 2 and min(rs.values()) > 1 and horizon**-0.5 < min(r["m_min"] for r in out["readings"].values())
    report(f"S6 r_s/r_half = {horizon:.0f} (m/M_Pl)^2 at the edge: " + ", ".join(f"{k} {v:.1f}" for k, v in rs.items())
           + f"; a weak field outside needs m < {horizon**-0.5:.3f} M_Pl.")
    report("RESULT: every reading of the record ties t0 to a Planckian carrier mass that nothing reads. The "
           "electron mass lump records nothing of 1/(4G). The route from the smallest closed event to t0 is "
           "closed (t0 is the boundary identification). No carrier mass makes the lump a point particle.")
    return out


if __name__ == "__main__":
    run()
