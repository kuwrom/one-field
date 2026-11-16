"""dust_gate_state.py: what the amount's state would need to pass the dust gate.

THE QUESTION. dust_gate.py closes the classical route: the carrier
carries no conserved current with a mean density in the empty field.
The remaining dynamical route is the amount's stress in a state, the
induced level.  For any ensemble of quanta the stress follows from the
quanta's dispersion and occupation, so the dust gate becomes a
condition on the state.  This probe states that condition.

THE BOUNDS (dark_sector.py, THE DUST GATE; arXiv:1601.05097, 99.7%):
|w| < 0.0024 and c_s^2 < 3.2e-6 at recombination.

STATEMENTS (the numbers below check each one).

  S1 (a massless amount cannot pass, in any state).  An isotropic
     ensemble of massless quanta has a traceless stress, w = 1/3
     exactly. A homogeneous massless mode gives w = 1. Frozen gradient
     energy gives w = -1/3.  A mixture crosses w = 0 only for an
     instant: radiation dilutes as a^-4 and gradients as a^-2, so at the
     crossing dw/dln a = -1/3, and w leaves the band within
     |Delta ln a| = 0.0072 of it.

  S2 (a massive amount passes only if cold).  For quanta of mass m with
     momenta p, w = <p^2/(3 E^2)> and, for a collisionless ensemble,
     c_s^2 is of order <v^2>/3.  The w bound needs v_rms < 0.085 c. The
     c_s^2 bound needs v_rms < 3.1e-3 c (about 930 km/s) at
     recombination.  A coherent massive mode passes for m much larger
     than the Hubble rate.

  S3 (never thermal).  A real scalar that shared the plasma's
     temperature and decoupled above every Standard Model threshold
     carries Delta N_eff = (4/7)(11/4)^(4/3) (g0/g_dec)^(4/3) = 0.027
     with g0 = 43/11 and g_dec = 106.75, so 182 of them carry 4.9,
     against N_eff = 2.99 +- 0.17 (Planck 2018, arXiv:1807.06209).  The
     amount couples through the common mode alone, where the
     scattering rate over the expansion rate is of order (T/M_Pl)^3.
     That ratio is below one at every temperature below the boundary,
     so the amount does not thermalize. Its occupation is an input of
     the state, not an output of thermal contact.

RESULT.  The induced route passes the gate only if the amount has a mass
and a velocity dispersion below 3e-3 c at recombination, or is a
coherent massive mode.  Mass and occupation are species data of the
amount, and the framework does not name them. The route reopens only
if a closure fixes them.  Light at the regulator (dark_sector.py, the
gravity ledger) means m far below M_b, not m = 0, so it does not
conflict with a cold amount.

Usage: python3 tests/probes/dust_gate_state.py     (instant)
"""

import math

W_MAX, CS2_MAX = 0.0024, 3.2e-6


def s1():
    w_rad, w_kin, w_grad = 1 / 3, 1.0, -1 / 3
    # mixture of radiation (a^-4) and frozen gradients (a^-2): w(x), x = rho_r/rho_g, dx/dln a = -2x
    w = lambda x: (x * w_rad + w_grad) / (x + 1)
    dwdlna_at_crossing = (w_rad - w_grad) * (-2 * 1.0) / (1 + 1.0)**2      # at x = 1
    window = W_MAX / abs(dwdlna_at_crossing)
    return dict(w_rad=w_rad, w_kin=w_kin, w_grad=w_grad, w_cross=w(1.0),
                dwdlna=dwdlna_at_crossing, window_dlna=window)


def s2():
    # monochromatic ensemble of speed v (units of c): w = v^2/3 exactly for p^2/(3E^2) = v^2/3
    v_w = math.sqrt(3 * W_MAX)
    v_cs = math.sqrt(3 * CS2_MAX)
    return dict(v_max_w=v_w, v_max_cs=v_cs, v_max_cs_kms=v_cs * 299792.458)


def s3(n=182, g0=43 / 11, g_dec=106.75, bound=0.17):
    per = (4 / 7) * (11 / 4)**(4 / 3) * (g0 / g_dec)**(4 / 3)
    return dict(per_scalar=per, total=n * per, sigmas=n * per / bound)


def run(report=print):
    a, b = s1(), s2()
    c = s3()
    assert abs(c["per_scalar"] - 0.027) < 5e-4 and 4.8 < c["total"] < 5.0 and c["sigmas"] > 25
    assert abs(a["w_cross"]) < 1e-12 and abs(a["dwdlna"] + 1 / 3) < 1e-12
    assert abs(a["window_dlna"] - 0.0072) < 1e-6
    assert abs(b["v_max_cs"] - 3.098e-3) < 1e-5 and abs(b["v_max_w"] - 0.0849) < 1e-4
    report(f"S1 massless amount: w = {a['w_rad']:.3f} (isotropic quanta), {a['w_kin']:.0f} (homogeneous mode), "
           f"{a['w_grad']:+.3f} (frozen gradients). A mixture crosses w = 0 with dw/dln a = {a['dwdlna']:+.3f} "
           f"and leaves the band within |d ln a| = {a['window_dlna']:.4f}")
    report(f"S2 massive amount: the w bound needs v_rms < {b['v_max_w']:.3f} c, the c_s^2 bound v_rms < "
           f"{b['v_max_cs']:.2e} c ({b['v_max_cs_kms']:.0f} km/s) at recombination")
    report(f"S3 never thermal: one decoupled real scalar carries Delta N_eff = {c['per_scalar']:.4f}, 182 carry "
           f"{c['total']:.2f}, {c['sigmas']:.0f} sigma of the Planck bound. Coupled through the common mode alone "
           "the amount does not thermalize below the boundary")
    report("RESULT: the induced route needs a mass and a cold state for the amount, species data the "
           "framework does not name. It reopens only if a closure fixes them.")
    return dict(s1=a, s2=b, s3=c)


if __name__ == "__main__":
    run()
