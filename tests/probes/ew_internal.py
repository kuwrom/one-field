"""
ew_internal.py: how much of the two declared EW imports can the
framework compute from its own state?

The two imports (world/refdata): dr_hat_W = 0.06937, rho_hat = 1.01016
(PDG 2024 Eq. 10.26).  The import-audit taxonomy classifies import
content as [1] data-driven measurement, [2] parameter dependence
(recalculable from framework values), [3] pure loop mathematics.
This probe computes the class-[2] content internally and states
exactly what remains in classes [1] and [3].

  rho_hat:  dominant content is the top-doublet term
            delta_rho_t = 3 G_F m_t^2 / (8 sqrt(2) pi^2),
            a function of two FRAMEWORK OUTPUTS (G_F from v_EW, m_t).
            Remainder: bosonic/higher loops (~0.0008) = class [3].

  dr_hat_W: dominant content is Delta_r0 = 1 - alpha(0)/alpha_hat(M_Z),
            driven by the DISPERSIVE hadronic vacuum polarization , 
            e+e- spectral DATA, class [1], same admissibility as PDG
            masses.  The framework's internal VP (free quarks +
            Lambda_conf cutoff) approximates it; the gap below
            measures exactly how far the internal computation reaches.
            Remainder: EW loop remainder (~0.0029) = class [3].

Usage: python3 ew_internal.py
"""

import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))


def _solved_web():
    """Run the derivation chain silently (once) and return the web."""
    with contextlib.redirect_stdout(io.StringIO()):
        import root
        import masses as _masses
        import mixing as _mixing
        import couplings as _couplings
        import gravity as _gravity
        R = root.derive()
        m = _masses.derive(R)
        mx = _mixing.derive(R, m)
        c = _couplings.derive(R, m)
        _gravity.derive(R, m, mx, c)
    return root.WEB


from root import PDG_EW, d10, d11, S_lepton  # noqa: E402


def alpha_em_at(mu_GeV, alpha_0, lepton_MeV, quark_MeV, M_Pl_GeV):
    """QED vacuum polarization running of alpha(0) to mu (internal VP:
    free fermions with the Lambda_conf cutoff for light quarks)."""
    Lambda_conf = (1.0 / d10) * M_Pl_GeV * math.exp(-S_lepton)
    Q_u_sq, Q_d_sq, N_c = (d10 / d11) ** 2, 1.0 / d11 ** 2, d11
    m_e, m_mu, m_tau = [m / 1e3 for m in lepton_MeV]
    m_u, m_d, m_s, m_c, m_b, m_t = [m / 1e3 for m in quark_MeV]
    fermions = [(m_e, 1, 1.0, False), (m_mu, 1, 1.0, False),
                (m_tau, 1, 1.0, False),
                (m_u, N_c, Q_u_sq, True), (m_d, N_c, Q_d_sq, True),
                (m_s, N_c, Q_d_sq, True), (m_c, N_c, Q_u_sq, True),
                (m_b, N_c, Q_d_sq, True)]
    d_inv = 0.0
    for mf, nc, qf2, is_q in fermions:
        if mf >= mu_GeV:
            continue
        mf_eff = Lambda_conf if (is_q and mf < Lambda_conf) else mf
        d_inv += (2.0 * nc * qf2) / (3.0 * math.pi) * math.log(mu_GeV / mf_eff)
    return 1.0 / (1.0 / alpha_0 - d_inv)


def ew_chain(alpha_phys, v_EW, sin2W):
    """M_W, M_Z from the web's sin2W node.  TWO declared imports:
    dr_hat_W, rho_hat (root.PDG_EW)."""
    A0 = math.sqrt(math.pi * alpha_phys) * v_EW
    M_W = A0 / (math.sqrt(sin2W) * math.sqrt(1.0 - PDG_EW['dr_hat_W']))
    M_Z = M_W / (math.sqrt(PDG_EW['rho_hat']) * math.sqrt(1.0 - sin2W))
    return {'sin2W': sin2W, 'A0': A0, 'M_W': M_W, 'M_Z': M_Z}


def run(report=print):
    s = _solved_web().state
    alpha = 1.0 / s["inv_alpha"]
    v = s["v_EW_GeV"]
    m_t = s["m_t"] / 1e3
    M_Pl = s["M_Pl_MeV"] / 1e3

    report("INTERNAL EW LOOP: retiring what can be retired")
    report("=" * 64)

    # G_F from the framework's own v_EW
    G_F = 1.0 / (math.sqrt(2.0) * v**2)
    report(f"  G_F (internal)      = {G_F:.7e}  "
           f"(measured 1.1663788e-05: "
           f"{1e6*(G_F/1.1663788e-5 - 1):+.2f} ppm)")

    # rho_hat: top-doublet term from framework outputs
    d_rho_t = 3.0 * G_F * m_t**2 / (8.0 * math.sqrt(2.0) * math.pi**2)
    rho_import = PDG_EW["rho_hat"]
    bosonic_rem = (rho_import - 1.0) - d_rho_t
    report(f"  delta_rho(top)      = {d_rho_t:+.5f}   [class 2: computed"
           f" from framework G_F, m_t]")
    report(f"  rho_hat import      = {rho_import:.5f} -> remainder "
           f"{bosonic_rem:+.5f} = bosonic/higher loops [class 3: pure"
           f" loop math]")

    # dr_hat_W: Delta_r0 from the internal VP vs the dispersive value
    ew = ew_chain(alpha, v, s["sin2W"])
    a_hat = alpha_em_at(
        ew["M_Z"], alpha,
        [s["m_e"], s["m_mu"], s["m_tau"]],
        [s["m_u"], s["m_d"], s["m_s"], s["m_c"], s["m_b"], s["m_t"]],
        M_Pl)
    dr0_int = 1.0 - alpha / a_hat
    dr0_disp = 0.06646                      # dispersive (audit value)
    dr_import = PDG_EW["dr_hat_W"]
    ew_rem = dr_import - dr0_disp           # ~0.0029, class [3]
    report(f"  Delta_r0 (internal VP) = {dr0_int:.5f}   vs dispersive "
           f"{dr0_disp:.5f}  (gap {dr0_int-dr0_disp:+.5f})")
    report(f"  dr_hat_W import     = {dr_import:.5f} -> EW remainder "
           f"{ew_rem:+.5f} [class 3]; hadronic VP itself is DATA"
           f" [class 1]")

    report("-" * 64)
    report("  VERDICT: rho_hat is retirable up to +0.0008 of pure loop")
    report("  mathematics; dr_hat_W is irreducibly data-bounded through")
    report("  the hadronic VP (class 1, same admissibility as the PDG")
    report("  masses), the internal VP reaches it to "
           f"{100*abs(dr0_int/dr0_disp-1):.1f}%.  The imports are")
    report("  hereby RECLASSIFIED (measurement + loop-math), not free")
    report("  parameters; full retirement needs a dispersive-grade")
    report("  internal VP, which is a data question, not a theory one.")
    return {"G_F": G_F, "d_rho_t": d_rho_t, "dr0_int": dr0_int}


if __name__ == "__main__":
    run()
