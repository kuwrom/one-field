"""
Dark sector: Ω_DM/Ω_b and cosmological density bookkeeping.

The matter ratio uses the algebraic bridge data. The density display uses
the observational H₀ anchor and distinguishes the Friedmann critical
density from the current dark-energy density.

Two structures are reported with different statuses:

1. Bridge self-interference → DM/baryon ratio
   The bridge (h=1, c_coset=0, D²=1) uses the single normalized response
   R(1)=1/(2π), the same factor
   that converts α_alg = π/512 to α(0) = 1/137.035999050.
     f_baryonic = 1/(2π)     →  Ω_DM/Ω_b = 2π−1 ≈ 5.283
   The registered matter ledger contains one bridge response; no separate
   materialization coefficient is introduced.

   Equivalently: Ω_m/Ω_b = 2π (total matter/baryon ratio).

   WHY THE SAME FACTOR GOVERNS BOTH DOMAINS: the 1/(2π) is not
   "borrowed" from the EM sector.  It is the bridge's canonical
   closed-loop response: the weight sum h(7,G₂) + h(26,F₄) = 1 fixes
   h, and the framework's chosen unit-normalized current basis defines
   R(h)=h/(2π) once for every consumer. Thus R(1)=1/(2π). It is not
   the raw integral of the singular separated-point two-point kernel. This
   fraction answers a single question: what fraction of the bridge
   amplitude closes locally (i.e. contributes to a particular
   channel's observable)?

   In the EM channel, "closes locally" means: reaches the physical
   coupling. At depth 1: 1/α ≈ (512/π)(1 − 1/(2π)). At web depth 4
   the self-consistent cubic gives 1/α(0) = 137.035999050 exactly
   (root.py).

   In the matter channel, "closes locally" means: materializes as
   localized (baryonic) matter.  The standing waves that form in the
   relative sector carve venting channels into the common mode
   (gravity.py).  Of the total substrate energy circulating through
   the bridge, the fraction 1/(2π) localizes as particle-like
   excitations (baryons). The remainder (2π−1)/(2π) stays as
   delocalized bridge coherence, coupling only through the common
   mode (PvP = 0 → gravitational coupling only → dark matter).

   Both readings use the same normalized response of the same operator.
   The bridge is one object connecting two faces, so no independent
   channel coefficient is introduced. The identity c_coset=0 removes an
   independent unitary coset Virasoro sector; the registered Web controls
   later responses.

   WHAT ENDS THIS.  Ω_DM/Ω_b = 2π−1 = 5.2832 is exact in the registered
   matter ledger, with no coefficient to adjust.  It sits
   between Planck (−1.3σ) and ACT (+1.2σ); both datasets converging
   away from 5.283 kills it, and CMB-S4-class data decides.
   The dark sector is bridge coherence, not a particle: it couples
   only through the common mode (PvP = 0), so all direct, indirect
   and collider dark-matter searches return null indefinitely.  Any
   confirmed non-gravitational dark-matter detection kills it.

2. Cosmological density scale
   Volovik equilibrium: ρ_vac = 0 at the self-bound point
   (imported thermodynamic relation).
   Friedmann critical-density identity:
     ρ_crit = (3/8π) M_Pl² H₀²
   Current dark-energy density:
     ρ_Λ = Ω_Λ ρ_crit
   Jacobson's thermodynamic derivation leaves Λ as an integration
   constant; CKN saturation is an identification hypothesis, not a
   derivation of the observed value.

Note on Z₃ coupling eigenvalues:
   The Z₃ substrate has coupling eigenvalues λ₀ = g₀(1+√d₁₀) and
   λ₁ = g₀(1−1/√d₁₀), satisfying λ₀+2λ₁ = d₁₁g₀ (Sugawara identity).
   These characterise COUPLING STRENGTHS in the Madelung system, not
   energy-density fractions.

Reference:
    Kahsay, Kibrom Kidane (2026). Trilogy.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import math

from . import wiring

from .root import (d10, d11,
                   h_bridge, c_coset,
                   M_Pl_GeV, RHO_CRIT_MAPPED, PDG_COSMO, pct)


# ═══════════════════════════════════════════════════════════════════════
#  Observational reference values (comparison only, not inputs)
#  Centralized in audit/reference.py, imported via root like every
#  other reference block, inside the blindness firewall's
#  perturbation set.
# ═══════════════════════════════════════════════════════════════════════

from .root import PLANCK_2018, ACT_DR6  # noqa: E402


# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive(R, grav_data):
    """
    Evaluate the algebraic matter ratio and cosmological density map.

    R         : dict from root.derive()
    grav_data : dict from gravity.derive()

    DAG POSITION: both parameters are accepted for the runner's
    uniform signature and neither is consumed; everything here needs
    only root-level constants (π, d₁₀, h_bridge, ρ_crit), so this
    module's place in the dependency graph is directly after root.
    run.py keeps it late for narrative order only.
    """
    print("\n" + "=" * 78)
    print("  DARK SECTOR: algebraic DM/baryon ratio + density bookkeeping")
    print("=" * 78)

    sd = math.sqrt(d10)

    # Z₃ coupling eigenvalue structure (algebraic) ########
    #
    # g₁/g₀ = 1/√d₁₀  (the G₂ identification on the Z₃ NLS substrate.
    # see masses.py substrate section, same g₀, g₁ as the BdG layer)
    # λ₀ = g₀(1+√d₁₀)   [common mode → geometry]
    # λ₁ = g₀(1−1/√d₁₀) [relative modes → matter]
    #
    # λ₀+2λ₁ = d₁₁g₀ exactly for d₁₀=2 (Sugawara identity d₁₁²−1=4d₁₀).

    lam0 = 1.0 + sd                    # common mode
    lam1 = 1.0 - 1.0/sd               # relative mode
    denom = lam0 + 2.0 * lam1          # = d₁₁ = 3 exactly
    assert abs(denom - float(d11)) < 1e-10, f"denom = {denom}, expected {d11}"

    print(f"\n  Z₃ coupling eigenvalue structure:")
    print(f"    g₁/g₀ = 1/√d₁₀ = 1/√{d10}  [G₂ identification]")
    print(f"    λ₀ = 1+√{d10} = {lam0:.6f}  [common → geometry]")
    print(f"    λ₁ = 1−1/√{d10} = {lam1:.6f}  [relative → matter]")
    print(f"    λ₀+2λ₁ = d₁₁ = {denom:.1f}  [Sugawara identity]")
    print(f"    (Coupling strengths, not energy-density fractions)")

    # Bridge self-interference → baryon/DM split #########
    #
    # h_bridge = 1 (dimension one); c_coset = 0 removes the coset sector
    # D²_local = 1 (Lagrangian condensation)
    # f_baryonic = h_bridge/(2π) = 1/(2π)
    # Ω_DM/Ω_b = (1 − f_b)/f_b = 2π − 1
    #
    # The registered matter ledger contains one bridge response.
    # No external inputs or adjustable coefficients enter.

    f_baryonic = wiring.normalized_current_response(h_bridge)  # 1/(2π)
    f_dark = 1.0 - f_baryonic
    DM_baryon_ratio = f_dark / f_baryonic              # = 2π − 1

    # Verify exact algebraic form
    assert abs(DM_baryon_ratio - (2.0 * math.pi - 1.0)) < 1e-10

    # the unread remainder as a graph node ############
    from .root import WEB, Ledger
    WEB["dark_ratio"] = Ledger("dark_ratio", 2.0 * math.pi - 1.0, "mul")
    WEB.solve()          # dark node joins the solved web state
    assert abs(WEB.state["dark_ratio"] / DM_baryon_ratio - 1.0) < 1e-14

    # Total matter/baryon ratio (algebraic identity)
    matter_baryon_ratio = 1.0 / f_baryonic             # = 2π

    print(f"\n  Bridge self-interference → matter split:")
    print(f"    h_bridge = {float(h_bridge)},  c_coset = {float(c_coset)}")
    print(f"    f_baryon = 1/(2π) = {f_baryonic:.6f}")
    print(f"    Ω_DM/Ω_b = 2π−1 = {DM_baryon_ratio:.6f}")
    print(f"    Ω_m/Ω_b  = 2π   = {matter_baryon_ratio:.6f}")
    print(f"    (Exact: c_coset = 0 closes the bridge response)")

    # Comparison with CMB observations ##############
    #
    # The DM/baryon MASS ratio Ω_c h²/Ω_b h² is H₀-independent and
    # measured directly from the CMB power spectrum.  No external
    # conversion factors needed, just a ratio of observables.

    Omega_c_h2_obs = PLANCK_2018['Omega_c_h2']
    Omega_b_h2_obs = PLANCK_2018['Omega_b_h2']
    sigma_c = PLANCK_2018['Omega_c_h2_err']
    sigma_b = PLANCK_2018['Omega_b_h2_err']

    DM_b_planck = Omega_c_h2_obs / Omega_b_h2_obs
    sigma_ratio_rel = math.sqrt((sigma_c / Omega_c_h2_obs)**2
                                + (sigma_b / Omega_b_h2_obs)**2)
    sigma_DM_b = DM_b_planck * sigma_ratio_rel

    # ACT DR6 ratio
    DM_b_ACT = ACT_DR6['Omega_c_h2'] / ACT_DR6['Omega_b_h2']
    sigma_ACT_rel = math.sqrt((ACT_DR6['Omega_c_h2_err'] / ACT_DR6['Omega_c_h2'])**2
                              + (ACT_DR6['Omega_b_h2_err'] / ACT_DR6['Omega_b_h2'])**2)
    sigma_DM_b_ACT = DM_b_ACT * sigma_ACT_rel

    # Pulls and errors
    pull_planck = (DM_baryon_ratio - DM_b_planck) / sigma_DM_b
    pull_ACT = (DM_baryon_ratio - DM_b_ACT) / sigma_DM_b_ACT
    err_planck = pct(DM_baryon_ratio, DM_b_planck)
    err_ACT = pct(DM_baryon_ratio, DM_b_ACT)

    # Total matter/baryon ratio from CMB
    m_b_planck = (Omega_c_h2_obs + Omega_b_h2_obs) / Omega_b_h2_obs
    err_m_b = pct(matter_baryon_ratio, m_b_planck)

    print(f"\n  Comparison with CMB observations (ratio-to-ratio, no conversion):")
    print(f"    Planck 2018:  Ω_c h²/Ω_b h² = {DM_b_planck:.4f} ± {sigma_DM_b:.4f}")
    print(f"    ACT DR6:      Ω_c h²/Ω_b h² = {DM_b_ACT:.4f} ± {sigma_DM_b_ACT:.4f}")
    print(f"    Prediction:   2π−1           = {DM_baryon_ratio:.4f}")
    print(f"    Planck pull:  {pull_planck:+.2f}σ  ({err_planck:+.2f}%)")
    print(f"    ACT pull:     {pull_ACT:+.2f}σ  ({err_ACT:+.2f}%)")
    print(f"    Sits between: {abs(DM_baryon_ratio - DM_b_planck):.3f} from Planck, {abs(DM_baryon_ratio - DM_b_ACT):.3f} from ACT")
    print(f"\n    Ω_m/Ω_b = 2π = {matter_baryon_ratio:.4f}  (Planck: {m_b_planck:.4f}, {err_m_b:+.2f}%)")

    # Cosmological density bookkeeping
    #
    # H₀ is an observational scale anchor.  The code evaluates the
    # Friedmann identity ρ_crit = (3/8π) M_Pl²H₀².
    #
    # Step 1, Volovik equilibrium:
    #   The imported Gibbs-Duhem relation gives ρ_vac(eq) = 0 for a
    #   self-bound condensate.
    #
    # Step 2, Jacobson-Clausius:
    #   Jacobson's local thermodynamic argument recovers Einstein's
    #   equation of state but leaves Lambda as an integration constant.
    #   It does not determine rho_Lambda.
    #
    # Step 3, CKN saturation:
    #   Saturating the UV/IR bound is an identification hypothesis.
    #
    rho_crit_mapped = RHO_CRIT_MAPPED
    rho_naive = M_Pl_GeV**4
    critical_density_suppression = math.log10(rho_crit_mapped / rho_naive)

    # Current epoch comparison
    Omega_Lambda = PDG_COSMO['Omega_Lambda']
    rho_Lambda_current = Omega_Lambda * rho_crit_mapped

    print(f"\n  Cosmological density bookkeeping:")
    print(f"    Step 1: Volovik equilibrium map, ρ_vac(eq) = 0")
    print(f"    Step 2: Jacobson leaves Lambda as an integration constant")
    print(f"    Step 3: CKN saturation is an identification hypothesis")
    print(f"    ρ_crit(H₀)       = {rho_crit_mapped:.2e} GeV⁴")
    print(f"    ρ_Λ(current)     = {rho_Lambda_current:.2e} GeV⁴  "
          f"(Ω_Λ = {Omega_Lambda})")
    print(f"    ρ_naive = M_Pl⁴ = {rho_naive:.1e} GeV⁴")
    print(f"    Critical-density ratio: 10^{critical_density_suppression:.0f}")

    # # Status ###########################

    print(f"\n  Status:")
    print(f"    Ω_DM/Ω_b = 2π−1:  {abs(err_planck):.2f}% from Planck ({abs(pull_planck):.2f}σ)  [algebraic]")
    print(f"    Dark-energy density status: current ρ_Λ uses measured Ω_Λ; "
          "CKN is an IH")

    return {
        # Z₃ coupling eigenvalue structure
        'lam0': lam0,
        'lam1': lam1,
        # Bridge self-interference (exact, c_coset=0)
        'f_baryonic': f_baryonic,
        'f_dark_matter': f_dark,
        'DM_baryon_ratio': DM_baryon_ratio,       # 2π−1
        'matter_baryon_ratio': matter_baryon_ratio, # 2π
        'DM_baryon_ratio_obs': DM_b_planck,
        'err_DM_baryon': err_planck,
        'pull_DM_baryon': pull_planck,
        # Cosmological densities
        'rho_crit_mapped': rho_crit_mapped,
        'rho_Lambda_current': rho_Lambda_current,
        'critical_density_suppression': critical_density_suppression,
    }
