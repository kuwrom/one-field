"""
Layer 6b -- Neutrino Masses from the F₄ Neutral Sector.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Structural prediction:
    The 26 of F₄ under Spin(8) decomposes as
        26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ 1 ⊕ 1
    The three 8's carry charged fermions (leptons, up-quarks, down-quarks)
    via Spin(8) triality.  The two singlets accommodate exactly two
    right-handed neutrinos.

    With only 2 right-handed neutrinos, the 3×3 light neutrino mass
    matrix from the Type-I seesaw has rank ≤ 2.  This forces one
    eigenvalue to vanish:
        m₁ = 0    (normal ordering)

    The remaining masses are fixed by the measured oscillation
    splittings alone -- no additional tunable inputs:
        m₂ = √(Δm²₂₁)
        m₃ = √(Δm²₃₁)

    This is a conditional prediction: it follows IF the F₄ branch
    of the emergence tree accommodates the neutral sector, which is
    the natural and minimal interpretation.

NuFit reference:
    Esteban et al., JHEP 2009, 178 (2020), arXiv:2007.14792;
    NuFit 6.0 (2024), http://www.nu-fit.org/
"""

import math
from .constants import NUFIT_OSCILLATION
from .formatting import H, S, box


def derive(alg: dict):
    """
    Derive neutrino masses from the F₄ singlet structure.

    Parameters
    ----------
    alg : dict from algebra.derive()  -- provides the 26 decomposition

    Returns
    -------
    dict with keys:
        m1_meV, m2_meV, m3_meV     : individual masses in meV
        sum_m_nu_meV               : Σmν in meV
        ordering                    : 'normal'
        n_RH                        : number of right-handed neutrinos (2)
        dm2_21, dm2_31              : mass-squared splittings used (eV²)
    """

    H("LAYER 6b:  NEUTRINO MASSES -- F₄ NEUTRAL SECTOR")

    # ── 6b.1  F₄ singlet structure ──────────────────────────────────

    S("6b.1  Right-handed neutrino count from F₄ decomposition")

    dim_F4_fund = alg['dim_F4_fund']  # = 26

    # Spin(8) triality decomposition of the 26 of F₄
    n_8v = 8   # vector     → charged leptons
    n_8s = 8   # spinor     → up-type quarks
    n_8c = 8   # co-spinor  → down-type quarks
    n_singlets = dim_F4_fund - n_8v - n_8s - n_8c  # = 2

    print(f"  F₄ fundamental: dim = {dim_F4_fund}")
    print(f"  Spin(8) decomposition:  26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ {n_singlets}·1")
    print(f"")
    print(f"  8_v  → charged leptons   (e, μ, τ)")
    print(f"  8_s  → up-type quarks    (u, c, t)")
    print(f"  8_c  → down-type quarks  (d, s, b)")
    print(f"  {n_singlets}·1  → right-handed neutrinos")
    print(f"")
    print(f"  Number of right-handed neutrinos: N_R = {n_singlets}")

    assert n_singlets == 2, f"Expected 2 singlets, got {n_singlets}"

    # ── 6b.2  Rank-2 seesaw → m₁ = 0 ───────────────────────────────

    S("6b.2  Rank-2 Type-I seesaw mechanism")

    print(f"  Type-I seesaw: M_light = −M_D · M_R⁻¹ · M_D^T")
    print(f"")
    print(f"  M_D is 3×{n_singlets} (3 LH flavours, {n_singlets} RH neutrinos)")
    print(f"  M_R is {n_singlets}×{n_singlets} (Majorana mass matrix)")
    print(f"  M_light is 3×3 with rank ≤ {n_singlets}")
    print(f"")
    print(f"  rank(M_light) ≤ {n_singlets}  →  at least {3 - n_singlets} zero eigenvalue(s)")
    print(f"  → m₁ = 0  (normal ordering)")

    # ── 6b.3  Neutrino masses from oscillation data ─────────────────

    S("6b.3  Masses from m₁ = 0 + oscillation splittings")

    dm2_21, dm2_21_err = NUFIT_OSCILLATION['dm2_21']  # eV²
    dm2_31, dm2_31_err = NUFIT_OSCILLATION['dm2_31']  # eV²

    # With m₁ = 0:
    #   Δm²₂₁ = m₂² − m₁² = m₂²  →  m₂ = √(Δm²₂₁)
    #   Δm²₃₁ = m₃² − m₁² = m₃²  →  m₃ = √(Δm²₃₁)
    m1 = 0.0
    m2 = math.sqrt(dm2_21)   # eV
    m3 = math.sqrt(dm2_31)   # eV

    m1_meV = m1 * 1e3
    m2_meV = m2 * 1e3
    m3_meV = m3 * 1e3
    sum_m_nu_meV = m1_meV + m2_meV + m3_meV

    # Propagate uncertainties
    m2_err = 0.5 * dm2_21_err / math.sqrt(dm2_21) * 1e3  # meV
    m3_err = 0.5 * dm2_31_err / math.sqrt(dm2_31) * 1e3  # meV
    sum_err = m2_err + m3_err  # conservative (linear)

    print(f"  Structural prediction:  m₁ = 0  (rank-2 seesaw)")
    print(f"")
    print(f"  NuFit 6.0 (Normal Ordering):")
    print(f"    Δm²₂₁ = {dm2_21:.2e} ± {dm2_21_err:.2e} eV²")
    print(f"    Δm²₃₁ = {dm2_31:.4e} ± {dm2_31_err:.4e} eV²")
    print(f"")
    print(f"  m₁ = 0  (predicted)")
    print(f"  m₂ = √(Δm²₂₁) = {m2_meV:.2f} ± {m2_err:.2f} meV")
    print(f"  m₃ = √(Δm²₃₁) = {m3_meV:.2f} ± {m3_err:.2f} meV")
    print(f"")
    print(f"  Σmν = {sum_m_nu_meV:.1f} ± {sum_err:.1f} meV")

    # ── 6b.4  Cosmological and experimental reach ───────────────────

    S("6b.4  Experimental reach")

    print(f"  Oscillation floor (NO, m₁=0): Σmν = {sum_m_nu_meV:.1f} meV")
    print(f"  DESI DR2 sensitivity (ΛCDM): ~58 meV  (current)")
    print(f"  Euclid + DESI projected:     ~40 meV  (target)")
    print(f"  KATRIN endpoint:             < 450 meV (direct)")
    print(f"")
    print(f"  The prediction sits at the oscillation floor:")
    print(f"  any cosmological bound Σmν < 58 meV would rule it out.")
    print(f"  Inverted ordering (m₃ ≈ 0, Σ ≈ 100 meV) is excluded by")
    print(f"  the algebraic structure (rank-2 seesaw forces m₁ = 0).")

    box([
        f"F₄ neutral sector → neutrino masses (m₁ = 0 predicted)",
        f"",
        f"  26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ 1 ⊕ 1",
        f"  2 singlets → 2 right-handed neutrinos",
        f"  Rank-2 Type-I seesaw → m₁ = 0  (normal ordering)",
        f"",
        f"  m₁ = 0           (structural prediction)",
        f"  m₂ = {m2_meV:.2f} meV    (from Δm²₂₁)",
        f"  m₃ = {m3_meV:.2f} meV   (from Δm²₃₁)",
        f"  Σmν = {sum_m_nu_meV:.1f} meV   (testable by DESI/Euclid)",
        f"",
        f"  Falsifiable: any bound Σmν < {sum_m_nu_meV:.0f} meV excludes",
        f"  the framework.  Inverted ordering is excluded.",
    ])

    return {
        'm1_meV': m1_meV,
        'm2_meV': m2_meV,
        'm3_meV': m3_meV,
        'sum_m_nu_meV': sum_m_nu_meV,
        'sum_m_nu_err_meV': sum_err,
        'ordering': 'normal',
        'n_RH': n_singlets,
        'dm2_21': dm2_21,
        'dm2_31': dm2_31,
    }
