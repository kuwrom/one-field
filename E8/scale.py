"""
Layer 1 -- Scale Emergence: v_EW from M_Pl.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The electroweak scale is an instanton suppression factor:
    v_EW = M_Pl × exp(−(9π²/2 − 6 + 15/512))

where 9π²/2 is the lepton instanton action (h∨(F₄)·π²/2),
6 = C₂(26,F₄) is the Casimir correction, and 15/512 is the
one-loop QED vertex correction from 30 modes × α/(2π):

    30 = 26 (full F₄ fundamental under Spin(9))
       + 4  (real Higgs-doublet scalar DOFs)

with α = π/512 derived from the Singh ratio (α_s/α_em = 16).

Physical motivation (readable without the papers)
──────────────────────────────────────────────────

WHY the electroweak scale is an instanton suppression of M_Pl:

  The Z₃-NLS substrate has a natural UV scale: the Planck mass M_Pl,
  set by the healing length ξ₀ = ℓ_Pl/2.  Standing-wave knots
  (particles) form at scales exponentially below M_Pl, because
  their stability requires tunnelling through a barrier whose
  height is set by the instanton action of the gauge sector.

  In the lepton branch (G₂ sector), the instanton action is
  S₀ = 2π/(b₀·α_{G₂}) where b₀ = 32/3 is the one-loop beta
  coefficient and α_{G₂}(M_Pl) = 1/(24π).  This gives
  S₀ = 2π × 24π/(32/3) = 9π²/2 ≈ 44.41.  The lepton confinement
  scale is Λ_conf = ½ M_Pl exp(−S₀) ≈ 314 MeV.

WHY the Casimir correction −C₂(26) = −6 appears:

  Quarks live in the 26 of F₄ (one emergence layer deeper than
  leptons, which sit in the identity representation).  In the
  instanton background, the functional determinant for a field
  in representation λ contributes exp(+C₂(λ)) to the effective
  action.  For the 26 of F₄, C₂ = 6 (from the Dynkin index
  ℓ(26)·dim(F₄)/dim(26) = 3×52/26 = 6).  This REDUCES the
  effective action: S_eff = S₀ − C₂(26) = 9π²/2 − 6, pushing
  the quark emergence scale UP relative to leptons.

WHY the 30-mode vertex correction 15/512 appears:

  The 't Hooft instanton vertex receives a one-loop contribution
  α/(2π) from each charged mode in the background.  In the F₄
  emergence layer:
    • 26 modes from the F₄ fundamental (all 26 contribute because
      F₄ treats the 26 irreducibly, including the Spin(9) singlet)
    • 4 real Higgs-doublet scalar DOFs (enter through the electroweak
      instanton sector governed by the F₄(1) fusion rule)
  Total: N = 30.  With α = π/512 (from the Singh ratio):
  δS = 30 × (π/512)/(2π) = 15/512 ≈ 0.029.

  The identity channel (C₂=0) gives Λ_conf.
  The 26 channel (C₂=6, + vertex correction) gives v_EW.
  The ratio v_EW/M_Pl = exp(−38.44) is the exponential of an
  algebraic constant, the hierarchy problem is resolved.
"""

import math
from .constants import (
    M_Pl_GeV, M_Pl_MeV, ALPHA_EM,
    N_VERTEX, N_F4_FUND, N_HIGGS,
)
from .formatting import H, S, box


def derive(alg: dict):
    """
    Derive the electroweak scale from the Planck mass.

    Parameters
    ----------
    alg : dict from algebra.derive()

    Returns
    -------
    dict with keys:
        S_0               : lepton instanton action
        S_quark_LO        : quark instanton action (Casimir only)
        S_quark            : full quark instanton action (+ 30-mode vertex)
        v_EW_pred_GeV      : predicted v_EW in GeV
        v_EW_pred_MeV      : predicted v_EW in MeV
        Lambda_conf_MeV    : lepton confinement scale in MeV
    """

    C2_fund_F4 = alg['C2_fund_F4']

    H("LAYER 1:  SCALE EMERGENCE -- v_EW FROM M_Pl")

    # ── 1.1  Instanton action ─────────────────────────────────────────

    S("1.1  Instanton action")

    # Lepton instanton action: S₀ = h∨(F₄)·π²/2 = 9π²/2
    # h∨(F₄) = 9 is the dual Coxeter number of F₄.
    # This is one-loop dimensional transmutation: S₀ = 2π/(b₀·α_{G₂}),
    # with b₀(G₂) = 32/3 and α_{G₂}(M_Pl) = 1/(24π).
    S_0 = 9.0 * math.pi**2 / 2.0          # h∨(F₄)·π²/2

    # Quark instanton action (leading order): subtract the Casimir of
    # the 26 of F₄.  The functional determinant in the instanton
    # background contributes exp(+C₂(λ)), reducing the effective action.
    S_quark_LO = S_0 - C2_fund_F4          # Casimir correction

    print(f"  S₀ (lepton) = h∨(F₄) · π²/2 = 9 · π²/2 = {S_0:.4f}")
    print(f"  S_quark(LO) = S₀ − C₂(26)   = {S_0:.4f} − {C2_fund_F4} = {S_quark_LO:.4f}")

    # ── 1.2  One-loop QED vertex correction (30 modes) ────────────────

    S("1.2  Instanton-vertex term (30 modes)")

    # 't Hooft instanton vertex correction: each charged mode in the
    # background contributes α/(2π) to the instanton effective action.
    #   26 modes: full F₄ fundamental (including Spin(9) singlet, because
    #             F₄ treats the 26 irreducibly, see docstring)
    #    4 modes: real Higgs doublet DOFs (enter through electroweak
    #             instanton sector governed by F₄(1) fusion rule)
    #   Total: 30 modes.  α = π/512 (from Singh ratio, see constants.py).
    N_vertex = N_VERTEX                  # = 30, derived from branching rules
    delta_S_QED = N_vertex * ALPHA_EM / (2.0 * math.pi)
    S_quark = S_quark_LO + delta_S_QED

    print(f"  26 of F₄ under Spin(9) = 1 ⊕ 9 ⊕ 16  (all {N_F4_FUND} modes)")
    print(f"  Higgs doublet: {N_HIGGS} real scalar DOFs")
    print(f"  Total vertex modes: {N_F4_FUND} + {N_HIGGS} = {N_vertex}")
    print(f"  α_em = π/512 = {ALPHA_EM:.8f}  (derived: Singh ratio)")
    print(f"  δS = {N_vertex} × α/(2π) = {N_vertex} × (π/512)/(2π) = 15/512 = {delta_S_QED:.10f}")
    print(f"  S_quark = {S_quark_LO:.4f} + {delta_S_QED:.10f} = {S_quark:.10f}")

    # Verification: fitted mode count
    v_EW_PDG_MeV = 246.22e3
    S_exact = math.log(M_Pl_MeV / v_EW_PDG_MeV)
    N_fitted = (S_exact - S_quark_LO) / (ALPHA_EM / (2.0 * math.pi))
    print(f"  Fitted mode count to hit PDG v_EW exactly = {N_fitted:.2f}")
    print(f"  Derived count = {N_vertex}  (deviation: {N_fitted - N_vertex:.2f} modes)")

    # ── 1.3  Lepton confinement scale ────────────────────────────────

    S("1.3  Lepton confinement scale")

    # Λ_conf = ½ M_Pl exp(−S₀): the scale where lepton standing-wave
    # knots first become stable.  The factor ½ is the substrate healing
    # length ξ₀ = ℓ_Pl/2.  This is the identity-channel (C₂=0) scale.
    Lambda_conf_MeV = 0.5 * M_Pl_MeV * math.exp(-S_0)
    print(f"  Λ_conf = ½ · M_Pl · exp(−S₀) = {Lambda_conf_MeV:.2f} MeV")

    # ── 1.4  Electroweak scale ────────────────────────────────────────

    S("1.4  Electroweak scale")

    v_EW_LO_MeV = M_Pl_MeV * math.exp(-S_quark_LO)
    v_EW_pred_MeV = M_Pl_MeV * math.exp(-S_quark)
    v_EW_pred_GeV = v_EW_pred_MeV / 1e3
    v_EW_PDG = 246.22

    err_LO = 100.0 * (v_EW_LO_MeV / 1e3 - v_EW_PDG) / v_EW_PDG
    err_NLO = 100.0 * (v_EW_pred_GeV - v_EW_PDG) / v_EW_PDG

    box([
        f"Casimir only:  v_EW = M_Pl × exp(−(9π²/2 − 6))",
        f"                    = {v_EW_LO_MeV/1e3:.2f} GeV   (PDG: {v_EW_PDG} GeV,  {err_LO:+.1f}%)",
        f"",
        f"Full:  v_EW = M_Pl × exp(−(9π²/2 − 6 + 15/512))",
        f"            = {v_EW_pred_GeV:.2f} GeV   (PDG: {v_EW_PDG} GeV,  {err_NLO:+.2f}%)",
        f"",
        f"Vertex term: 30 × (π/512)/(2π) = 15/512 = {delta_S_QED:.10f}",
        f"  26 (F₄ fund) + 4 (Higgs doublet) = 30 modes",
        f"Improvement: {err_LO:+.1f}% → {err_NLO:+.2f}%",
    ])

    return {
        'S_0': S_0,
        'S_quark_LO': S_quark_LO,
        'S_quark': S_quark,
        'v_EW_pred_GeV': v_EW_pred_GeV,
        'v_EW_pred_MeV': v_EW_pred_MeV,
        'Lambda_conf_MeV': Lambda_conf_MeV,
    }
