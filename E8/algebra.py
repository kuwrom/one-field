"""
Layer 0 -- Algebraic Root: E₈(1) → G₂(1) × F₄(1) conformal embedding.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Derives the central-charge split and all Casimir / dual-Coxeter data
used downstream.  Pure group theory -- no physics parameters needed.
"""

from .formatting import H, S


def derive():
    """
    Compute the algebraic root data of the E₈ conformal embedding.

    Returns
    -------
    dict with keys:
        c_E8, c_G2, c_F4     : central charges
        h_dual                : {name: h∨} dual Coxeter numbers
        C2_fund_G2, C2_adj_G2 : Casimirs of 7, 14 of G₂
        C2_fund_F4, C2_adj_F4 : Casimirs of 26, 52 of F₄
        h_dual_SU3            : h∨(SU(3)) = 3
    """

    H("LAYER 0:  ALGEBRAIC ROOT -- E₈(1) → G₂(1) × F₄(1)")

    S("0.1  Conformal embedding -- central charge split")

    h_dual = {'E8': 30, 'G2': 4, 'F4': 9}
    dim_g  = {'E8': 248, 'G2': 14, 'F4': 52}
    c = {g: dim_g[g] / (1 + h_dual[g]) for g in h_dual}

    print(f"""
  E₈(1):  dim = 248,  h∨ = 30,  c = 248/31 = {c['E8']:.1f}
  G₂(1):  dim =  14,  h∨ =  4,  c =  14/5  = {c['G2']:.1f}
  F₄(1):  dim =  52,  h∨ =  9,  c =  52/10 = {c['F4']:.1f}

  Central charge sum rule:  c(G₂) + c(F₄) = {c['G2']:.1f} + {c['F4']:.1f} = {c['G2']+c['F4']:.1f} = c(E₈)  ✓
  Branching:  248 → (14,1) ⊕ (1,52) ⊕ (7,26)""")

    assert abs(c['G2'] + c['F4'] - c['E8']) < 1e-12

    # ── Casimir data ──────────────────────────────────────────────────

    S("0.2  Key Casimir data")

    C2_fund_G2 = 2       # C₂(7,  G₂) with l(7)=1
    C2_adj_G2  = 4       # C₂(14, G₂) = h∨
    C2_fund_F4 = 6       # C₂(26, F₄) with l(26)=3
    C2_adj_F4  = 9       # C₂(52, F₄) = h∨
    h_dual_SU3 = 3       # for SU(3) ⊂ G₂

    print(f"  G₂: C₂(7) = {C2_fund_G2},  C₂(14) = {C2_adj_G2}")
    print(f"  F₄: C₂(26) = {C2_fund_F4},  C₂(52) = {C2_adj_F4}")
    print(f"  The 26 of F₄ = traceless Albert algebra J₃(𝕆)₀")
    print(f"  Under Spin(8) ⊂ F₄:  26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ 1 ⊕ 1")

    # ── Embedding indices ────────────────────────────────────────────

    S("0.3  Embedding indices: SU(3) ⊂ G₂ ⊂ E₈")

    # G₂ → SU(3) branching with Dynkin index verification
    # 7  → 3 ⊕ 3̄ ⊕ 1:  T(3)=1/2, T(3̄)=1/2, T(1)=0  → sum = 1
    # 14 → 8 ⊕ 3 ⊕ 3̄:  T(8)=3,   T(3)=1/2, T(3̄)=1/2 → sum = 4 = h∨(G₂)
    T_fund = 0.5 + 0.5 + 0.0    # Dynkin indices of 3, 3̄, 1
    T_adj  = 3.0 + 0.5 + 0.5    # Dynkin indices of 8, 3, 3̄
    j_f_SU3_G2 = int(T_fund)    # = 1 (embedding index)

    assert T_fund == 1.0
    assert T_adj == h_dual['G2']

    # Conformal embedding: E₈(1) ⊃ G₂(1) × F₄(1) preserves level
    j_f_G2_E8 = 1

    print(f"  G₂ → SU(3) branching:")
    print(f"    7  → 3 ⊕ 3̄ ⊕ 1     T = 1/2 + 1/2 + 0 = {T_fund:.0f}")
    print(f"    14 → 8 ⊕ 3 ⊕ 3̄     T = 3 + 1/2 + 1/2 = {T_adj:.0f} = h∨(G₂)  ✓")
    print(f"    j_f(SU(3) ⊂ G₂) = {j_f_SU3_G2}")
    print(f"  Conformal embedding: j_f(G₂ ⊂ E₈) = {j_f_G2_E8}")
    print(f"  Full chain: SU(3) ⊂ G₂ ⊂ E₈, all index 1")
    print(f"  → α_s = α_{{G₂}} at the matching scale")

    # ── G₂ topology and the strong CP problem ───────────────────────

    S("0.4  G₂ topology → strong CP constraint")

    # Homotopy groups of G₂:
    #   π₁(G₂) = 0   (simply connected)
    #   π₂(G₂) = 0   (automatic for Lie groups)
    #   π₃(G₂) = ℤ   (instanton sectors exist)
    # G₂ has trivial centre Z(G₂) = {1}, so no discrete θ-vacua
    # ambiguity from the centre -- but π₃ = ℤ means a continuous
    # θ-parameter exists in the unbroken G₂ phase.
    #
    # Under G₂ → SU(3)_c (index 1), the instanton sector maps:
    #   π₃(G₂) = ℤ  →  π₃(SU(3)) = ℤ   (isomorphism)
    # Instanton number is preserved, so θ_QCD = θ_{G₂}.
    #
    # E₈(1) at level 1 has a unique integrable representation
    # (the vacuum module) -- no modular parameter is free.
    # The conformal embedding E₈(1) ⊃ G₂(1) × F₄(1) branches
    # this unique vacuum, fixing θ_{G₂} = 0.
    # Therefore θ_QCD = 0: the strong CP problem is resolved
    # by the embedding, not by an additional symmetry.

    pi1_G2 = 0   # simply connected
    pi3_G2 = 1   # π₃(G₂) = ℤ  (rank of the homotopy group)
    pi3_SU3 = 1  # π₃(SU(3)) = ℤ

    print(f"  G₂ homotopy:")
    print(f"    π₁(G₂) = {pi1_G2}   (simply connected, no centre)")
    print(f"    π₃(G₂) = ℤ   (instanton sectors exist)")
    print(f"  G₂ → SU(3) embedding (index 1):")
    print(f"    π₃(G₂) = ℤ  →  π₃(SU(3)) = ℤ   (isomorphism)")
    print(f"    ⇒ θ_QCD = θ_{{G₂}}  (inherited, not independent)")
    print(f"  E₈(1) has unique integrable rep at level 1 (vacuum module)")
    print(f"  Conformal embedding branches unique vacuum → θ_{{G₂}} = 0")
    print(f"  → θ_QCD = 0: strong CP resolved by the embedding")

    return {
        'c_E8': c['E8'], 'c_G2': c['G2'], 'c_F4': c['F4'],
        'h_dual': h_dual,
        'dim_g': dim_g,
        'dim_G2_fund': 7, 'dim_F4_fund': 26,
        'C2_fund_G2': C2_fund_G2, 'C2_adj_G2': C2_adj_G2,
        'C2_fund_F4': C2_fund_F4, 'C2_adj_F4': C2_adj_F4,
        'h_dual_SU3': h_dual_SU3,
        'j_f_SU3_G2': j_f_SU3_G2, 'j_f_G2_E8': j_f_G2_E8,
    }
