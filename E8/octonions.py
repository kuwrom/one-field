"""
Layer 4b -- Octonionic Clebsch-Gordan Verification.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Verifies the octonionic G₂/SU(3) Clebsch-Gordan coefficients that
fix the Brannen amplitude ratio B/A = √2 and hence the Koide
relation Q = 2/3.

The calculation follows Supporting Math Sections 5-6:
  1. Build the octonion multiplication table on Im(𝕆) ≅ ℝ⁷
  2. Define the holomorphic basis z_α = (e_{2α-1} - i e_{2α})/√2
  3. Compute the singlet channel:   z_α × z̄_β = i δ_{αβ} e₇  →  |C₁| = 1
  4. Compute the antisymmetric channel: z_α × z_β = √2 ε_{αβγ} z̄_γ  →  |C₃̄| = √2
  5. Channel weight: W(3̄) = 1/d(3) = 1/2
  6. Amplitude ratio: B/A = 2|C₃̄|/(|C₁|·d(3)) = 2√2/2 = √2

No free parameters.  The ratio |C₃̄|/|C₁| = √2 is frame-independent
(both maps are unique up to phase).
"""

import math
import numpy as np
from .formatting import H, S


# ═══════════════════════════════════════════════════════════════════════
#  Octonion multiplication table
# ═══════════════════════════════════════════════════════════════════════

# Fano-plane triples (i,j,k) with e_i × e_j = e_k (cyclic).
# Convention: e₇ is the SU(3)-invariant direction.
# The six-plane (e₁,...,e₆) carries the 3 ⊕ 3̄ of SU(3).
FANO_TRIPLES = [
    (1, 2, 7),
    (3, 4, 7),
    (5, 6, 7),
    (1, 3, 5),
    (1, 4, 6),  # note: e₁ × e₄ = -e₆ (orientation)
    (2, 3, 6),  # note: e₂ × e₃ = -e₆ (orientation)
    (2, 4, 5),  # note: e₂ × e₄ = -e₅ (orientation)
]

# Signs: for triple (i,j,k), e_i × e_j = +e_k or -e_k.
# The sign is +1 for cyclic order on the Fano plane, -1 for anti-cyclic.
FANO_SIGNS = [
    +1,  # e₁ × e₂ = +e₇
    +1,  # e₃ × e₄ = +e₇
    +1,  # e₅ × e₆ = +e₇
    +1,  # e₁ × e₃ = +e₅
    -1,  # e₁ × e₄ = -e₆
    -1,  # e₂ × e₃ = -e₆
    -1,  # e₂ × e₄ = -e₅
]


def _cross_complex(a, b):
    """
    Cross product on Im(𝕆) ≅ ℝ⁷ extended to complex coefficients.

    For each Fano triple (i,j,k) with sign s:
      e_i × e_j = s * e_k,  e_j × e_k = s * e_i,  e_k × e_i = s * e_j
    (cyclic with same sign; antisymmetric in the two arguments)
    """
    result = np.zeros(7, dtype=complex)
    for (i, j, k), sign in zip(FANO_TRIPLES, FANO_SIGNS):
        ii, jj, kk = i - 1, j - 1, k - 1
        result[kk] += sign * (a[ii] * b[jj] - a[jj] * b[ii])
        result[ii] += sign * (a[jj] * b[kk] - a[kk] * b[jj])
        result[jj] += sign * (a[kk] * b[ii] - a[ii] * b[kk])
    return result


def derive(wzw_data: dict):
    """
    Verify the octonionic CG coefficients and derive B/A = √2.

    Parameters
    ----------
    wzw_data : dict from wzw.derive(), providing d(1,0) = 2

    Returns
    -------
    dict with keys:
        C1, C3bar     : |C₁| and |C₃̄| CG coefficients
        BA_ratio      : derived B/A = √2
        Q0            : derived Q₀ = 2/3
    """

    d_fund = wzw_data['d10']  # quantum dimension d(1,0) = 2

    H("LAYER 4b:  OCTONIONIC CLEBSCH-GORDAN VERIFICATION")

    # ── 4b.1  Holomorphic basis ──────────────────────────────────────

    S("4b.1  Holomorphic basis on the SU(3) six-plane")

    # z_α = (e_{2α-1} - i e_{2α}) / √2,  α = 1,2,3
    # z̄_α = (e_{2α-1} + i e_{2α}) / √2
    # e₇ is the SU(3)-invariant direction.

    z = np.zeros((3, 7), dtype=complex)
    zbar = np.zeros((3, 7), dtype=complex)
    for alpha in range(3):
        z[alpha, 2 * alpha] = 1.0 / math.sqrt(2)       # e_{2α-1} coefficient
        z[alpha, 2 * alpha + 1] = -1j / math.sqrt(2)    # -i e_{2α} coefficient
        zbar[alpha] = z[alpha].conj()

    e7 = np.zeros(7, dtype=complex)
    e7[6] = 1.0

    print(f"  z_α = (e_{{2α-1}} − i e_{{2α}}) / √2,  α = 1,2,3")
    print(f"  e₇ = SU(3)-invariant unit imaginary octonion")

    # ── 4b.2  Singlet channel: 3 ⊗ 3̄ → 1 ────────────────────────────

    S("4b.2  Singlet channel: 3 ⊗ 3̄ → 1")

    print(f"  Computing z_α × z̄_β for all α,β:")
    C1_values = []
    for alpha in range(3):
        for beta in range(3):
            prod = _cross_complex(z[alpha], zbar[beta])
            # Should be i δ_{αβ} e₇
            e7_coeff = prod[6]
            off_e7 = np.max(np.abs(prod[:6]))
            if alpha == beta:
                print(f"    z_{alpha+1} × z̄_{beta+1} = ({e7_coeff:.4f}) e₇"
                      f"   [expected: i e₇ = {1j}]   off-e₇ residual: {off_e7:.2e}")
                C1_values.append(abs(e7_coeff))
            else:
                print(f"    z_{alpha+1} × z̄_{beta+1} = ({e7_coeff:.4f}) e₇"
                      f"   [expected: 0]   off-e₇ residual: {off_e7:.2e}")

    C1 = np.mean(C1_values)
    assert abs(C1 - 1.0) < 1e-10, f"|C₁| = {C1}, expected 1"
    print(f"\n  |C₁| = {C1:.6f}  ✓")

    # ── 4b.3  Antisymmetric channel: 3 ⊗ 3 → 3̄ ─────────────────────

    S("4b.3  Antisymmetric channel: 3 ⊗ 3 → 3̄")

    print(f"  Computing z_α × z_β for α ≠ β:")
    C3bar_values = []

    # z₁ × z₂ should give √2 z̄₃
    pairs = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]  # (α, β, γ) with ε_{αβγ} = +1
    for alpha, beta, gamma in pairs:
        prod = _cross_complex(z[alpha], z[beta])
        # Project onto z̄_γ: overlap = <prod, z̄_γ> using standard inner product
        overlap = np.vdot(zbar[gamma], prod)  # conjugate-linear in first arg
        norm_zbar = np.sqrt(np.vdot(zbar[gamma], zbar[gamma]).real)
        coeff_abs = abs(overlap) / norm_zbar
        print(f"    z_{alpha+1} × z_{beta+1}:  coefficient of z̄_{gamma+1}"
              f" = {overlap/norm_zbar:.6f}   |coeff| = {coeff_abs:.6f}"
              f"   [expected: √2 = {math.sqrt(2):.6f}]")
        C3bar_values.append(coeff_abs)

    C3bar = np.mean(C3bar_values)
    assert abs(C3bar - math.sqrt(2)) < 1e-10, f"|C₃̄| = {C3bar}, expected √2"
    print(f"\n  |C₃̄| = {C3bar:.6f} = √2  ✓")

    # ── 4b.4  Frame independence ─────────────────────────────────────

    S("4b.4  Frame independence")

    ratio_CG = C3bar / C1
    print(f"  |C₃̄| / |C₁| = {ratio_CG:.6f} = √2")
    print(f"  Both maps (3⊗3̄→1 and 3⊗3→3̄) are unique up to phase,")
    print(f"  so the absolute ratio √2 is frame-independent.")

    # ── 4b.5  Channel weight and B/A derivation ──────────────────────

    S("4b.5  Channel weight and amplitude ratio")

    W_fund = 1.0 / d_fund  # = 1/2
    print(f"  Quantum dimension d(3) = {d_fund}")
    print(f"  Channel weight W(3̄) = 1/d(3) = {W_fund:.6f}")

    # B/A = 2|C₃̄| · W(3̄) / |C₁|
    #     = 2 × √2 × (1/2) / 1
    #     = √2
    BA_derived = 2.0 * C3bar * W_fund / C1
    print(f"\n  B/A = 2|C₃̄| · W(3̄) / |C₁|")
    print(f"      = 2 × {C3bar:.4f} × {W_fund:.4f} / {C1:.4f}")
    print(f"      = {BA_derived:.6f}")
    assert abs(BA_derived - math.sqrt(2)) < 1e-10
    print(f"      = √2  ✓")

    # ── 4b.6  Koide parameter ────────────────────────────────────────

    S("4b.6  Koide parameter Q₀")

    Q0 = 1.0 / 3.0 + BA_derived**2 / 6.0
    print(f"  Q = 1/3 + |B/A|²/6 = 1/3 + 2/6 = {Q0:.10f}")
    assert abs(Q0 - 2.0 / 3.0) < 1e-10
    print(f"  Q₀ = 2/3  ✓  (Koide relation)")

    return {
        'C1': C1,
        'C3bar': C3bar,
        'BA_ratio': BA_derived,
        'Q0': Q0,
    }
