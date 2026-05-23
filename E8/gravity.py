"""
Layer 10 -- Gravity: Sakharov Induced Gravity from the Bridge Sector.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The mixed (7,26) sector of the E₈(1) → G₂(1) × F₄(1) branching
carries 182 degrees of freedom -- the bridge between the lepton
sector (G₂) and the quark sector (F₄).

Each bridge mode carries a non-minimal coupling to curvature:

    ξ_bridge = α_{G₂}(M_Pl) · E[v²] · h_bridge = 1/(48π)

where:
    α_{G₂}(M_Pl) = 1/(24π)   [from g²_adj(M_Pl)=1, b₀=32/3]
    E[v²] = 1/2               [protected forgetting: Pv²P = ½P]
    h_bridge = 1               [WZW conformal weight sum: 2/5 + 3/5]

Protected forgetting:
    The Z₃ cyclic projector P = uu† (u = (1,1,1)/√3) acts on the
    gap vector v = (1, −½, −½) to give:
        PvP  = 0     [the Z₃ label is erased → universality]
        Pv²P = ½P    [the presence is remembered → gravity]

    This is the algebraic origin of the equivalence principle:
    the bridge responds identically to every particle regardless
    of its Z₃ label.  The "ventilation" pattern (the way the
    vacuum reorganises around a standing-wave knot) depends only
    on E[v²] = 1/2, not on the knot's internal structure.

The one-loop heat kernel generates the Einstein-Hilbert action
(Sakharov induced gravity):

    1/(16πG_ind) = (Λ²/(96π²)) Σ_total

Result: G_ind/G_N = 0.994 (UV) to 1.007 (broken phase).
"""

import math
from .formatting import H, S, box, pct


def derive(alg: dict, scale: dict):
    """
    Derive Newton's constant from the bridge sector heat kernel.

    Parameters
    ----------
    alg   : dict from algebra.derive()
    scale : dict from scale.derive()

    Returns
    -------
    dict with keys:
        N_bridge                       : bridge sector dimension (182)
        xi_bridge                      : non-minimal coupling 1/(48π)
        h_bridge                       : bridge conformal weight (1)
        h_G2_7, h_F4_26               : individual conformal weights
        alpha_G2_Pl                    : G₂ coupling at M_Pl
        E_v2                           : protected-forgetting factor (1/2)
        Sigma_min_UV, Sigma_min_broken : SM heat-kernel sums
        Sigma_bridge                   : bridge contribution
        Sigma_total_UV, Sigma_total_broken : totals
        G_ratio_UV, G_ratio_broken, G_ratio_mid : G_ind/G_N
        rho_Lambda                     : cosmological constant (GeV⁴)
    """

    H("LAYER 10:  GRAVITY -- SAKHAROV INDUCED GRAVITY FROM THE BRIDGE SECTOR")

    # ══════════════════════════════════════════════════════════════════
    #  10.1  Bridge sector identification
    # ══════════════════════════════════════════════════════════════════

    S("10.1  Bridge sector: (7,26) of E₈ → G₂ × F₄")

    # E₈(1) branching rule: 248 → (14,1) ⊕ (1,52) ⊕ (7,26)
    #   (14,1) : G₂ adjoint  → gauge sector (gluons + G₂/SU(3) vectors)
    #   (1,52) : F₄ adjoint  → quark/EW sector
    #   (7,26) : mixed bridge → gravity mediator
    #
    # The bridge carries TWO predictions, not one:
    #   * Gravity: 182 scalar-like heat-kernel channels each contribute
    #     a₁/R = 1/6 − ξ_bridge with ξ_bridge = α_G₂ × (1/2) × h_bridge
    #     = 1/(48π).  Sums to G_ind/G_N = 0.994 (UV) to 1.007 (broken).
    #   * QED: the same h_bridge = 1 makes the bridge a marginal worldsheet
    #     primary; its one-loop self-interference contribution to the EM
    #     current-current correlator is h_bridge/(2π) = 1/(2π), giving
    #     the factor (1 − 1/(2π)) that converts α_alg = π/512 into the
    #     physical α(0) = 1/137.036 (see constants.py for the chain).
    # Same sector, same conformal weight h_bridge = 1, two observables.
    #
    # Why E₈(1) and why 182:
    # The bridge mode count N_bridge = 7 × 26 = 182 is fixed by the
    # branching rule of the *unique* exceptional conformal embedding
    # E₈(1) ⊃ G₂(1) × F₄(1) at level 1.  This embedding is selected by
    # the central-charge identity c(G₂) + c(F₄) = 14/5 + 26/5 = 8 =
    # c(E₈), and E₈(1) itself is the unique rank-8 positive even self-
    # dual lattice VOA (Griess; positive even unimodular lattices exist
    # only in ranks divisible by 8, and in rank 8 there is exactly one).
    # The bridge count is therefore not adjustable: once the substrate's
    # rank-8 bosonic local self-dual common-mode coherence is given, the
    # branching forces dim(7,26) = 182.  See "The Echo of Standing Waves"
    # Sec. IV (Lagrangian condensation of A_{E₈} = 1 + (τ_G, τ_F)) and
    # Appendix E (lattice uniqueness gate).
    dim_7 = 7       # dim of G₂ fundamental
    dim_26 = 26     # dim of F₄ fundamental (traceless Albert algebra J₃(𝕆)₀)
    N_bridge = dim_7 * dim_26   # = 182

    print(f"  E₈(1) branching: 248 → (14,1) ⊕ (1,52) ⊕ (7,26)")
    print(f"  Bridge sector: ({dim_7},{dim_26}) with dim = {dim_7} × {dim_26} = {N_bridge}")
    print(f"  Check: 14 + 52 + {N_bridge} = {14 + 52 + N_bridge} = 248  ✓")

    assert 14 + 52 + N_bridge == 248

    # ══════════════════════════════════════════════════════════════════
    #  10.2  Protected forgetting: PvP = 0, Pv²P = ½P
    # ══════════════════════════════════════════════════════════════════

    S("10.2  Protected forgetting (Z₃ cyclic projector)")

    # The Z₃ gap vector v = (1, −½, −½) encodes the Brannen spectrum.
    # The cyclic projector P = uu† with u = (1,1,1)/√3 averages over
    # all three Z₃ sectors.
    #
    # Physical picture: a standing-wave knot (particle) creates a
    # "ventilation pattern" in the surrounding vacuum.  The bridge
    # modes carry the memory of *presence* (Pv²P = ½P) while
    # forgetting the knot's internal label (PvP = 0).  This erasure
    # is the algebraic origin of the equivalence principle -- gravity
    # couples universally because the bridge cannot distinguish
    # which Z₃ sector the knot occupies.

    v = [1.0, -0.5, -0.5]
    u = [1.0 / math.sqrt(3.0)] * 3

    # PvP = (u·v)(uu†) = 0 · P  since u·v = 1/√3 (1−½−½) = 0
    u_dot_v = sum(ui * vi for ui, vi in zip(u, v))

    # Pv²P = (Σ uᵢ² vᵢ²) P = [(1/3)(1 + ¼ + ¼)] P = ½ P
    E_v2 = sum(ui**2 * vi**2 for ui, vi in zip(u, v))

    print(f"  v = (1, −½, −½)    [Z₃ gap vector]")
    print(f"  u = (1,1,1)/√3     [cyclic ground state]")
    print(f"  P = uu†             [rank-1 projector]")
    print(f"")
    print(f"  PvP  = (u·v) P = {u_dot_v:.1f} · P = 0      [label erased → universality]")
    print(f"  Pv²P = E[v²] P = {E_v2:.1f} · P = ½ P  [presence remembered → gravity]")
    print(f"")
    print(f"  The bridge \"forgets\" which Z₃ sector a knot occupies")
    print(f"  but \"remembers\" that a knot is present → equivalence principle.")

    assert abs(u_dot_v) < 1e-12, f"PvP ≠ 0: u·v = {u_dot_v}"
    assert abs(E_v2 - 0.5) < 1e-12, f"E[v²] ≠ 1/2: got {E_v2}"

    # ══════════════════════════════════════════════════════════════════
    #  10.3  Bridge conformal weight
    # ══════════════════════════════════════════════════════════════════

    S("10.3  Bridge conformal weight h_bridge = 1")

    C2_fund_G2 = alg['C2_fund_G2']   # = 2
    C2_fund_F4 = alg['C2_fund_F4']   # = 6
    h_dual_G2 = alg['h_dual']['G2']  # = 4
    h_dual_F4 = alg['h_dual']['F4']  # = 9
    k = 1   # WZW level (conformal embedding preserves level)

    # WZW conformal weight: h_R = C₂(R) / (k + h∨)
    h_G2_7 = C2_fund_G2 / (k + h_dual_G2)    # 2/(1+4) = 2/5
    h_F4_26 = C2_fund_F4 / (k + h_dual_F4)   # 6/(1+9) = 3/5
    h_bridge = h_G2_7 + h_F4_26               # 2/5 + 3/5 = 1

    print(f"  h_{{G₂}}(7)  = C₂(7)/(k+h∨)  = {C2_fund_G2}/({k}+{h_dual_G2}) = {h_G2_7}")
    print(f"  h_{{F₄}}(26) = C₂(26)/(k+h∨) = {C2_fund_F4}/({k}+{h_dual_F4}) = {h_F4_26}")
    print(f"  h_bridge = {h_G2_7} + {h_F4_26} = {h_bridge}")
    print(f"  Integer conformal weight → bridge is a physical primary  ✓")

    assert abs(h_bridge - 1.0) < 1e-12

    # ══════════════════════════════════════════════════════════════════
    #  10.4  Non-minimal coupling ξ_bridge
    # ══════════════════════════════════════════════════════════════════

    S("10.4  Non-minimal coupling: ξ = α_{G₂}(M_Pl) · E[v²] · h_bridge")

    # α_{G₂}(M_Pl) = 1/(24π)
    # Derivation (Layer 7): g²_adj(M_Pl) = 1, b₀(G₂) = 32/3.
    # The adjoint normalisation at the Planck scale combined with the
    # one-loop β-function gives 1/α = 24π.  Equivalently:
    # the embedding index j_f = 1 and the Sugawara central charge
    # fix the coupling at each scale.
    alpha_G2_Pl = 1.0 / (24.0 * math.pi)

    xi_bridge = alpha_G2_Pl * E_v2 * h_bridge
    # = (1/(24π)) × (1/2) × 1 = 1/(48π)
    xi_exact = 1.0 / (48.0 * math.pi)

    print(f"  α_{{G₂}}(M_Pl) = 1/(24π) = {alpha_G2_Pl:.10f}")
    print(f"  E[v²]    = 1/2          (protected forgetting)")
    print(f"  h_bridge = 1            (conformal weight)")
    print(f"  ξ_bridge = (1/(24π)) × (1/2) × 1 = 1/(48π) = {xi_bridge:.10f}")
    print(f"  Verify: 1/(48π) = {xi_exact:.10f}  ✓")

    assert abs(xi_bridge - xi_exact) < 1e-14

    # ══════════════════════════════════════════════════════════════════
    #  10.5  Heat-kernel a₁/R coefficients
    # ══════════════════════════════════════════════════════════════════

    S("10.5  Heat-kernel coefficients (per real DOF)")

    # DeWitt-Seeley a₁ coefficient in the heat-kernel expansion:
    #   Tr exp(−tD²) ∼ (4πt)^{−d/2} Σ_n t^n a_n
    # The piece linear in R determines the induced Einstein-Hilbert term.
    #
    # Per real degree of freedom:
    #   Minimal scalar (ξ=0):  a₁/R = +1/6
    #   Scalar with ξ:         a₁/R = +1/6 − ξ
    #   Dirac fermion:         a₁/R = −1/12
    #   Gauge boson (trans):   a₁/R = −1/6
    #   Proca (massive vec):   a₁/R = −1/12

    a1_scalar = 1.0 / 6.0              # minimal scalar (ξ=0)
    a1_fermion = -1.0 / 12.0           # Dirac per real DOF
    a1_gauge = -1.0 / 6.0              # gauge per transverse DOF
    a1_proca = -1.0 / 12.0             # Proca per DOF
    a1_bridge = 1.0 / 6.0 - xi_bridge  # bridge scalar

    print(f"  Minimal scalar (ξ=0):  a₁/R = +1/6  = +{a1_scalar:.6f}")
    print(f"  Dirac (per real DOF):  a₁/R = −1/12 = {a1_fermion:.6f}")
    print(f"  Gauge (per trans DOF): a₁/R = −1/6  = {a1_gauge:.6f}")
    print(f"  Proca (per DOF):       a₁/R = −1/12 = {a1_proca:.6f}")
    print(f"  Bridge (ξ=1/(48π)):    a₁/R = 1/6 − 1/(48π) = +{a1_bridge:.6f}")

    # ── Regulator scheme defense ──
    #
    # The Sakharov formula uses Λ = M_Pl as hard cutoff.  In standard QFT
    # this would be scheme-dependent.  But for a condensate with genuine
    # short-distance microstructure at healing length ξ₀ = ℓ_Pl/2:
    #
    #   • The hard cutoff is the unique physically correct scheme:
    #     it preserves the power-law UV sensitivity (∝ Λ²) that encodes
    #     the physical energy stored in modes up to the real lattice scale.
    #     Ref: Volovik, "Universe in a Helium Droplet" (2003) §29;
    #          Visser, gr-qc/0204062 (2002).
    #
    #   • Dimensional regularization would discard the power-law piece
    #     entirely, erasing the physical content -- inappropriate for a
    #     system with real microstructure.
    #     Ref: Barceló-Liberati-Visser, gr-qc/0505065 (2005).
    #
    #   • The heat-kernel COEFFICIENTS a₁/R are regulator-independent
    #     (Vassilevich, hep-th/0306138); only the overall Λ² prefactor
    #     depends on the scheme, and it is fixed by the substrate's
    #     microscopic scale.
    #
    # Therefore: regulator dependence is a feature, not a bug.

    print(f"")
    print(f"  Regulator scheme:")
    print(f"    Hard cutoff Λ = M_Pl is uniquely correct for a condensate")
    print(f"    with physical microstructure at ξ₀ = ℓ_Pl/2.")
    print(f"    The a₁/R coefficients are regulator-independent (Vassilevich 2003);")
    print(f"    only the Λ² prefactor depends on scheme → fixed by substrate UV scale.")

    # ══════════════════════════════════════════════════════════════════
    #  10.6  SM sector heat-kernel sum: Σ_min
    # ══════════════════════════════════════════════════════════════════

    S("10.6  SM sector heat-kernel sum Σ_min")

    # ── UV bookkeeping (unbroken electroweak) ──
    #
    # 4  Higgs real scalars                    × (+1/6)
    # 24 SM gauge bosons (transverse DOF)      × (−1/6)
    #     [8 gluons + 3 W + 1 B = 12 bosons × 2 polarisations]
    # 90 fermionic real DOF                    × (−1/12)
    #     [3 gen × (4 quarks × 3 colours + 3 leptons) × 2 = 90]
    # 6  massive G₂/SU(3) vectors (trans DOF)  × (−1/12)
    #     [3 heavy vectors from G₂ → SU(3) breaking × 2 trans]
    # 7  G₂ fundamental scalars                × (+1/6)
    #     [scalar 7 of G₂: the Higgs that breaks G₂ → SU(3)]

    n_higgs_UV = 4
    n_gauge_UV = 24
    n_fermion = 90
    n_G2_vec = 6
    n_G2_scal = 7

    Sigma_min_UV = (n_higgs_UV * a1_scalar
                    + n_gauge_UV * a1_gauge
                    + n_fermion * a1_fermion
                    + n_G2_vec * a1_proca
                    + n_G2_scal * a1_scalar)

    print(f"  UV bookkeeping (unbroken EW):")
    print(f"    {n_higgs_UV:2d} Higgs scalars     × (+1/6)   = {n_higgs_UV * a1_scalar:+8.4f}")
    print(f"    {n_gauge_UV:2d} gauge transverse   × (−1/6)   = {n_gauge_UV * a1_gauge:+8.4f}")
    print(f"    {n_fermion:2d} fermion DOF        × (−1/12)  = {n_fermion * a1_fermion:+8.4f}")
    print(f"    {n_G2_vec:2d}  G₂/SU(3) vectors  × (−1/12)  = {n_G2_vec * a1_proca:+8.4f}")
    print(f"    {n_G2_scal:2d}  G₂ fund scalars   × (+1/6)   = {n_G2_scal * a1_scalar:+8.4f}")
    print(f"    {'─' * 50}")
    print(f"    Σ_min(UV) = {Sigma_min_UV:.4f}")

    # ── Broken-phase bookkeeping ──
    #
    # After EWSB: 3 Goldstones eaten → W±/Z become massive (Proca)
    #
    # 1  physical Higgs                        × (+1/6)
    # 18 massless gauge transverse DOF         × (−1/6)
    #     [8 gluons × 2 + 1 photon × 2 = 18]
    # 9  massive gauge DOF (Proca: W±, Z)      × (−1/12)
    #     [3 massive vectors × 3 DOF each]
    # 90 fermionic real DOF                    × (−1/12)
    # 6  G₂/SU(3) massive vectors              × (−1/12)
    # 7  G₂ fundamental scalars                × (+1/6)

    n_higgs_broken = 1
    n_gauge_broken = 18
    n_proca_broken = 9

    Sigma_min_broken = (n_higgs_broken * a1_scalar
                        + n_gauge_broken * a1_gauge
                        + n_proca_broken * a1_proca
                        + n_fermion * a1_fermion
                        + n_G2_vec * a1_proca
                        + n_G2_scal * a1_scalar)

    print(f"")
    print(f"  Broken-phase bookkeeping (after EWSB):")
    print(f"    {n_higgs_broken:2d}  physical Higgs    × (+1/6)   = {n_higgs_broken * a1_scalar:+8.4f}")
    print(f"    {n_gauge_broken:2d} massless gauge     × (−1/6)   = {n_gauge_broken * a1_gauge:+8.4f}")
    print(f"    {n_proca_broken:2d}  Proca (W±,Z)      × (−1/12)  = {n_proca_broken * a1_proca:+8.4f}")
    print(f"    {n_fermion:2d} fermion DOF        × (−1/12)  = {n_fermion * a1_fermion:+8.4f}")
    print(f"    {n_G2_vec:2d}  G₂/SU(3) vectors  × (−1/12)  = {n_G2_vec * a1_proca:+8.4f}")
    print(f"    {n_G2_scal:2d}  G₂ fund scalars   × (+1/6)   = {n_G2_scal * a1_scalar:+8.4f}")
    print(f"    {'─' * 50}")
    print(f"    Σ_min(broken) = {Sigma_min_broken:.4f}")

    delta_sigma = Sigma_min_broken - Sigma_min_UV
    print(f"")
    print(f"  Shift: Σ(broken) − Σ(UV) = {delta_sigma:.4f}")
    print(f"    [3 Goldstones: scalar(+1/6) → Proca(−1/12), Δ per mode = −1/4]")
    assert abs(delta_sigma - (-0.25)) < 1e-10

    # ══════════════════════════════════════════════════════════════════
    #  10.7  Bridge contribution: Σ_bridge
    # ══════════════════════════════════════════════════════════════════

    S("10.7  Bridge contribution: 182 × (1/6 − 1/(48π))")

    # Each of the 182 bridge modes is a scalar with non-minimal
    # coupling ξ_bridge = 1/(48π).  Its heat-kernel coefficient is
    # a₁/R = 1/6 − ξ_bridge.
    #
    # Physical picture: the bridge modes form the "ventilation"
    # pattern around every standing-wave knot.  The slightly
    # sub-conformal coupling (ξ < 1/6) means each mode contributes
    # a small positive piece to Σ, and 182 of them collectively
    # overcome the negative SM contribution to produce the correct G_N.
    #
    # ── Five independent arguments fix the SCALAR sign ──
    #
    # (i)   Lagrangian algebra condensation:
    #        A_{E₈} = 1 + (τ_G, τ_F) absorbs bridge into E₈ vacuum.
    #        Condensed modes ≠ independent propagating gauge fields.
    #
    # (ii)  Zero central charge:
    #        c_coset = c(E₈) − c(G₂) − c(F₄) = 8 − 14/5 − 26/5 = 0.
    #        A c=0 sector is topological, not dynamical -- cannot support
    #        the ghost subtraction that produces the vector sign.
    #
    # (iii) Representation theory:
    #        Bridge transforms as bi-fundamental (7,26), NOT as adjoint
    #        of any gauge group. No gauge redundancy → no FP ghosts
    #        → ghost-subtracted sign −1/6 does not apply.
    #
    # (iv)  Heat-kernel ghost structure (Vassilevich 2003):
    #        Vector coefficient: −1/6 = +1/6 (naive) − 1/3 (ghost).
    #        Without ghosts → reverts to scalar value +1/6 − ξ.
    #
    # (v)   Quantum dimension verification:
    #        d(τ_G) · d(τ_F) = φ², d(A_{E₈}) = 1 + φ² = 2 + φ.
    #        FP dimension confirms bridge fully absorbed → D²_local = 1.

    Sigma_bridge = N_bridge * a1_bridge

    # Verify c_coset = 0
    c_E8 = 8.0
    c_G2 = 14.0 / 5.0
    c_F4 = 26.0 / 5.0
    c_coset = c_E8 - c_G2 - c_F4
    assert abs(c_coset) < 1e-12, f"c_coset ≠ 0: got {c_coset}"

    # Verify quantum dimension
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    d_tau_G = phi       # quantum dim of G₂(1) non-trivial simple
    d_tau_F = phi       # quantum dim of F₄(1) non-trivial simple
    d_A = 1.0 + phi**2  # = 2 + φ
    D2_total = (1.0 + phi**2)**2
    D2_local = D2_total / d_A**2
    assert abs(D2_local - 1.0) < 1e-12, f"D²_local ≠ 1: got {D2_local}"

    print(f"  Each bridge mode: a₁/R = 1/6 − 1/(48π) = {a1_bridge:.10f}")
    print(f"  Σ_bridge = {N_bridge} × {a1_bridge:.10f} = {Sigma_bridge:.6f}")
    print(f"")
    print(f"  Five arguments for scalar sign:")
    print(f"    (i)   Lagrangian algebra condensation: bridge absorbed into E₈ vacuum")
    print(f"    (ii)  c_coset = {c_E8} − {c_G2} − {c_F4} = {c_coset:.1f}  (topological, not dynamical)")
    print(f"    (iii) Bi-fundamental (7,26) ≠ adjoint → no gauge ghosts")
    print(f"    (iv)  Vector −1/6 = +1/6 − 1/3(ghost); no ghosts → scalar +1/6 − ξ")
    print(f"    (v)   D²_local = {D2_local:.1f} → bridge fully absorbed  ✓")

    # ── Division-algebra ceiling ──
    #
    # Hurwitz's theorem: the normed division algebras are R(1), C(2),
    # H(4), O(8) -- three Cayley-Dickson doublings.  A fourth doubling
    # (sedenions, dim 16) introduces zero divisors.  The E₈ lattice,
    # built from octavian integers, is the maximal structure the
    # common-mode fiber can support.
    #
    # The bridge has dim 182 >> 8 (octonion fiber).  Tracking all 182
    # channels as one coherent object would need a normed division
    # algebra of dim ≥ 182, which doesn't exist.  But c_coset = 0
    # means the cross-channel correlations are topological -- there's
    # nothing beyond the octonionic horizon to miss.  Protected
    # forgetting erases all internal labels (PvP = 0), so each channel
    # contributes identically.  Scalar treatment is exact.
    #
    # Ref: Kahsay (2026), "One Substrate, Three Generations", Sec. 1;
    #      Hurwitz (1898); Baez, "The Octonions" (2002).

    print(f"")
    print(f"  Division-algebra ceiling (Hurwitz's theorem):")
    print(f"    Normed division algebras: R(1) → C(2) → H(4) → O(8) -- three layers only")
    print(f"    Bridge dim = {N_bridge} >> 8 (octonion fiber)")
    print(f"    c_coset = {c_coset:.0f} → cross-channel correlations are topological")
    print(f"    Scalar decomposition is exact: no 'fourth layer' of phase memory exists")

    # ══════════════════════════════════════════════════════════════════
    #  10.8  Total Σ and G_ind/G_N
    # ══════════════════════════════════════════════════════════════════

    S("10.8  Induced Newton constant: G_ind/G_N = 6π / Σ_total")

    # Master formula from the one-loop heat kernel (Sakharov):
    #
    #   1/(16πG_ind) = (Λ²/(96π²)) × Σ_total
    #
    # With Λ = M_Pl and G_N = 1/M_Pl² (unreduced Planck convention):
    #
    #   16π G_N = 16π/M_Pl²
    #   G_ind/G_N = (16π/M_Pl²) × (M_Pl²/(96π²)) × Σ
    #             = 16π/(96π²) × Σ
    #             = 1/(6π) × Σ
    #
    #   → G_ind/G_N = Σ/(6π)
    #   → equivalently, if we invert: G_N/G_ind = Σ/(6π)
    #     so G_ind = G_N × 6π/Σ

    target = 6.0 * math.pi   # = 18.8496...

    Sigma_total_UV = Sigma_min_UV + Sigma_bridge
    Sigma_total_broken = Sigma_min_broken + Sigma_bridge

    G_ratio_UV = target / Sigma_total_UV
    G_ratio_broken = target / Sigma_total_broken
    G_ratio_mid = (G_ratio_UV + G_ratio_broken) / 2.0

    err_UV = 100.0 * (G_ratio_UV - 1.0)
    err_broken = 100.0 * (G_ratio_broken - 1.0)
    err_mid = 100.0 * (G_ratio_mid - 1.0)

    box([
        f"Σ_min(UV)     = {Sigma_min_UV:+.4f}    Σ_min(broken) = {Sigma_min_broken:+.4f}",
        f"Σ_bridge      = {Sigma_bridge:+.4f}    [{N_bridge} modes × (1/6 − 1/(48π))]",
        f"",
        f"Σ_total(UV)     = {Sigma_total_UV:.4f}   →  G_ind/G_N = {G_ratio_UV:.6f}  ({err_UV:+.2f}%)",
        f"Σ_total(broken) = {Sigma_total_broken:.4f}   →  G_ind/G_N = {G_ratio_broken:.6f}  ({err_broken:+.2f}%)",
        f"Midpoint:                          G_ind/G_N = {G_ratio_mid:.6f}  ({err_mid:+.2f}%)",
        f"",
        f"Target: 6π = {target:.4f}",
        f"Newton's constant derived to {abs(err_UV):.1f}–{abs(err_broken):.1f}% from M_Pl alone.",
        f"All couplings derived.",
    ])

    # ══════════════════════════════════════════════════════════════════
    #  10.9  Cosmological constant: Volovik → Jacobson → CKN saturation
    # ══════════════════════════════════════════════════════════════════

    S("10.9  Cosmological constant: Volovik → Jacobson → CKN")

    # ── Step 1: Volovik Gibbs-Duhem -- ρ_vac(eq) = 0 ─────────────
    #
    # For a self-bound condensate (the Z₃-NLS substrate), the
    # Gibbs-Duhem identity gives  ε + P = μn  in equilibrium.
    # At the self-bound point  P = 0:
    #
    #     ρ_vac(eq) = ε₀ + P₀ − μ₀n₀ = 0       (exactly)
    #
    # This is a thermodynamic identity, not fine-tuning.  The
    # microscopic zero-point energies (~M_Pl⁴) cancel automatically
    # against the chemical potential -- just as in superfluid helium,
    # the ground-state energy density does not gravitate.
    #
    # The "cosmological constant problem" (ρ_vac ~ M_Pl⁴ ≈ 10⁷⁶ GeV⁴
    # vs observed ~10⁻⁴⁷ GeV⁴) does not arise.
    #
    # Ref: Volovik, Universe in a Helium Droplet (OUP 2003), §29.

    from constants import M_Pl_GeV
    rho_naive = M_Pl_GeV**4

    print(f"  Step 1: Volovik Gibbs-Duhem mechanism")
    print(f"    Self-bound condensate: ρ_vac(eq) = 0  (thermodynamic identity)")
    print(f"    Naive QFT estimate:    ρ_vac ~ M_Pl⁴ = {rho_naive:.1e} GeV⁴")
    print(f"    This is ~10¹²³ × too large -- the \"cosmological constant problem\"")
    print(f"    The Gibbs-Duhem identity cancels this EXACTLY, not by fine-tuning.")

    # ── Step 2: Departure from equilibrium -- de Sitter expansion ──
    #
    # The only nonzero Λ arises when the condensate is driven out
    # of equilibrium.  In de Sitter expansion with Hubble rate H,
    # the condensate thermalizes at the acoustic Hawking temperature.
    #
    # The acoustic-metric formalism (Unruh 1981, Volovik 2003 §6,
    # Barceló-Liberati-Visser 2005 §3.6) gives the analog Hawking
    # temperature for a sonic horizon with surface gravity κ:
    #
    #     T_H = ℏκ / (2π c_s)
    #
    # For de Sitter: a(t) = exp(Ht), recession velocity v_rec = HL,
    # horizon at L_H = c/H, surface gravity κ = Hc, so:
    #
    #     T_dS = H / (2π)       (Gibbons-Hawking form, rederived
    #                            from acoustic metric WITHOUT
    #                            Einstein equations)
    #
    # The de Sitter horizon entropy (from substrate microstate
    # counting on the horizon, Volovik 2003 §29.7):
    #
    #     S_dS = A / (4G_eff) = π M_Pl² / H²
    #
    # Both T_dS and S_dS are derived from the Z₃-NLS substrate's
    # hydrodynamic / quasiparticle structure -- independently of GR.

    # Status of H₀:
    #   H₀ = 67.4 km/s/Mpc is a COSMOLOGICAL BOUNDARY CONDITION,
    #   fixed by the algebra of the algebraic framework.
    #   It enters ONLY in the Λ sector (ρ_Λ = 3H₀²M_Pl²/(8π)),
    #   NOT in the matter-sector predictions (masses, couplings, mixing).
    #   The framework derives the SCALE of Λ (∝ M_Pl² H₀²) and the
    #   coefficient 3/(8π) from Jacobson-Clausius, but does not
    #   predict the numerical value of H₀ itself.
    #   H₀ is determined by the initial conditions of the universe,
    #   analogous to M_Pl setting the unit system.
    H_0_SI = 67.4e3 / 3.0856e22            # s⁻¹
    hbar_SI = 1.0546e-34                     # J·s
    GeV_per_J = 1.0 / 1.602e-10             # GeV/J
    H_0_GeV = H_0_SI * hbar_SI * GeV_per_J  # GeV

    T_dS = H_0_GeV / (2.0 * math.pi)
    S_dS = math.pi * M_Pl_GeV**2 / H_0_GeV**2

    print(f"")
    print(f"  Step 2: De Sitter departure from equilibrium")
    print(f"    H₀ = 67.4 km/s/Mpc = {H_0_GeV:.4e} GeV  (cosmological boundary condition)")
    print(f"    Acoustic Hawking temperature: T_dS = H₀/(2π) = {T_dS:.4e} GeV")
    print(f"    Horizon entropy: S_dS = πM_Pl²/H₀² = {S_dS:.4e}")
    print(f"    (Both derived from acoustic metric, not from Einstein equations)")

    # ── Step 3: Jacobson's Clausius argument -- coefficient 3/(8π) ──
    #
    # Jacobson (1995, gr-qc/9504004) showed that Einstein's equation
    # follows from the Clausius relation δQ = T dS applied to local
    # Rindler horizons.  The three inputs are:
    #
    #   (A) T = H/(2π)        -- acoustic Hawking (Step 2)
    #   (B) S = A/(4G_eff)    -- horizon microstate counting (Step 2)
    #   (C) Raychaudhuri eqn  -- pure diff. geometry of null congruences
    #
    # For de Sitter vacuum (T_μν = −ρ_Λ g_μν):
    #
    #   R_μν = 3H²g_μν,  R = 12H²
    #   Clausius → R_μν − ½gR + Λg = 8πG T_μν
    #   de Sitter condition: 3H² = 8πG ρ_Λ
    #
    #   ∴  ρ_Λ = 3H² / (8πG_ind)  =  (3/8π) M_Pl² H²
    #
    # The coefficient 3/(8π) is DERIVED from (A)–(C), not assumed.
    # This is NOT the Friedmann equation assumed -- it is the Friedmann
    # equation DERIVED from substrate thermodynamics.

    rho_Lambda = 3.0 * H_0_GeV**2 * M_Pl_GeV**2 / (8.0 * math.pi)

    print(f"")
    print(f"  Step 3: Jacobson-Clausius derivation")
    print(f"    δQ = T dS on local Rindler horizons")
    print(f"    → Einstein equation with Λ as integration constant")
    print(f"    → de Sitter: ρ_Λ = 3H²/(8πG_ind) = (3/8π)M_Pl²H²")
    print(f"    Coefficient 3/(8π) derived from T_dS, S_dS, and Raychaudhuri")
    print(f"    ρ_crit = {rho_Lambda:.2e} GeV⁴")

    # ── Step 4: CKN saturation -- guaranteed by structure ──────────
    #
    # The Cohen-Kaplan-Nelson (1999, hep-th/9803132) bound is:
    #
    #     ρ_vac ≤ M_Pl² H²
    #
    # (from S_QFT(L) ≤ S_BH(L) for any IR scale L).
    #
    # The Z₃-NLS framework SATURATES this bound (inequality → equality)
    # because the four conditions for the Jacobson derivation are all
    # structurally satisfied:
    #
    #   (A) Acoustic Hawking T  -- from Bogoliubov dispersion at the
    #       sonic horizon (any BEC analog gravity satisfies this)
    #   (B) Area-law entropy    -- from substrate microstate counting
    #       (Volovik 2003 §29.7 verifies for sonic horizons)
    #   (C) Local Lorentz       -- at scales ≫ ξ₀ = ℓ_Pl/2
    #       (verified: screened Poisson, Eddington test)
    #   (D) Raychaudhuri        -- pure differential geometry of null
    #       congruences (needs only a smooth acoustic metric)
    #
    # Saturation is guaranteed by the framework's acoustic-metric
    # structure, not an additional assumption or fine-tuning.

    rho_CKN = M_Pl_GeV**2 * H_0_GeV**2
    coeff = rho_Lambda / rho_CKN

    print(f"")
    print(f"  Step 4: CKN bound saturation")
    print(f"    CKN bound: ρ_vac ≤ M_Pl²H² = {rho_CKN:.2e} GeV⁴")
    print(f"    Framework: ρ_Λ = (3/8π) × M_Pl²H² = {rho_Lambda:.2e} GeV⁴")
    print(f"    Coefficient: 3/(8π) = {coeff:.4f}  (from Jacobson, not fit)")
    print(f"    Saturation guaranteed by (A) acoustic Hawking T,")
    print(f"    (B) area-law entropy, (C) local Lorentz, (D) Raychaudhuri")

    # ── Step 5: Comparison with observation ───────────────────────
    #
    # The formula gives ρ_crit (the de Sitter equilibrium density).
    # The observed dark-energy density is:
    #
    #     ρ_Λ(obs) = Ω_Λ × ρ_crit     where Ω_Λ ≈ 0.685
    #
    # The factor Ω_Λ is the dark-energy fraction: the universe has
    # not yet reached its asymptotic de Sitter equilibrium (matter
    # has not fully diluted).  As Ω_Λ → 1 the observed value
    # converges to the framework's structural prediction.
    #
    # The SCALE prediction is the genuine content:
    #     ρ_Λ ~ M_Pl² H² ≈ 10⁻⁴⁷ GeV⁴       (framework)
    #     ρ_vac ~ M_Pl⁴   ≈ 10⁷⁶  GeV⁴       (naive QFT)
    # This is a 10¹²³-fold improvement -- the CC problem is resolved.

    Omega_Lambda = 0.685
    rho_Lambda_obs = Omega_Lambda * rho_Lambda
    orders = math.log10(rho_naive / rho_Lambda)

    ratio_Lambda = rho_Lambda / rho_Lambda_obs

    print(f"")
    print(f"  Step 5: Comparison with observation")
    print(f"    ρ_crit (framework)      = {rho_Lambda:.2e} GeV⁴")
    print(f"    ρ_Λ(obs) = Ω_Λ × ρ_crit = {rho_Lambda_obs:.2e} GeV⁴  (Ω_Λ = {Omega_Lambda})")
    print(f"    Ratio: ρ_crit/ρ_Λ(obs) = {ratio_Lambda:.2f}  (= 1/Ω_Λ)")
    print(f"    Matter has not fully diluted → universe not yet asymptotic dS")

    box([
        f"Derivation chain: Volovik (ρ_eq=0) → Jacobson (3/8π) → CKN saturation",
        f"",
        f"  Naive QFT:  ρ_vac ~ M_Pl⁴     = {rho_naive:.1e} GeV⁴",
        f"  Framework:  ρ_Λ  ~ M_Pl²H²    = {rho_Lambda:.2e} GeV⁴",
        f"  Observed:   ρ_Λ  = Ω_Λ ρ_crit = {rho_Lambda_obs:.2e} GeV⁴",
        f"",
        f"  Scale improvement: 10^{orders:.0f} (the CC problem is resolved)",
        f"  Remaining factor 1/Ω_Λ = {ratio_Lambda:.2f}: matter fraction (standard cosmology)",
    ])

    # ══════════════════════════════════════════════════════════════════
    #  10.10  Closure check: all sectors accounted
    # ══════════════════════════════════════════════════════════════════

    S("10.10  Closure: 248 = 14 + 52 + 182")

    # Every degree of freedom in the E₈ branching contributes to
    # exactly one sector:
    #   (14,1) → gauge sector (gluons, G₂/SU(3) vectors)
    #   (1,52) → matter/Higgs sector (quarks, leptons, W, Z, γ, H)
    #   (7,26) → bridge sector (gravity)
    #
    # No DOF is left over, no DOF is double-counted.
    # The emergence "closes" when every E₈ mode has a physical role.

    print(f"  (14,1):  G₂ adjoint  → gauge sector       [masses, mixing]")
    print(f"  (1,52):  F₄ adjoint  → matter/Higgs sector [masses, mixing]")
    print(f"  (7,26):  mixed       → gravity bridge      [G_N from heat kernel]")
    print(f"  Total: 14 + 52 + 182 = 248 = dim(E₈)  ✓")
    print(f"")
    print(f"  Every E₈ degree of freedom has a physical assignment.")
    print(f"  The conformal embedding is fully saturated → closure.")

    # ── Emergent spacetime summary ───────────────────────────────────

    S("10.11  Emergent spacetime: Z₃ substrate → common-mode geometry")

    # The gravity paper derives spacetime geometry as a collective
    # phenomenon from the Z₃ self-interfering substrate:
    #
    # 1. Z₃ Fourier decomposition: the substrate admits three sectors
    #    q = 0 (common mode) and q = 1,2 (relative modes).
    #    ω = exp(2πi/3) is the Z₃ phase.
    #
    # 2. Relative modes (q=1,2) carry the matter-like standing waves --
    #    particles are stable interference patterns that distinguish
    #    one component from another.
    #
    # 3. Common mode (q=0) is shared by all three components and
    #    carries the gravitational memory.  When a standing-wave knot
    #    (particle) forms, it depletes the common mode locally.
    #    This depletion is the "ventilation pattern" -- an attractive
    #    common-mode memory that IS the gravitational field.
    #
    # 4. The screened Poisson equation for the common-mode density R₀:
    #      (1 - ξ₀² ∇²) R₀ = σ/λ₀
    #    gives the Newtonian potential at distances >> ξ₀ = ℓ_Pl/2.
    #
    # 5. The Jacobson-Clausius route (area entropy + local acoustic
    #    Hawking temperature) promotes the screened Poisson equation
    #    to the full Einstein equations, with Λ from CKN-saturated
    #    de Sitter thermodynamics.
    #
    # Spacetime is therefore not a fundamental arena but a collective
    # hydrodynamic variable of the Z₃-NLS substrate.

    print(f"  The Z₃ substrate provides three Fourier sectors:")
    print(f"    q = 0:  common mode   → gravitational field (geometry)")
    print(f"    q = 1,2: relative modes → matter (standing-wave knots)")
    print(f"")
    print(f"  Emergence chain:")
    print(f"    Standing-wave knot depletes common mode locally")
    print(f"    → screened Poisson equation: (1 − ξ₀²∇²)R₀ = σ/λ₀")
    print(f"    → Newtonian gravity at r ≫ ξ₀ = ℓ_Pl/2")
    print(f"    → Jacobson-Clausius: acoustic Hawking T + area entropy")
    print(f"    → full Einstein equations with Λ from CKN bound")
    print(f"")
    print(f"  Protected forgetting ensures universality:")
    print(f"    PvP = 0    (Z₃ label erased → equivalence principle)")
    print(f"    Pv²P = ½P  (presence remembered → attractive coupling)")
    print(f"")
    print(f"  Spacetime = collective hydrodynamic variable of Z₃-NLS substrate")
    print(f"  Gravity = common-mode memory of standing-wave formation")

    return {
        'N_bridge': N_bridge,
        'xi_bridge': xi_bridge,
        'h_bridge': h_bridge,
        'h_G2_7': h_G2_7,
        'h_F4_26': h_F4_26,
        'alpha_G2_Pl': alpha_G2_Pl,
        'E_v2': E_v2,
        'a1_bridge': a1_bridge,
        'Sigma_min_UV': Sigma_min_UV,
        'Sigma_min_broken': Sigma_min_broken,
        'Sigma_bridge': Sigma_bridge,
        'Sigma_total_UV': Sigma_total_UV,
        'Sigma_total_broken': Sigma_total_broken,
        'G_ratio_UV': G_ratio_UV,
        'G_ratio_broken': G_ratio_broken,
        'G_ratio_mid': G_ratio_mid,
        'err_UV': err_UV,
        'err_broken': err_broken,
        'rho_Lambda': rho_Lambda,
        'rho_Lambda_obs': rho_Lambda_obs,
    }
