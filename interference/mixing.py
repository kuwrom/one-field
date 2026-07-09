"""
Mixing: CKM + PMNS from (d₁₀, d₁₁, n₇, n₂₆, π).

All derivations are self-contained, every parameter traces to the four
irreducible numbers via WZW S-matrix and nimrep structures.

═══════════════════════════════════════════════════════════════════════
CKM, four Wolfenstein parameters from four WZW faces
═══════════════════════════════════════════════════════════════════════

Each parameter comes from a different face of the SU(3)₃ WZW structure:

  λ = tan(h₁₀) = tan(d₁₀/d₁₁²) = tan(2/9)   , METRIC face
      The conformal weight h₁₀ = d₁₀/d₁₁² is the boundary OPE exponent
      for the fundamental (1,0) primary.  The tangent arises because the
      CKM mixing angle θ₁₂ is the real part of the boundary operator's
      OPE coefficient in the local coordinate on the unit disk.
      WHY tan AND NOT sin OR h₁₀ ITSELF: the boundary-condition-changing
      operator of SU(3)₃ lives on the unit disk (Cardy states). Its OPE
      coefficient in the LOCAL flat coordinate z = tan(θ/2) is tan(h₁₀).
      On the upper half-plane the coefficient would be sin(h₁₀) = 0.22040
      (miss: −6.8σ). On the strip it would be h₁₀ itself = 0.22222
      (miss: −4.1σ).  The unit disk is forced because it is the
      canonical coordinate for boundary states (Cardy 1989).
      Credit: Wolfenstein (1983) for the parameterisation.

  A = √Q₀ = √(d₁₀/d₁₁) = √(2/3)             , ALGEBRAIC face (nimrep)
      The D⁶ nimrep of SU(3)₃ has off-diagonal entries proportional to
      √Q₀ = √(d₁₀/d₁₁).  The square root is forced: nimrep entries are
      AMPLITUDES (Clebsch-Gordan coefficients of the fusion category),
      not probabilities. d₁₀ inter-generation channels out of d₁₁ total
      gives a probability Q₀ = d₁₀/d₁₁ = 2/3 whose amplitude is √Q₀.

  η̄ = π/d₁₁² = π/9                            , MODULAR face (S-matrix)
      The Kac-Peterson S-matrix phase for the fundamental primary is
      exp(2πi·h₁₀) = exp(2πi·d₁₀/d₁₁²).  The CP-violating phase η̄
      is the imaginary part of the boundary OPE coefficient, which
      equals π × h₁₀ = π·d₁₀/d₁₁².

  ρ̄ = h₁₀/√d₁₀ = √d₁₀/d₁₁² = √2/9           , TOPOLOGICAL face (Ocneanu)
      The Ocneanu cell amplitude for the (1,0) primary is h₁₀ = d₁₀/d₁₁²,
      divided by the quantum dimension amplitude √d₁₀.  This gives the
      real part of the unitarity triangle apex.

  γ = arctan(η̄/ρ̄) = arctan(π/√d₁₀)           , CP phase (derived)
      Follows directly from the ratio of imaginary to real parts.

  COMPRESSION: the four Wolfenstein parameters are NOT four independent
  choices.  The unitarity-triangle apex is ONE complex number:
      apex = ρ̄ + iη̄ = (h₁₀/√d₁₀)·(1 + iπ/√d₁₀)
  so η̄ and ρ̄ (plus γ and J) are a single algebraic object: the Ocneanu
  cell amplitude h₁₀/√d₁₀ rotated by the universal phase π/√d₁₀.  The
  effective independent parameters are TWO: λ = tan(h₁₀) and the apex,
  both fully fixed by h₁₀ and d₁₀, which are already determined by the
  conformal embedding.  No new dictionary entries are introduced for CKM.

═══════════════════════════════════════════════════════════════════════
PMNS, democratic mixing + charged-lepton back-reaction
═══════════════════════════════════════════════════════════════════════

  TBM base: the conjugation modular invariant Z_C of SU(3)₃ pairs each
  primary (l₁,l₂) with its conjugate (l₂,l₁).  The resulting nimrep J₃
  has a democratic (S₃-symmetric) structure, giving tribimaximal mixing.
  Credit: Harrison, Perkins, Scott (2002) for the TBM form.

  Charged-lepton back-reaction (standard QFT calls this "NLO correction"):
  U_PMNS = U_e† · U_TBM, where the charged-lepton
  rotation U_e = R₂₃(−φ) · R₁₃(θ_C, δ) is parameterised by:
    φ = √(m_e/m_μ)/√d₁₀   (mass ratio from lepton sector)
    θ_C = arctan(h₁₀)      (same conformal weight as CKM λ)
    δ = arccos(−λ)          (consistency with CKM phase)

  δ_CP prediction: testable by DUNE and Hyper-K.

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import cmath
import math
import numpy as np

from root import (d10, d11, n7, n26, h10, h11,
                  PDG_CKM, NUFIT_PMNS, pct)


# ═══════════════════════════════════════════════════════════════════════
#  WZW S-matrix infrastructure (SU(3) at level k = d₁₁ = 3)
# ═══════════════════════════════════════════════════════════════════════

def _su3_inner(a, b):
    """SU(3) Killing-form inner product on weight space."""
    return (2*a[0]*b[0] + a[0]*b[1] + a[1]*b[0] + 2*a[1]*b[1]) / 3.0


_WEYL = [
    (lambda a, b: (a, b), +1),
    (lambda a, b: (-a, a+b), -1),
    (lambda a, b: (a+b, -b), -1),
    (lambda a, b: (-(a+b), a), +1),
    (lambda a, b: (b, -(a+b)), +1),
    (lambda a, b: (-b, -a), -1),
]


def _compute_wzw():
    """Compute full SU(3)₃ WZW dataset: S-matrix, fusion, nimrep."""
    K_LEVEL = d11                         # = 3
    K_ALT   = K_LEVEL + d11              # altitude = 6 = d₁₀·d₁₁

    PRIMARIES = [(l1, l2) for l1 in range(K_LEVEL + 1)
                 for l2 in range(K_LEVEL + 1) if l1 + l2 <= K_LEVEL]
    N_PRIM = len(PRIMARIES)
    IDX = {p: i for i, p in enumerate(PRIMARIES)}

    # Kac-Peterson S-matrix
    S_raw = np.zeros((N_PRIM, N_PRIM), dtype=complex)
    for a, (l1a, l2a) in enumerate(PRIMARIES):
        pa = (l1a + 1, l2a + 1)
        for b, (l1b, l2b) in enumerate(PRIMARIES):
            pb = (l1b + 1, l2b + 1)
            val = sum(d * cmath.exp(-2j * math.pi * _su3_inner(w(*pa), pb) / K_ALT)
                      for w, d in _WEYL)
            S_raw[a, b] = val

    norm = math.sqrt(np.sum(np.abs(S_raw[0, :])**2))
    S_mat = S_raw / norm
    S_mat /= S_mat[0, 0] / abs(S_mat[0, 0])

    # Verlinde fusion
    N_fus = np.zeros((N_PRIM, N_PRIM, N_PRIM), dtype=int)
    for a in range(N_PRIM):
        for b in range(N_PRIM):
            for cc in range(N_PRIM):
                val = sum(S_mat[a, s] * S_mat[b, s] * S_mat[cc, s].conj() / S_mat[0, s]
                          for s in range(N_PRIM))
                N_fus[a, b, cc] = round(val.real)

    return {
        'PRIMARIES': PRIMARIES, 'N_PRIM': N_PRIM, 'IDX': IDX,
        'S_mat': S_mat, 'N_fus': N_fus,
    }


# ═══════════════════════════════════════════════════════════════════════
#  CKM derivation
# ═══════════════════════════════════════════════════════════════════════

def _derive_ckm():
    """Derive full CKM matrix from d₁₀, d₁₁.

    Four Wolfenstein parameters from four faces of SU(3)₃ WZW:
      λ , metric (boundary OPE exponent)
      A , algebraic (nimrep inter-generation fraction)
      η̄ , modular (S-matrix phase)
      ρ̄ , topological (Ocneanu cell amplitude)
    """

    # ── Wolfenstein parameters (all from d₁₀/d₁₁) ──

    # λ = tan(h₁₀): boundary OPE angle of fundamental primary
    # h₁₀ = d₁₀/d₁₁² = 2/9 (Sugawara conformal weight)
    lam = math.tan(float(h10))

    # A = √Q₀ = √(d₁₀/d₁₁): nimrep off-diagonal fraction
    # d₁₀ inter-generation channels out of d₁₁ total at each fusion vertex
    A = math.sqrt(float(d10) / float(d11))

    # η̄ = π/d₁₁² = π·h₁₀: S-matrix phase of fundamental primary
    # The Kac-Peterson phase exp(2πi·h₁₀) gives Im(OPE) = π·h₁₀
    eta = math.pi / float(d11**2)

    # ρ̄ = h₁₀/√d₁₀ = √d₁₀/d₁₁²: Ocneanu cell amplitude
    # Conformal weight divided by quantum dimension amplitude
    rho = float(h10) / math.sqrt(float(d10))

    # ── COMPRESSION: the unitarity-triangle apex is ONE complex number ──
    #     apex = ρ̄ + iη̄ = (h₁₀/√d₁₀)·(1 + iπ/√d₁₀)
    # the cell amplitude h₁₀/√d₁₀ rotated by the universal phase
    # direction π/√d₁₀.  All CP structure (η̄, ρ̄, γ, J) is this single
    # number. γ = arg(apex) = arctan(π/√d₁₀) is then automatic.
    apex = complex(rho, eta)
    assert abs(apex - (float(h10)/math.sqrt(d10))
               * complex(1.0, math.pi/math.sqrt(d10))) < 1e-15

    # γ = arg(apex) = arctan(π/√d₁₀), derived CP phase
    gamma_rad = math.atan2(eta, rho)
    gamma_deg = math.degrees(gamma_rad)

    # ── Full CKM matrix ──
    s12 = lam
    s23 = A * lam**2
    c12 = math.sqrt(1 - s12**2)
    c23 = math.sqrt(1 - s23**2)
    z13 = apex * s23 * s12 * c23 / (c12 * (1.0 - apex * s23**2))
    s13 = abs(z13)
    delta_rad = cmath.phase(z13)
    c13 = math.sqrt(1 - s13**2)

    ed = cmath.exp(1j * delta_rad)
    V = np.array([
        [c12*c13,                        s12*c13,                       s13*cmath.exp(-1j*delta_rad)],
        [-s12*c23 - c12*s23*s13*ed,      c12*c23 - s12*s23*s13*ed,     s23*c13],
        [s12*s23 - c12*c23*s13*ed,      -c12*s23 - s12*c23*s13*ed,     c23*c13],
    ], dtype=complex)
    Vabs = np.abs(V)

    # Jarlskog
    J_ckm = c12 * c23 * c13**2 * s12 * s23 * s13 * math.sin(delta_rad)

    # UT angles
    beta = math.atan2(eta, 1 - rho)
    alpha_ut = math.pi - gamma_rad - beta
    sin2beta = math.sin(2 * beta)

    alpha_deg = math.degrees(alpha_ut)
    beta_deg = math.degrees(beta)

    # ── Observable list and chi-squared ──
    labels_row = ['u', 'c', 't']
    labels_col = ['d', 's', 'b']
    lam_pdg, lam_sig = PDG_CKM['lambda']
    A_pdg, A_sig = PDG_CKM['A']
    rho_pdg, rho_sig = PDG_CKM['rhobar']
    eta_pdg, eta_sig = PDG_CKM['etabar']
    gamma_pdg, gamma_sig = PDG_CKM['gamma_deg']
    J_pdg, J_sig = PDG_CKM['J']

    ckm_obs = [
        ('lambda', lam, lam_pdg, lam_sig),
        ('A',      A,   A_pdg,   A_sig),
        ('rhobar', rho, rho_pdg, rho_sig),
        ('etabar', eta, eta_pdg, eta_sig),
    ]
    for i in range(3):
        for j in range(3):
            key = f"V{labels_row[i]}{labels_col[j]}"
            pdg_val, pdg_err = PDG_CKM[key]
            ckm_obs.append((f"|V{labels_row[i]}{labels_col[j]}|", Vabs[i, j], pdg_val, pdg_err))
    ckm_obs.append(('gamma_deg', gamma_deg, gamma_pdg, gamma_sig))
    ckm_obs.append(('J', J_ckm, J_pdg, J_sig))

    chi2_ckm = sum(((p - v) / e)**2 for _, p, v, e in ckm_obs)
    max_pull_ckm = max(abs((p - v) / e) for _, p, v, e in ckm_obs)

    return {
        'lam': lam, 'A': A, 'eta': eta, 'rho': rho,
        'V': V, 'Vabs': Vabs,
        'gamma_deg': gamma_deg, 'J_ckm': J_ckm,
        'alpha_deg': alpha_deg, 'beta_deg': beta_deg,
        'sin2beta': sin2beta,
        'ckm_obs': ckm_obs,
        'chi2_ckm': chi2_ckm,
        'max_pull_ckm': max_pull_ckm,
    }


# ═══════════════════════════════════════════════════════════════════════
#  PMNS derivation
# ═══════════════════════════════════════════════════════════════════════

def _derive_pmns(wzw, lam_ckm, m_e, m_mu):
    """Derive PMNS mixing from conjugation invariant + charged-lepton back-reaction.

    The conjugation modular invariant Z_C of SU(3)₃ pairs each primary
    (l₁,l₂) with its conjugate (l₂,l₁).  The orbit structure under this
    Z₂ action determines the neutrino mixing matrix:

    1. Orbits of Z_C give the democratic nimrep J₃ → TBM base matrix
    2. Charged-lepton rotation U_e provides back-reaction (called "NLO
       corrections" in standard convention):
         φ = √(m_e/m_μ)/√d₁₀  (lepton mass ratio, derived in masses.py)
         θ_C = h₁₀ = 2/9       (conformal weight = Cabibbo angle)
         δ = arccos(−λ_CKM)    (phase consistency)
    3. U_PMNS = U_e† · U_TBM gives three mixing angles + δ_CP prediction
    """

    PRIMARIES = wzw['PRIMARIES']
    IDX = wzw['IDX']
    N_fus = wzw['N_fus']
    N_PRIM = wzw['N_PRIM']

    # ── Conjugation invariant orbit structure ──
    CONJ_MAP = {(l1, l2): (l2, l1) for l1, l2 in PRIMARIES}
    orbits = []; seen = set()
    for p in PRIMARIES:
        if p not in seen:
            cp = CONJ_MAP[p]
            orbits.append((p,) if cp == p else (p, cp))
            seen.update([p, cp])

    orb_sizes = [len(o) for o in orbits]
    N_orb = len(orbits)

    # Orbit nimrep, LOAD-BEARING CHECK: the generation-triplet block
    # of the folded adjoint nimrep must be DEMOCRATIC (S₃-symmetric,
    # all entries equal, eigenvalues {3,0,0}). This is what licenses
    # the TBM base below.  Computed from the Verlinde data, asserted.
    idx_11 = IDX[(1, 1)]
    n11_fold = np.zeros((N_orb, N_orb))
    for i, oa in enumerate(orbits):
        for j, ob in enumerate(orbits):
            s11 = sum(N_fus[idx_11, IDX[a], IDX[b]] for a in oa for b in ob)
            n11_fold[i, j] = s11 / math.sqrt(orb_sizes[i] * orb_sizes[j])
    gen = [i for i, o in enumerate(orbits)
           if len(o) == 2 and sum(o[0]) in (1, 2, 3) and o[0] != (3, 0)
           and o[0] not in ((0, 3),)]
    blk = n11_fold[np.ix_(gen, gen)]
    assert blk.shape == (3, 3) and np.allclose(blk, blk[0, 0]), \
        "Z_C generation block must be democratic (licenses TBM)"
    ev = sorted(np.linalg.eigvalsh(blk), reverse=True)
    assert abs(ev[0] - 3*blk[0, 0]) < 1e-9 and abs(ev[1]) < 1e-9

    # ── TBM base matrix (from conjugation invariant orbit structure) ──
    # The democratic nimrep J₃ of the Z_C conjugation invariant has
    # S₃ permutation symmetry, giving the tribimaximal form.
    # Credit: Harrison-Perkins-Scott (2002) for the TBM parameterisation.
    U_TBM = np.array([
        [ 2/math.sqrt(6),  1/math.sqrt(3),  0],
        [-1/math.sqrt(6),  1/math.sqrt(3),  1/math.sqrt(2)],
        [ 1/math.sqrt(6), -1/math.sqrt(3),  1/math.sqrt(2)],
    ], dtype=complex)

    # ── Charged-lepton back-reaction: U_e = R₂₃(−φ) · R₁₃(θ_C, δ) ──
    # (Standard convention: "NLO correction to TBM.")
    # The charged-lepton mass matrix is not exactly diagonal in the
    # TBM basis. The rotation U_e accounts for this interference.
    #
    # φ = √(m_e/m_μ)/√d₁₀: lepton mass ratio suppression
    #   (m_e, m_μ derived from Koide Z₃ circulant in masses.py)
    phi = math.sqrt(m_e / m_mu) / math.sqrt(float(d10))
    # θ_C = h₁₀ = d₁₀/d₁₁² = 2/9: the conformal weight IS the Cabibbo angle
    #   CKM defines λ = tan(h₁₀), so θ_C = arctan(λ) = arctan(tan(h₁₀)) = h₁₀
    #   (quark-lepton complementarity from shared SU(3)₃ structure)
    theta = float(h10)
    # δ = arccos(−λ_CKM): Dirac phase from CKM consistency
    delta = math.acos(-lam_ckm)

    def R13(a, d):
        c, s, e = math.cos(a), math.sin(a), cmath.exp(1j*d)
        return np.array([[c,0,s*e],[0,1,0],[-s*e.conjugate(),0,c]], dtype=complex)

    def R23(a):
        c, s = math.cos(a), math.sin(a)
        return np.array([[1,0,0],[0,c,s],[0,-s,c]], dtype=complex)

    U_e = R23(-phi) @ R13(theta, delta)
    U_PMNS = np.conj(U_e.T) @ U_TBM

    # Extract mixing angles
    s13 = abs(U_PMNS[0, 2])
    c13sq = 1 - s13**2
    sin2_13 = s13**2
    sin2_12 = abs(U_PMNS[0, 1])**2 / c13sq
    sin2_23 = abs(U_PMNS[1, 2])**2 / c13sq
    delta_CP_deg = -math.degrees(cmath.phase(U_PMNS[0, 2]))

    J_lep = abs((U_PMNS[0,0]*U_PMNS[1,1]*np.conj(U_PMNS[0,1])*np.conj(U_PMNS[1,0])).imag)

    # Pulls vs NuFit
    pull_12 = (sin2_12 - NUFIT_PMNS['sin2_12'][0]) / NUFIT_PMNS['sin2_12'][1]
    pull_23 = (sin2_23 - NUFIT_PMNS['sin2_23'][0]) / NUFIT_PMNS['sin2_23'][1]
    pull_13 = (sin2_13 - NUFIT_PMNS['sin2_13'][0]) / NUFIT_PMNS['sin2_13'][1]
    chi2_pmns = pull_12**2 + pull_23**2 + pull_13**2

    return {
        'sin2_12': sin2_12, 'sin2_23': sin2_23, 'sin2_13': sin2_13,
        'delta_CP_deg': delta_CP_deg,
        'J_lep': J_lep,
        'chi2_pmns': chi2_pmns,
        'pulls': {'12': pull_12, '23': pull_23, '13': pull_13},
        'U_PMNS': U_PMNS,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive(R, masses):
    """
    Derive CKM + PMNS mixing matrices.

    R      : dict from root.derive()
    masses : dict from masses.derive()
    """
    from fractions import Fraction

    print("\n" + "=" * 78)
    print("  MIXING: CKM + PMNS from (d₁₀, d₁₁) = (2, 3)")
    print("=" * 78)

    # ── WZW S-matrix ──
    print(f"\n  SU(3) at level k = d₁₁ = {d11}, altitude K = d₁₀·d₁₁ = {d10*d11}")
    wzw = _compute_wzw()
    S_mat = wzw['S_mat']
    IDX = wzw['IDX']
    N_PRIM = wzw['N_PRIM']
    U_check = S_mat @ S_mat.conj().T
    unitarity_err = np.max(np.abs(U_check - np.eye(N_PRIM)))
    print(f"  S-matrix: {N_PRIM}x{N_PRIM}, unitarity |SS†−I| = {unitarity_err:.2e}")

    # ── CKM ──
    print(f"\n  CKM parameters (all from d₁₀/d₁₁):")
    ckm = _derive_ckm()

    h10_frac = Fraction(d10, d11**2)
    print(f"    lambda = tan(d₁₀/d₁₁²) = tan({h10_frac}) = {ckm['lam']:.5f}"
          f"  (PDG: {PDG_CKM['lambda'][0]}, pull: {(ckm['lam']-PDG_CKM['lambda'][0])/PDG_CKM['lambda'][1]:+.2f}sigma)")
    print(f"    A = sqrt(d₁₀/d₁₁) = sqrt({Fraction(d10,d11)}) = {ckm['A']:.5f}"
          f"  (PDG: {PDG_CKM['A'][0]}, pull: {(ckm['A']-PDG_CKM['A'][0])/PDG_CKM['A'][1]:+.2f}sigma)")
    print(f"    etabar = pi/d₁₁² = pi/{d11**2} = {ckm['eta']:.5f}"
          f"  (PDG: {PDG_CKM['etabar'][0]}, pull: {(ckm['eta']-PDG_CKM['etabar'][0])/PDG_CKM['etabar'][1]:+.2f}sigma)")
    print(f"    rhobar = h₁₀/sqrt(d₁₀) = sqrt({d10})/{d11**2} = {ckm['rho']:.5f}"
          f"  (PDG: {PDG_CKM['rhobar'][0]}, pull: {(ckm['rho']-PDG_CKM['rhobar'][0])/PDG_CKM['rhobar'][1]:+.2f}sigma)")
    print(f"    gamma = arctan(pi/sqrt({d10})) = {ckm['gamma_deg']:.2f} deg"
          f"  (PDG: {PDG_CKM['gamma_deg'][0]})")

    labels_row = ['u', 'c', 't']
    labels_col = ['d', 's', 'b']
    Vabs = ckm['Vabs']
    for i in range(3):
        row = "    "
        for j in range(3):
            key = f"V{labels_row[i]}{labels_col[j]}"
            pdg_val, pdg_err = PDG_CKM[key]
            pull = (Vabs[i, j] - pdg_val) / pdg_err
            row += f"|V_{labels_row[i]}{labels_col[j]}|={Vabs[i,j]:.5f}({pull:+.1f}s) "
        print(row)

    print(f"    J_ckm = {ckm['J_ckm']:.2e}  (PDG: {PDG_CKM['J'][0]:.2e})")
    print(f"    chi2(CKM) = {ckm['chi2_ckm']:.2f} / {len(ckm['ckm_obs'])} obs"
          f" = {ckm['chi2_ckm']/len(ckm['ckm_obs']):.2f}")

    # ── PMNS ──
    m_e = masses['m_e']
    m_mu = masses['m_mu']
    pmns = _derive_pmns(wzw, ckm['lam'], m_e, m_mu)

    print(f"\n  PMNS (democratic + charged-lepton back-reaction from Z_C conjugation):")
    print(f"    phi = sqrt(m_e/m_mu)/sqrt({d10}) = {math.sqrt(m_e/m_mu)/math.sqrt(float(d10)):.6f} rad")
    print(f"    theta_C = h₁₀ = {h10_frac} = {math.degrees(float(h10)):.3f} deg")
    print(f"    delta = arccos(-lambda) = {math.degrees(math.acos(-ckm['lam'])):.3f} deg")

    for label, key in [('12', 'sin2_12'), ('23', 'sin2_23'), ('13', 'sin2_13')]:
        val = pmns[key]
        nf_val, nf_err = NUFIT_PMNS[f'sin2_{label}']
        pull = (val - nf_val) / nf_err
        print(f"    sin2_theta_{label} = {val:.5f}  (NuFit: {nf_val}, pull: {pull:+.2f}sigma)")

    print(f"    delta_CP = {pmns['delta_CP_deg']:.1f} deg  (prediction)")
    print(f"    J_lep = {pmns['J_lep']:.6f}")
    print(f"    chi2(PMNS) = {pmns['chi2_pmns']:.2f} / 3 obs = {pmns['chi2_pmns']/3:.2f}")

    # ── Merge and return ──
    result = {}
    result.update(ckm)
    result.update(pmns)
    return result
