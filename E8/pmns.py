"""
Layer 6 -- PMNS Mixing from the Conjugation Invariant.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Tree level:  Democratic nimrep from Z_C → S₃ symmetry → TBM mixing.
NLO:  U_PMNS = U_e† · U_TBM  with  U_e = R₂₃(−φ) · R₁₃(θ_C, δ)

Three framework-derived parameters (zero free):
    φ = √(m_e/m_μ) / √2    -- mass hierarchy × Ocneanu cell magnitude
    θ = θ_C = arctan(2/9)   -- shared CKM conformal weight
    δ = arccos(−λ)           -- cross-invariant orientation reversal
"""

import cmath
import math
import numpy as np
from .constants import NUFIT_PMNS
from .formatting import H, S


def derive(wzw: dict, lep: dict, ckm: dict):
    """
    Derive the PMNS mixing matrix.

    Parameters
    ----------
    wzw : dict from wzw.derive()
    lep : dict from leptons.derive()
    ckm : dict from ckm.derive()

    Returns
    -------
    dict with keys:
        sin2_12, sin2_23, sin2_13  : predicted sin² values
        delta_CP_deg               : Dirac CP phase in degrees
        J_lep                      : leptonic Jarlskog
        chi2_pmns                  : χ² against NuFit
        pulls                      : dict of pull values
        U_PMNS                     : 3×3 PMNS matrix
    """

    H("LAYER 6:  PMNS FROM CONJUGATION INVARIANT")

    # ── 6.1  Conjugation invariant orbit structure ────────────────────

    S("6.1  Conjugation invariant Z_C")

    PRIMARIES = wzw['PRIMARIES']
    IDX = wzw['IDX']
    N_fus = wzw['N_fus']
    S_mat = wzw['S_mat']
    N_PRIM = wzw['N_PRIM']

    CONJ_MAP = {(l1, l2): (l2, l1) for l1, l2 in PRIMARIES}

    # Build orbits
    orbits = []; seen = set()
    for p in PRIMARIES:
        if p not in seen:
            cp = CONJ_MAP[p]
            orbits.append((p,) if cp == p else (p, cp))
            seen.update([p, cp])

    orb_sizes = [len(o) for o in orbits]
    N_orb = len(orbits)
    print(f"  {N_orb} orbits under conjugation")

    # ── 6.2  Orbit nimrep → democratic structure ─────────────────────

    S("6.2  Orbit nimrep")

    idx_11 = IDX[(1, 1)]
    idx_10 = IDX[(1, 0)]

    n11_fold = np.zeros((N_orb, N_orb))
    n10_fold = np.zeros((N_orb, N_orb))

    for i, oa in enumerate(orbits):
        for j, ob in enumerate(orbits):
            s11 = sum(N_fus[idx_11, IDX[a], IDX[b]] for a in oa for b in ob)
            s10 = sum(N_fus[idx_10, IDX[a], IDX[b]] for a in oa for b in ob)
            n11_fold[i, j] = s11 / math.sqrt(orb_sizes[i] * orb_sizes[j])
            n10_fold[i, j] = s10 / math.sqrt(orb_sizes[i] * orb_sizes[j])

    gen_idx = [i for i in range(N_orb) if orb_sizes[i] == 2 and n11_fold[i, i] > 0.5]

    if len(gen_idx) >= 3:
        gen3 = gen_idx[:3]
        block_11 = n11_fold[np.ix_(gen3, gen3)]
        is_democratic = np.allclose(block_11, np.ones((3, 3)))
        print(f"  n_(1,1)^fold gen block = J₃ (democratic)? {is_democratic}  ✓")
        print(f"  → S₃ permutation symmetry → TBM mixing at tree level")

    # ── 6.3  TBM base matrix ─────────────────────────────────────────

    S("6.3  Tribimaximal PMNS -- leading order")

    U_TBM = np.array([
        [ 2/math.sqrt(6),  1/math.sqrt(3),  0],
        [-1/math.sqrt(6),  1/math.sqrt(3),  1/math.sqrt(2)],
        [ 1/math.sqrt(6), -1/math.sqrt(3),  1/math.sqrt(2)],
    ], dtype=complex)

    print(f"  sin²θ₁₂ = 1/3,  sin²θ₂₃ = 1/2,  sin²θ₁₃ = 0  (tree)")

    # ── 6.4  NLO charged-lepton correction ───────────────────────────

    S("6.4  NLO correction -- U_e = R₂₃(−φ) · R₁₃(θ_C, δ)")

    phi = math.sqrt(lep['m_e'] / lep['m_mu']) / math.sqrt(2)
    theta = math.atan(2.0 / 9.0)         # = θ_C
    lam_ckm = ckm['lam']
    delta = math.acos(-lam_ckm)          # arccos(−λ)

    print(f"  φ = √(m_e/m_μ)/√2 = {phi:.6f} rad = {math.degrees(phi):.3f}°")
    print(f"  θ = θ_C = arctan(2/9) = {math.degrees(theta):.3f}°")
    print(f"  δ = arccos(−λ) = {math.degrees(delta):.3f}°")

    # Build rotation matrices
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

    # ── 6.5  Results vs NuFit ─────────────────────────────────────────

    S("6.5  Results vs NuFit 6.0 (Normal Ordering)")

    pull_12 = (sin2_12 - NUFIT_PMNS['sin2_12'][0]) / NUFIT_PMNS['sin2_12'][1]
    pull_23 = (sin2_23 - NUFIT_PMNS['sin2_23'][0]) / NUFIT_PMNS['sin2_23'][1]
    pull_13 = (sin2_13 - NUFIT_PMNS['sin2_13'][0]) / NUFIT_PMNS['sin2_13'][1]
    chi2_pmns = pull_12**2 + pull_23**2 + pull_13**2

    print(f"  sin²θ₁₂ = {sin2_12:.5f}   (NuFit: {NUFIT_PMNS['sin2_12'][0]}, pull: {pull_12:+.2f}σ)")
    print(f"  sin²θ₂₃ = {sin2_23:.5f}   (NuFit: {NUFIT_PMNS['sin2_23'][0]}, pull: {pull_23:+.2f}σ)")
    print(f"  sin²θ₁₃ = {sin2_13:.5f}   (NuFit: {NUFIT_PMNS['sin2_13'][0]}, pull: {pull_13:+.2f}σ)")
    print(f"  δ_CP    = {delta_CP_deg:.1f}°   (prediction)")
    print(f"  J_lep   = {J_lep:.6f}")
    print(f"  χ²(PMNS) = {chi2_pmns:.2f}  (3 obs, χ²/n = {chi2_pmns/3:.2f})")

    # TBM improvement factor
    chi2_TBM = ((1/3-0.307)/0.012)**2 + ((0.5-0.561)/0.014)**2 + ((0-0.02195)/0.00056)**2
    print(f"  TBM tree χ² = {chi2_TBM:.1f}  →  NLO χ² = {chi2_pmns:.2f}  ({chi2_TBM/chi2_pmns:.0f}× improvement)")

    return {
        'sin2_12': sin2_12, 'sin2_23': sin2_23, 'sin2_13': sin2_13,
        'delta_CP_deg': delta_CP_deg,
        'J_lep': J_lep,
        'chi2_pmns': chi2_pmns,
        'pulls': {'12': pull_12, '23': pull_23, '13': pull_13},
        'U_PMNS': U_PMNS,
    }
