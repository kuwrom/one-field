"""
Layer 5 -- CKM Matrix from Four WZW Structures.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Four Wolfenstein parameters, each from a distinct mathematical structure:
    λ = tan(h₁₀)       = tan(2/9)     ← metric   (conformal weight)
    A = √(n₁₁_off/d₁₁) = √(2/3)      ← algebraic (nimrep block)
    η̄ = arg(S_{f,f})    = π/9          ← modular  (S-matrix phase)
    ρ̄ = h₁₀ · |W|       = √2/9         ← topological (Ocneanu cells)

Physical motivation (readable without the papers)
──────────────────────────────────────────────────

The D⁽⁶⁾ modular invariant is the unique boundary extension of
SU(3)₃ selected by the Z₃ standing-wave condition.  It describes
how the three generations couple across the boundary of the WZW
theory (i.e., how they mix).  This single algebraic object has
four distinct mathematical faces, each of which determines one
Wolfenstein parameter:

λ = tan(2/9) — METRIC FACE (conformal weights):
    The conformal weight h(1,0) = 2/9 measures the scaling dimension
    of the lightest charged primary.  In the boundary OPE (Cardy/BPPZ
    formalism), chiral and anti-chiral sectors couple, and the Cabibbo
    mixing parameter is the ratio of boundary one-point functions:
    sin(h(1,0))/cos(h(1,0)) = tan(2/9).

A = √(2/3) — ALGEBRAIC FACE (nimrep):
    The D⁽⁶⁾ nimrep matrix n(1,1) has 2×2 blocks of the form [[1,2],[2,1]]
    for each generation pair.  The off-diagonal entry n_off = 2 counts the
    number of boundary channels connecting different generations.  The
    quantum dimension d(1,1) = 3 is the total number of channels.  The
    inter-generation fraction is 2/3, and A is its amplitude: √(2/3).

η̄ = π/9 — MODULAR FACE (S-matrix phase):
    The Kac–Peterson S-matrix of SU(3)₃ encodes the modular
    transformation τ → −1/τ.  The fundamental-fundamental entry
    S_{(1,0),(1,0)} has phase 11π/18 (from the Weyl orbit sum).
    After subtracting the standard phase normalisation π/2:
    η̄ = 11π/18 − π/2 = π/9.

ρ̄ = √2/9 — TOPOLOGICAL FACE (Ocneanu cells):
    The conformal weight h(1,0) = 2/9 supplies the length of the
    fundamental boundary step.  The unique Ocneanu cell (the
    topological datum of the D⁽⁶⁾ invariant, |W| = 1/√2) projects
    that step onto the real CP coordinate:
    ρ̄ = h(1,0) × |W| = (2/9) × (1/√2) = √2/9.

The CP phase γ = arctan(η̄/ρ̄) = arctan(π/√2) = 65.76° is then
the angle between the modular face and the topological face.
"""

import cmath
import math
import numpy as np
from .constants import PDG_CKM
from .formatting import H, S


def derive(wzw: dict):
    """
    Derive the full CKM matrix from SU(3)₃ WZW data.

    Parameters
    ----------
    wzw : dict from wzw.derive()

    Returns
    -------
    dict with keys:
        lam, A, eta, rho            : Wolfenstein parameters
        V, Vabs                     : CKM matrix and |V|
        gamma_deg, J_ckm            : CP phase and Jarlskog
        alpha_deg, beta_deg         : UT angles
        sin2beta                    : sin(2β)
        ckm_obs                     : list of (name, pred, pdg, err) tuples
        chi2_ckm, max_pull_ckm      : χ² statistics
    """

    h10 = wzw['h10']

    H("LAYER 5:  CKM MATRIX -- FOUR STRUCTURES → FOUR PARAMETERS")

    # The four Wolfenstein parameters below are NOT independent fits.
    # Each comes from a DISTINCT mathematical structure on SU(3)_3 / D^(6),
    # set by the conformal embedding.  See One Substrate, Secs. 9.1-9.4.

    # ── 5.1  λ = tan(h₁₀) ────────────────────────────────────────────
    # The Cabibbo angle from the METRIC face of the D⁽⁶⁾ invariant:
    # in the boundary OPE, chiral and anti-chiral sectors couple with
    # amplitudes sin(h) and cos(h), so the mixing ratio is tan(h).

    S("5.1  λ = tan(h₁₀) -- METRIC")

    # λ is the boundary-OPE phase tangent for the SU(3)_3 fundamental.
    # h_10 = 2/9 is the conformal weight (fixed by the conformal embedding);
    # tan of it is the leading Cabibbo-like rotation at the D^(6) boundary.
    lam = math.tan(h10)
    lam_pdg, lam_sig = PDG_CKM['lambda']
    print(f"  λ = tan(2/9) = {lam:.5f}   (PDG: {lam_pdg} ± {lam_sig}, pull: {(lam-lam_pdg)/lam_sig:+.2f}σ)")

    # ── 5.2  A = √(2/3) ──────────────────────────────────────────────
    # The hierarchy parameter from the ALGEBRAIC face:
    # n(1,1) nimrep blocks are [[1,2],[2,1]] for each generation pair.
    # Off-diagonal 2 out of quantum dimension 3 → fraction 2/3.
    # A is the amplitude (square root) of this fraction.

    S("5.2  A = √(2/3) -- ALGEBRAIC")

    # A = sqrt(off-diagonal fraction of the n_11 nimrep block of D^(6)).
    # The D^(6) boundary has 6 sites organised as 3 generation pairs, each
    # carrying the n_11 block [[1,2],[2,1]] (off-diag = 4, total = 6 → 2/3).
    A = math.sqrt(2.0 / 3.0)
    A_pdg, A_sig = PDG_CKM['A']
    print(f"  A = √(2/3) = {A:.5f}   (PDG: {A_pdg} ± {A_sig}, pull: {(A-A_pdg)/A_sig:+.2f}σ)")

    # ── 5.3  η̄ = π/9 ─────────────────────────────────────────────────
    # The CP-violating imaginary coordinate from the MODULAR face:
    # the Kac–Peterson S-matrix entry S_{(1,0),(1,0)} has phase 11π/18
    # from the Weyl orbit sum (computed in wzw.py).  After the standard
    # phase normalisation (subtract π/2): η̄ = 11π/18 − π/2 = π/9.

    S("5.3  η̄ = π/9 -- MODULAR (exact)")

    # η̄ is the argument of the Kac-Peterson S-matrix element S_{f,f}
    # for the SU(3)_3 fundamental primary.  Exact algebraic value.
    eta = math.pi / 9.0
    eta_pdg, eta_sig = PDG_CKM['etabar']
    print(f"  η̄ = π/9 = {eta:.5f}   (PDG: {eta_pdg} ± {eta_sig}, pull: {(eta-eta_pdg)/eta_sig:+.2f}σ)")

    # ── 5.4  ρ̄ = √2/9 ────────────────────────────────────────────────
    # The CP-violating real coordinate from the TOPOLOGICAL face:
    # h(1,0) = 2/9 is the boundary step length (conformal weight),
    # and |W| = 1/√2 is the unique Ocneanu cell amplitude of D⁽⁶⁾.
    # The Ocneanu cell projects the conformal step onto the real axis:
    # ρ̄ = h(1,0) × |W| = (2/9)/(√2) = √2/9.

    S("5.4  ρ̄ = h₁₀/√2 = √2/9 -- TOPOLOGICAL")

    # ρ̄ = h_10 * |W| where |W| = 1/sqrt(2) is the Ocneanu cell norm.
    # Type-I unitarity of the D^(6) cells fixes |W| with zero physical DOF.
    rho = h10 / math.sqrt(2)
    rho_pdg, rho_sig = PDG_CKM['rhobar']
    print(f"  ρ̄ = √2/9 = {rho:.5f}   (PDG: {rho_pdg} ± {rho_sig}, pull: {(rho-rho_pdg)/rho_sig:+.2f}σ)")

    # ── 5.5  CP phase γ ──────────────────────────────────────────────
    # γ is the angle between the modular face (η̄ = π/9) and the
    # topological face (ρ̄ = √2/9) of the same D⁽⁶⁾ algebra.
    # γ = arctan(η̄/ρ̄) = arctan((π/9)/(√2/9)) = arctan(π/√2) = 65.8°.

    S("5.5  γ = arctan(π/√2)")

    gamma_rad = math.atan2(eta, rho)
    gamma_deg = math.degrees(gamma_rad)
    gamma_pdg, gamma_sig = PDG_CKM['gamma_deg']
    print(f"  γ = {gamma_deg:.2f}°   (PDG: {gamma_pdg}° ± {gamma_sig}°, pull: {(gamma_deg-gamma_pdg)/gamma_sig:+.2f}σ)")

    # ── 5.6  Full CKM matrix ─────────────────────────────────────────

    S("5.6  Full CKM matrix")

    s12 = lam
    s23 = A * lam**2
    c12 = math.sqrt(1 - s12**2)
    c23 = math.sqrt(1 - s23**2)
    Rbar = complex(rho, eta)
    z13 = Rbar * s23 * s12 * c23 / (c12 * (1.0 - Rbar * s23**2))
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

    labels_row = ['u', 'c', 't']
    labels_col = ['d', 's', 'b']
    for i in range(3):
        row = "  "
        for j in range(3):
            key = f"V{labels_row[i]}{labels_col[j]}"
            pdg_val, pdg_err = PDG_CKM[key]
            pull = (Vabs[i, j] - pdg_val) / pdg_err
            row += f"  |V_{labels_row[i]}{labels_col[j]}| = {Vabs[i,j]:.5f} ({pull:+.2f}σ)"
        print(row)

    # Jarlskog
    J_ckm = c12 * c23 * c13**2 * s12 * s23 * s13 * math.sin(delta_rad)
    J_pdg, J_sig = PDG_CKM['J']

    # UT angles
    beta = math.atan2(eta, 1 - rho)
    alpha_ut = math.pi - gamma_rad - beta
    sin2beta = math.sin(2 * beta)

    alpha_deg = math.degrees(alpha_ut)
    beta_deg = math.degrees(beta)

    print(f"\n  J = {J_ckm:.2e}  (PDG: {J_pdg:.2e}, pull: {(J_ckm-J_pdg)/J_sig:+.2f}σ)")
    print(f"  α = {alpha_deg:.1f}°,  β = {beta_deg:.1f}°,  γ = {gamma_deg:.1f}°")
    print(f"  sin2β = {sin2beta:.3f}  (PDG: {PDG_CKM['sin2beta'][0]})")

    # ── Build observable list and χ² ──────────────────────────────────

    ckm_obs = [
        ('λ',       lam,        lam_pdg,      lam_sig),
        ('A',       A,          A_pdg,        A_sig),
        ('ρ̄',      rho,        rho_pdg,      rho_sig),
        ('η̄',      eta,        eta_pdg,      eta_sig),
    ]
    for i in range(3):
        for j in range(3):
            key = f"V{labels_row[i]}{labels_col[j]}"
            pdg_val, pdg_err = PDG_CKM[key]
            ckm_obs.append((f"|V{labels_row[i]}{labels_col[j]}|", Vabs[i,j], pdg_val, pdg_err))
    ckm_obs.append(('γ (°)', gamma_deg, gamma_pdg, gamma_sig))
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
