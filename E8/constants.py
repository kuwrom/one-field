"""
Physical constants and experimental reference data.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The sole dimensional input is the Planck mass M_Pl, which sets
the unit system.  All dimensionless couplings are derived:

    alpha_{G_2}(v_EW) = pi/32       (exact WZW/Sugawara identity)
    alpha_s/alpha_em  = 16          (Singh, arXiv:2603.28810,
                                     from J_3(C x O) Jordan eigenvalues:
                                     charge-trace factor 8/3 times
                                     broken-phase dilution factor 6)
    alpha_em(algebra) = pi/512      (= (pi/32)/16, emergence scale)

    The physical fine-structure constant is a second algebraic identity:
    1/alpha(0) = (512/pi)(1 − 1/(2pi))
               = 256(2pi − 1)/pi^2
               = 137.036  (to 0.0003%)

    The factor (1 − 1/(2pi)) is an algebraic consequence of the
    conformal embedding; no perturbative dressing is applied.
    The instanton action uses alpha = pi/512; lepton pole masses
    use alpha(0).

The framework therefore has zero dimensionless free parameters.
Everything else below is PDG / NuFit reference data for comparison.
"""

import math

# ═══════════════════════════════════════════════════════════════════════
#  SOLE INPUT (dimensional -- sets the unit system)
# ═══════════════════════════════════════════════════════════════════════

M_Pl_GeV = 1.22089e19          # Planck mass  [GeV]
M_Pl_MeV = M_Pl_GeV * 1e3     # Planck mass  [MeV]

# ═══════════════════════════════════════════════════════════════════════
#  DERIVED COUPLINGS: two-scale alpha_em
# ═══════════════════════════════════════════════════════════════════════
#
#  Step 1 -- Algebraic coupling at the emergence scale:
#    alpha_{G_2}(v_EW) = pi/32     (Layer 7, exact WZW/Sugawara)
#    alpha_s(v_EW) = alpha_{G_2}   (embedding index 1)
#    alpha_s/alpha_em = 16          (Singh: (8/3) x 6 from J_3(C x O))
#    => alpha_em(algebraic) = pi/512
#
#  This is the coupling at the conformal embedding scale.
#  It enters the instanton vertex correction: delta_S = 30 × (pi/512)/(2pi) = 15/512.
#
#  Step 2 -- Physical fine-structure constant (algebraic identity):
#    1/alpha(0) = (512/pi) × (1 - 1/(2pi))
#               = 256(2pi - 1)/pi^2
#               = 137.0364  (expt: 137.036, agreement 0.0003%)
#
#    The (1 - 1/(2pi)) factor is an algebraic consequence of the
#    conformal embedding structure -- not an external perturbative
#    correction.  This alpha(0) enters the lepton pole-mass vertex
#    correction.

# -- Singh ratio alpha_s / alpha_em = 16 --
# Derived from the charge-trace structure of J_3(C x O).
# Two factors:
#
#   (a) Charge-trace factor = sum of Q^2 for one SM generation:
#         u-quark: Q = 2/3, N_c = 3  ->  3 x (2/3)^2 = 4/3
#         d-quark: Q = 1/3, N_c = 3  ->  3 x (1/3)^2 = 1/3
#         e-lepton: Q = 1             ->  1
#         neutrino: Q = 0             ->  0
#         Total = 4/3 + 1/3 + 1 = 8/3
#
#   (b) Broken-phase dilution factor = C_2(26, F_4) = 6
#       The six octonionic ladder-operator directions in the
#       Albert algebra J_3(O) that are frozen out at low energy.
#
#   Product: (8/3) x 6 = 16
#
# Ref: Singh, arXiv:2603.28810; independently verified in
#      Kahsay (2026), "One Substrate, Three Generations", Sec. 5.
Q2_u = 3.0 * (2.0/3.0)**2                           # = 4/3
Q2_d = 3.0 * (1.0/3.0)**2                           # = 1/3
Q2_e = 1.0 * 1.0**2                                 # = 1
Q2_nu = 0.0                                         # = 0
CHARGE_TRACE = Q2_u + Q2_d + Q2_e + Q2_nu           # = 8/3
C2_FUND_F4 = 6                                      # Casimir C_2(26, F_4)
SINGH_RATIO = CHARGE_TRACE * C2_FUND_F4             # = (8/3) x 6 = 16

ALPHA_G2_vEW = math.pi / 32.0        # exact WZW identity
ALPHA_EM = ALPHA_G2_vEW / SINGH_RATIO  # = pi/512 (emergence scale)

# Physical fine-structure constant (algebraic identity)
#   1/alpha(0) = (512/pi)(1 - 1/(2pi)) = 256(2pi-1)/pi^2
ALPHA_PHYS = math.pi**2 / (256.0 * (2.0 * math.pi - 1.0))

# Derived QED factors
QED_FACTOR = 1.0 - ALPHA_PHYS / (2.0 * math.pi)   # pole-mass vertex correction
QED_FACTOR_BARE = 1.0 - ALPHA_EM / (2.0 * math.pi)   # for reference (emergence scale)

# ===================================================================
#  DERIVED STRUCTURAL CONSTANTS (computed, not hardcoded)
# ===================================================================
#
#  These were previously written as bare numbers (2/9, sqrt(2), 30).
#  Here each is derived from the algebra so the code is end-to-end.

# -- Brannen angle theta = 2/9 --
# Sugawara formula for SU(3) at WZW level k = 3:
#   h_{(1,0)} = C_2(fund, SU(3)) / (k + h_dual(SU(3)))
#             = (4/3) / (3 + 3) = 2/9
# The conformal weight of the fundamental primary field of SU(3)_3
# fixes the phase in the Z_3 gap triplet.
# Ref: Kahsay (2026), "The Innocent Lepton", Sec. 4;
#      Di Francesco et al., "Conformal Field Theory" (1997), Eq. 15.60.
K_SU3 = 3                                          # WZW level (from conformal embedding)
H_DUAL_SU3 = 3                                     # dual Coxeter number of SU(3)
C2_FUND_SU3 = 4.0 / 3.0                            # quadratic Casimir C_2(3, SU(3))
THETA_BRANNEN = C2_FUND_SU3 / (K_SU3 + H_DUAL_SU3) # = (4/3)/6 = 2/9

# -- B/A ratio = sqrt(2) --
# From the Z_3 Bogoliubov-Koide identity (exact for any Z_3 spectrum):
#   Q = 1/3 + B^2/(6A^2)
# Corollary: Q = 2/3 iff |B/A| = sqrt(2).
# The G_2 cyclic BdG spectrum satisfies Q = 2/3 because the
# quantum dimension d_{(1,0)} = 2 of SU(3)_3 forces this ratio.
# Ref: Kahsay (2026), "The Koide Relation as an Exact Z_3
#      Bogoliubov Identity", Corollary 1;
#      Kahsay (2026), "The Innocent Lepton", Sec. 3.
BA_RATIO = math.sqrt(2.0)                          # |B/A| = sqrt(2)

# -- Vertex mode count N = 30 --
# F_4 fundamental under Spin(9): 26 = 1 + 9 + 16 (all modes contribute)
# plus the Higgs doublet: 2 complex = 4 real scalar DOFs.
# Ref: Kahsay (2026), "One Substrate, Three Generations", Sec. 3.
N_F4_FUND = 26                                     # dim of F_4 fundamental
N_HIGGS = 4                                        # real DOFs of Higgs doublet
N_VERTEX = N_F4_FUND + N_HIGGS                     # = 30

# ═══════════════════════════════════════════════════════════════════════
#  PDG 2025 REFERENCE VALUES  (comparison only -- not inputs)
# ═══════════════════════════════════════════════════════════════════════

PDG_MASSES = {                 # all in MeV
    'e': 0.51100, 'mu': 105.658, 'tau': 1776.93,
    'u': 2.16, 'c': 1273.0, 't': 172570.0,
    'd': 4.70, 's': 93.5,  'b': 4183.0,
}

PDG_CKM = {                   # Wolfenstein parameters + |V_ij| + UT angles
    'lambda': (0.22501, 0.00068),
    'A':      (0.826,   0.016),
    'rhobar': (0.1591,  0.0094),
    'etabar': (0.3523,  0.0073),
    'Vud': (0.97435, 0.00016),  'Vus': (0.22501, 0.00068),  'Vub': (0.00373, 0.00009),
    'Vcd': (0.22487, 0.00068),  'Vcs': (0.97349, 0.00016),  'Vcb': (0.04183, 0.00079),
    'Vtd': (0.00858, 0.00019),  'Vts': (0.04111, 0.00077),  'Vtb': (0.99912, 0.00003),
    'alpha_deg': (84.1, 4.5),
    'beta_deg':  (22.6, 0.5),
    'gamma_deg': (65.7, 3.0),
    'J':         (3.12e-5, 0.13e-5),
    'sin2beta':  (0.709, 0.011),
}

# NuFit 6.0 (Sep 2024), Normal Ordering
NUFIT_PMNS = {
    'sin2_12': (0.307, 0.012),
    'sin2_23': (0.561, 0.014),
    'sin2_13': (0.02195, 0.00056),
}

# NuFit 6.0 (Sep 2024), Normal Ordering -- mass-squared splittings
# Ref: Esteban et al., JHEP 2009, 178 (2020), arXiv:2007.14792;
#      NuFit 6.0 (2024), http://www.nu-fit.org/
NUFIT_OSCILLATION = {
    'dm2_21': (7.49e-5, 0.20e-5),     # Δm²₂₁ in eV²  (solar)
    'dm2_31': (2.513e-3, 0.027e-3),   # Δm²₃₁ in eV²  (atmospheric, NO)
}
