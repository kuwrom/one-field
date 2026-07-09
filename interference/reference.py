"""
reference.py: the ONE home of measured reference data.

Doctrine (enforced by tests/test_blindness.py):
  - Everything here is REFERENCE DATA (pulls, error columns)
    except the explicitly declared entries: PDG_EW['dr_hat_W'] and
    PDG_EW['rho_hat'] (cited radiative imports)
    and PDG_EW['M_Z'] (the alpha_s evaluation coordinate, shown
    immaterial: a 5 sigma shift moves alpha_s by < 3e-5).
  - The quarantine test perturbs every non-declared value at once and
    asserts the whole derivation chain re-runs identically to 1e-12.
  - The AST audit freezes every source location that reads these
    dicts.  A new read anywhere fails CI until admitted consciously.

The dimensional anchor m_e and the Lambda-sector reference H_0 live
in root.py (they are declared INPUTS, not comparisons).  root.py
re-exports these dicts for package consumers; both module attributes
bind the same dict objects.
"""

# External comparison labels only; they are not framework stages.
# Leptons: values conventionally tabulated as propagator-pole masses.
# Quarks: PDG 2024 conventional coordinates --
#   u, d, s: MS-bar at mu = 2 GeV.  c, b: MS-bar at mu = m.
#   t: direct kinematic mass reference.
# The algebra produces one value per fermion. It does not run to a
# scale.  The PDG coordinate each value matches is an observable
# mapping, not a second derivation.

PDG_MASSES = {
    'e': 0.51099895069, 'mu': 105.6583755, 'tau': 1776.93,
    'u': 2.16, 'c': 1273.0, 't': 172570.0,      # MeV
    'd': 4.70, 's': 93.5,  'b': 4183.0,
}

PDG_CKM = {
    'lambda': (0.22501, 0.00068),
    'A':      (0.826,   0.016),
    'rhobar': (0.1591,  0.0094),
    'etabar': (0.3523,  0.0073),
    'Vud': (0.97435, 0.00016),  'Vus': (0.22501, 0.00068),  'Vub': (0.00373, 0.00009),
    'Vcd': (0.22487, 0.00068),  'Vcs': (0.97349, 0.00016),  'Vcb': (0.04183, 0.00079),
    'Vtd': (0.00858, 0.00019),  'Vts': (0.04111, 0.00077),  'Vtb': (0.99912, 0.00003),
    'alpha_deg': (84.1, 4.5),  'beta_deg': (22.6, 0.5),  'gamma_deg': (65.7, 3.0),
    'J': (3.12e-5, 0.13e-5),  'sin2beta': (0.709, 0.011),
}
# The Wolfenstein parameters, nine magnitudes, and J are correlated outputs
# of the PDG global fit.  alpha, beta, gamma, and sin(2 beta) are separately
# quoted phase constraints, so this mixed block has no common covariance
# model.  Several PDG errors are asymmetric; the scalar stored here retains
# the positive-side uncertainty (rounded where shown).  mixing.py therefore
# reports approximate symmetric marginal pulls and a descriptive diagonal
# pull score, never a chi-square.

# Frozen NuFit 6.0 snapshot, normal ordering, IC19 without SK atmospheric data.
# This is intentionally a historical comparison coordinate; it is not a
# claim that 6.0 remains the latest NuFit release.
# The "IC19 without SK atmospheric data" column is the discriminator:
# sin2_23 = 0.561 there vs 0.470 with SK atmospheric.
# The stored errors are symmetric approximations to NuFit's asymmetric
# one-sigma intervals; correlations are omitted, so the displayed diagonal
# pull score is descriptive rather than a chi-square.
NUFIT_PMNS = {
    'sin2_12': (0.307, 0.012),
    'sin2_23': (0.561, 0.014),
    'sin2_13': (0.02195, 0.00056),
}

# References only: no interference/ code consumes these
# numbers (gravity.py's ν block is structural, no Δm² input; the
# normal-ordering benchmark combines them with m_1 = 0 for display).
NUFIT_OSCILLATION = {
    'dm2_21': (7.49e-5, 0.20e-5),
    'dm2_31': (2.513e-3, 0.027e-3),
}


# EW-sector comparison values and declared imports
# (PDG 2024 EW review, Erler-Freitas: Eq. 10.26 and the s-hat fit table)
PDG_EW = {
    'M_Z': 91.1876, 'M_W': 80.3692,
    'sin2W_MSbar': 0.23129, 'sin2W_err': 0.00004,   # global SM fit 2024
    'm_H': 125.20, 'alpha_s_MZ': 0.1180,
    # MS-bar radiative backreaction imports (DECLARED, cited).
    # (Standard QFT: "radiative corrections". Here they are backreactions
    # whose inputs are all re-derived from the framework.)
    # ZERO PARAMETER DISCIPLINE: these are functions of (G_F, m_t,
    # m_H, alpha-hat(M_Z)). Every ingredient is RECALCULATED from the
    # framework's own predictions and verified compatible within the
    # quoted errors (executable audit in gravity.py).  The only
    # data-driven ingredient is the dispersive hadronic VP inside
    # Delta-r0, a spectral measurement of the web.
    'dr_hat_W': 0.06937,   # ± 0.00006   PDG 2024 Eq. (10.26) text
    'rho_hat':  1.01016,   # ± 0.00009   PDG 2024 (incl. bosonic loops)
}

# Cosmological comparison data (single source)
PDG_COSMO = {
    'eta_B': 6.12e-10,     # baryon asymmetry, Planck 2018
    'Omega_Lambda': 0.685, # dark energy fraction, Planck 2018
}

# Proton-stability observational bound (single source; consumed by
# tests/probes/substrate_depth.py, the CONTINGENCY_SURFACE floor).
# Reference datum: sets a measured floor, steers no prediction.
SUPER_K = {
    'tau_p_yr': 1.6e34,        # Super-K partial lifetime bound
    'yr_s': 3.1557e7,          # Julian year, seconds
    'N_baryon': 1.0e80,        # baryons in the observable universe
    't_Hubble_s': 4.35e17,     # ~13.8 Gyr
}

# CMB density-fit comparison blocks (single source, consumed by
# dark_sector.py).  References only.
# Planck 2018 TT,TE,EE+lowE+lensing (Table 2, arXiv:1807.06209);
# the ratio Omega_c h^2 / Omega_b h^2 is H0-independent.
PLANCK_2018 = {
    'Omega_b_h2':     0.02237,     # ± 0.00015
    'Omega_b_h2_err': 0.00015,
    'Omega_c_h2':     0.1200,      # ± 0.0012
    'Omega_c_h2_err': 0.0012,
    'H_0_km_s_Mpc':   67.4,        # ± 0.5 (CC scale display only)
}
# ACT DR6 + CMB lensing + DESI DR1 (arXiv:2503.14452, March 2025)
ACT_DR6 = {
    'Omega_b_h2':     0.0226,      # ± 0.0001
    'Omega_b_h2_err': 0.0001,
    'Omega_c_h2':     0.118,       # ± 0.001
    'Omega_c_h2_err': 0.001,
}
