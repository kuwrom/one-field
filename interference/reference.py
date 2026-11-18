"""
reference.py: the ONE home of measured reference data.

Rules (enforced by tests/audit/test_blindness.py):
  - Everything here is REFERENCE DATA (pulls, error columns), except
    the explicitly declared entries.  These are PDG_EW['dr_hat_W'] and
    PDG_EW['rho_hat'] (cited radiative imports), and PDG_EW['M_Z']
    (the alpha_s evaluation coordinate, shown not to matter: a 5 sigma
    shift moves alpha_s by < 3e-5).
  - The quarantine test perturbs every nondeclared value at once and
    asserts that the whole derivation chain reruns identically to 1e-12.
  - The AST audit freezes every source location that reads these
    dicts.  A new read anywhere fails CI until it is deliberately
    admitted.

The anchor m_e is in root.py because it sets the unit.  H_0 is there
too, as the comparison scale of the critical density.
root.py reexports these dicts for package users.  Both module
attributes refer to the same dict objects.
"""

# External comparison labels only.  They are not framework stages.
# Leptons: values conventionally tabulated as propagator pole masses.
# Quarks: PDG 2024 conventional coordinates:
#   u, d, s: MS-bar at mu = 2 GeV.  c, b: MS-bar at mu = m.
#   t: direct kinematic mass reference.
# Selection rule: quark masses and quark mass ratios are the PDG values.
# A ratio whose lattice and dispersive determinations differ carries
# both, with the lattice value from FLAG for Nf = 2+1+1
# (PDG_LIGHT_RATIOS).  The algebra produces one value per fermion.  It
# does not run to a scale.  The PDG coordinate each value is compared
# at is a mapping to the observable, not a second derivation.

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
# model.  Several PDG errors are asymmetric.  The scalar stored here keeps
# the positive side uncertainty (rounded where shown).  mixing.py therefore
# reports approximate symmetric marginal pulls and a descriptive diagonal
# pull score, never a chi-square.

# Frozen NuFit 6.0 snapshot, normal ordering, IC19 without SK atmospheric data.
# This is a deliberately fixed comparison coordinate.  It makes no
# claim about which NuFit release is the latest.
# The "IC19 without SK atmospheric data" column is the discriminator:
# sin2_23 = 0.561 there vs 0.470 with SK atmospheric.
# The stored errors are symmetric approximations to NuFit's asymmetric
# one sigma intervals. Correlations are omitted, so the displayed diagonal
# pull score is descriptive rather than a chi-square.
NUFIT_PMNS = {
    'sin2_12': (0.307, 0.012),
    'sin2_23': (0.561, 0.014),
    'sin2_13': (0.02195, 0.00056),
}

# References only: no interference/ code uses these numbers.
# gravity.py's ν block is structural and takes no Δm² input.  The
# normal ordering benchmark combines them with m_1 = 0 for display.
NUFIT_OSCILLATION = {
    'dm2_21': (7.49e-5, 0.20e-5),
    'dm2_31': (2.513e-3, 0.027e-3),
}


# EW sector comparison values and declared imports
# (PDG 2024 EW review, Erler-Freitas: Eq. 10.26 and the s-hat fit table)
# Scheme invariant light quark ratios (comparison side only).
PDG_LIGHT_RATIOS = {
    'mu_over_md': (0.473, 0.017),        # PDG 2024
    'ms_over_mud': (27.30, 0.08),        # PDG 2024
    'Q_ellipse_dispersive': (22.1, 0.7), # eta -> 3 pi dispersive
    'Q_ellipse_lattice': (22.5, 0.5),    # FLAG 2024 Eq. 52, Nf = 2+1+1 (Nf = 2+1: 23.3(5), Eq. 51)
}

PDG_EW = {
    # M_W: the LHC-TeV MW Working Group world average without CDF II.
    'M_Z': 91.1876, 'M_W': 80.3625, 'M_W_err': 0.0077,
    'G_F': 1.1663788e-5, 'G_F_err': 0.0000006e-5,   # GeV^-2, PDG 2024. The one G_F in the repo
    'sin2W_MSbar': 0.23129, 'sin2W_err': 0.00004,   # global SM fit 2024
    'm_H': 125.20, 'alpha_s_MZ': 0.1180, 'alpha_s_MZ_err': 0.0009,
    # MS-bar radiative backreaction imports (DECLARED, cited).
    # (Backreactions whose inputs are all rederived from the framework.)
    # ZERO PARAMETER DISCIPLINE: these are functions of (G_F, m_t,
    # m_H, alpha hat(M_Z)).  G_F and m_t are RECALCULATED from the
    # framework's own predictions and checked to be compatible within
    # the quoted errors.  The coupling alpha hat(M_Z) is classed as
    # measurement data (the dispersive hadronic VP).  The framework
    # does not produce the Higgs mass m_H.  It enters the cited loop
    # remainders as an import at its measured value (executable audit
    # in higgs_ew.py).
    'dr_hat_W': 0.06937,   # ± 0.00006   PDG 2024 Eq. (10.26) text
    'rho_hat':  1.01016,   # ± 0.00009   PDG 2024 (incl. bosonic loops)
    # comparison values of the import audit (higgs_ew.py), display only:
    'rho_hat_err': 0.00009,       # PDG 2024 Eq. (10.26)
    'm_t_EW_fit': 172.61,         # GeV, the m_t inside the EW fit's rho hat
    'inv_alpha_hat_MZ': 127.929,  # dispersive 1/alpha hat(M_Z), MS-bar
}

# 1/alpha(0) determinations, comparison only: (value, standard uncertainty)
ALPHA_INV_REF = {
    'Cs_recoil': (137.035999046, 27e-9),   # Berkeley 2018
    'Rb_recoil': (137.035999206, 11e-9),   # LKB 2020
    'CODATA':    (137.035999177, 21e-9),   # CODATA 2022 blend
}

# Cosmological comparison data (single source)
PDG_COSMO = {
    'Omega_Lambda': 0.685, # dark energy fraction, Planck 2018
}

# Proton stability observational bound (single source, used by
# tests/probes/substrate_depth.py for the CONTINGENCY_SURFACE floor).
# Reference datum: it sets a measured floor and affects no prediction.
SUPER_K = {
    'tau_p_yr': 1.6e34,        # Super-K partial lifetime bound
    'yr_s': 3.1557e7,          # Julian year, seconds
    'N_baryon': 1.0e80,        # baryons in the observable universe
    't_Hubble_s': 4.35e17,     # ~13.8 Gyr
}

# CMB density fit comparison blocks (single source, used by
# dark_sector.py).  References only.
# Planck 2018 TT,TE,EE+lowE+lensing (Table 2, arXiv:1807.06209).
# The ratio Omega_c h^2 / Omega_b h^2 is independent of H0.
PLANCK_2018 = {
    'Omega_b_h2':     0.02237,     # ± 0.00015
    'Omega_b_h2_err': 0.00015,
    'Omega_c_h2':     0.1200,      # ± 0.0012
    'Omega_c_h2_err': 0.0012,
    'H_0_km_s_Mpc':   67.4,        # ± 0.5 (CC scale display only)
}
# ACT DR6 + CMB lensing + DESI DR1 (arXiv:2503.14452)
ACT_DR6 = {
    'Omega_b_h2':     0.0226,      # ± 0.0001
    'Omega_b_h2_err': 0.0001,
    'Omega_c_h2':     0.118,       # ± 0.001
    'Omega_c_h2_err': 0.001,
}
