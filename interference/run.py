#!/usr/bin/env python3
"""
End-to-end runner: every prediction from (d₁₀, d₁₁, n₇, n₂₆) + pi + m_e
(the electron anchors the scale; M_Pl and G are outputs).

Calls root -> masses -> mixing -> couplings -> gravity (incl. Higgs,
baryogenesis, neutrinos) in dependency order, then prints a compact
scorecard matching the v3_release format.

Usage:
    python run.py

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import sys
import os

# Ensure the package directory is on the import path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from root import PDG_MASSES, pct

# ═══════════════════════════════════════════════════════════════════════
#  Step 1: Algebraic root
# ═══════════════════════════════════════════════════════════════════════
import root
R = root.derive()

# ═══════════════════════════════════════════════════════════════════════
#  Step 2: Masses (9 fermions)
# ═══════════════════════════════════════════════════════════════════════
import masses
masses_data = masses.derive(R)

# ═══════════════════════════════════════════════════════════════════════
#  Step 3: Mixing (CKM + PMNS)
# ═══════════════════════════════════════════════════════════════════════
import mixing
mixing_data = mixing.derive(R, masses_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 4: Couplings (alpha_s + alpha(0))
# ═══════════════════════════════════════════════════════════════════════
import couplings
couplings_data = couplings.derive(R, masses_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 5: Gravity + Higgs + Baryogenesis + Neutrinos
# ═══════════════════════════════════════════════════════════════════════
import gravity
grav_data = gravity.derive(R, masses_data, mixing_data, couplings_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 6: Octonionic CG verification
# ═══════════════════════════════════════════════════════════════════════
import octonions
octonion_data = octonions.derive()

# ═══════════════════════════════════════════════════════════════════════
#  Step 6b: Generation word lemma (base(n) = boundary-walk counts)
# ═══════════════════════════════════════════════════════════════════════
import words
words_data = words.derive()

# ═══════════════════════════════════════════════════════════════════════
#  Step 7: Dark sector (DM/baryon ratio + CC scale)
# ═══════════════════════════════════════════════════════════════════════
import dark_sector
dark_data = dark_sector.derive(R, grav_data)


# ═══════════════════════════════════════════════════════════════════════
#  GLOBAL SCORECARD
# ═══════════════════════════════════════════════════════════════════════
#
# The structural and polynomial identities previously certified by
# proofs.py (d_11^2-1=4 d_10, C_2 factorisations, h_bridge=1, c_coset=0,
# the Fraction ratios 38/9, 217/18, 1165/12, 289/432, 32/27, 649/972,
# sin^2 theta_W = 3/13, N_vertex=30, etc.) are now asserted directly in
# test_interference.py and by inline asserts inside each module's derive()
# function.  The test suite is the certificate.

n_checks = 35   # 9 masses + 15 CKM + 3 PMNS + 1 alpha_s + 1 Higgs + 1 G_N + 1 eta_B + 3 EW + 1 DM/baryon
n_structural = 3  # m₁=0, normal ordering, CC scale 10⁻¹²³

print("\n" + "=" * 78)
print(f"  GLOBAL SCORECARD, {n_checks} NUMERICAL CHECKS, ZERO FIT PARAMETERS")
print("=" * 78)

# ── Mass sector ──
print(f"\n  --- Mass sector (9 predictions) ---")

mass_map = {
    'e':   masses_data['m_e'],   'mu':  masses_data['m_mu'],  'tau': masses_data['m_tau'],
    'u':   masses_data['m_u'],   'd':   masses_data['m_d'],   's':   masses_data['m_s'],
    'c':   masses_data['m_c'],   'b':   masses_data['m_b'],   't':   masses_data['m_t'],
}

methods = {
    'e': 'Z₃ Koide (Q₀=d₁₀/d₁₁)', 'mu': 'Z₃ Koide (Q₀=d₁₀/d₁₁)', 'tau': 'Z₃ Koide (Q₀=d₁₀/d₁₁)',
    'u': '(d₁₀²+h₁₀) m_e', 'd': '(d₁₁²+h₁₀) m_e', 's': 'Q₀+h₁₀/K³ + √(d₁₀³Q₀²/d₁₁) bridge',
    'c': '(d₁₀²d₁₁+δ) m_μ', 'b': 'Q₀+h₁₁/K³ Koide', 't': '(d₁₁⁴+d₁₀⁴+1/2K) m_τ',
}

print(f"\n  {'Particle':>8s}  {'Predicted':>12s}  {'PDG':>12s}  {'Error':>8s}  Method")
print(f"  {'─'*8}  {'─'*12}  {'─'*12}  {'─'*8}  {'─'*30}")

n_sub1, n_sub2, max_err = 0, 0, 0.0
for name in ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']:
    pred = mass_map[name]
    pdg = PDG_MASSES[name]
    err = pct(pred, pdg)
    if abs(err) <= 1: n_sub1 += 1
    if abs(err) <= 2: n_sub2 += 1
    max_err = max(max_err, abs(err))
    print(f"  {name:>8s}  {pred:12.4f}  {pdg:12.4f}  {err:+7.2f}%  {methods[name]}")

print(f"\n  {n_sub1}/9 within 1%,  {n_sub2}/9 within 2%,  max error {max_err:.1f}%")

# ── CKM sector ──
print(f"\n  --- CKM sector (15 observables) ---")
chi2_ckm = mixing_data['chi2_ckm']
n_ckm = len(mixing_data['ckm_obs'])
print(f"  chi2 = {chi2_ckm:.2f},  n = {n_ckm},  chi2/n = {chi2_ckm/n_ckm:.2f}")
print(f"  Max pull: {mixing_data['max_pull_ckm']:.2f}sigma")

# ── PMNS sector ──
print(f"\n  --- PMNS sector (3 observables + delta_CP prediction) ---")
chi2_pmns = mixing_data['chi2_pmns']
print(f"  chi2 = {chi2_pmns:.2f},  chi2/n = {chi2_pmns/3:.2f}")
print(f"  delta_CP = {mixing_data['delta_CP_deg']:.1f} deg  (testable by DUNE/Hyper-K)")

# ── Neutrino sector ──
print(f"\n  --- Neutrino sector (structural predictions: m1=0, ordering) ---")
print(f"  Structural: m1 = 0  (rank-{grav_data['n_RH']} seesaw from {grav_data['n_RH']} RH neutrinos in F4)")
print(f"  Ordering:   {grav_data['ordering']}  (m1 < m2 < m3)")
print(f"  Testable:   KATRIN endpoint, 0νββ non-observation, cosmological Σmν")

# ── Strong coupling ──
print(f"\n  --- Strong coupling ---")
print(f"  alpha_s(M_Z) = {couplings_data['alpha_s_MZ_thresh']:.4f}  ({couplings_data['err_thresh']:+.2f}%)")

# ── Gravity sector ──
print(f"\n  --- Gravity sector ---")
print(f"  Bridge: {grav_data['N_bridge']} modes, xi = 1/(48pi)")
print(f"  G_ind/G_N (UV)     = {grav_data['G_ratio_UV']:.6f}  ({grav_data['err_UV']:+.2f}%)")
print(f"  G_ind/G_N (broken) = {grav_data['G_ratio_broken']:.6f}  ({grav_data['err_broken']:+.2f}%)")
print(f"  G_ind/G_N (mid)    = {grav_data['G_ratio_mid']:.6f}")

# ── Higgs ──
print(f"\n  --- Higgs mass ---")
print(f"  lambda(M_Pl) = -{grav_data['delta_bridge']:.6f}  (F4 fusion + bridge threshold)")
print(f"  m_H = {grav_data['mH_pred']:.2f} GeV  (expt: 125.20 +/- 0.11)")

# ── Electroweak ──
print(f"\n  --- Electroweak sector (3 predictions: sin²θ_W, M_Z, M_W) ---")
sin2W = grav_data['sin2W']
print(f"  sin²θ_W = 3/13 + h₇(α/2π) = {sin2W:.6f}  (PDG 2024 fit: 0.23129(4), pull {(sin2W-0.23129)/4e-5:+.2f}σ)")
print(f"  M_Z = {grav_data['M_Z_derived']:.4f} GeV  (PDG: {grav_data['M_Z_PDG']}, {pct(grav_data['M_Z_derived'], grav_data['M_Z_PDG']):+.4f}%)")
print(f"  M_W = {grav_data['M_W_derived']:.4f} GeV  (PDG: {grav_data['M_W_PDG']}(13), pull {(grav_data['M_W_derived']-grav_data['M_W_PDG'])/0.0133:+.2f}σ)")

# ── Baryogenesis ──
print(f"\n  --- Baryogenesis ---")
b_pct = abs(1 - grav_data['ratio_B']) * 100
print(f"  eta_B = {grav_data['dim_G2_fund']} * J_lep * exp(-2pi²) = {grav_data['eta_B']:.4e}")
print(f"  Planck: {grav_data['eta_B_obs']:.2e}  ({b_pct:.2f}% agreement)")

# ── Dark sector ──
print(f"\n  --- Dark sector (1 prediction + 1 structural) ---")
DM_b = dark_data['DM_baryon_ratio']
DM_b_obs = dark_data['DM_baryon_ratio_obs']
err_DM_b = dark_data['err_DM_baryon']
pull_DM_b = dark_data['pull_DM_baryon']
print(f"  Ω_DM/Ω_b = 2π−1 = {DM_b:.5f}  (Planck: {DM_b_obs:.4f}, {err_DM_b:+.2f}%, {pull_DM_b:+.2f}σ)")
print(f"  CC scale: ρ_Λ ~ M_Pl²H₀² → 10^{dark_data['CC_suppression']:.0f} suppression  [structural]")


# ═══════════════════════════════════════════════════════════════════════
#  THE EMERGENCE TREE
# ═══════════════════════════════════════════════════════════════════════

m_c = masses_data['m_c']
m_t = masses_data['m_t']
m_b = masses_data['m_b']
m_s = masses_data['m_s']
p12 = mixing_data['pulls']['12']
p23 = mixing_data['pulls']['23']
p13 = mixing_data['pulls']['13']
a_MZ = couplings_data['alpha_s_MZ_thresh']
e_th = couplings_data['err_thresh']
M_V = couplings_data['M_V_derived']

print("\n" + "=" * 78)
print("  THE EMERGENCE TREE")
print("=" * 78)

tree = f"""
  E8(1) -> G2(1) x F4(1)    [conformal embedding, c=8]
    |
    +-- G2 sector -> leptons
    |     B/A = sqrt(2), theta = 2/9, Q = 2/3
    |     (e, mu, tau): 3 masses
    |
    +-- F4 sector -> quarks + v_EW + Higgs
    |     v_EW = M_Pl * exp(-(9pi²/2 - 6 + 15/512))
    |     30 modes = 26 (F4 fund) + 4 (Higgs DOFs)
    |     sin2W = d11/(d10²+d11²) = 3/13 = {sin2W:.5f}
    |     M_Z = {grav_data['M_Z_derived']:.2f} GeV,  M_W = {grav_data['M_W_derived']:.2f} GeV  (derived)
    |     alpha(emergence) = pi/512,  alpha(0) = 1/137.035999050
    |       ^ bridge self-interference: h=1, D²=1, c=0
    |     m_u = (38/9)m_e,  m_d = (83/9)m_e
    |
    |     Gen 1: m_u = (4+2/9) m_e
    |     Gen 2: m_c = (217/18) m_mu = {m_c:.1f} MeV
    |     Gen 3: m_t = (1165/12) m_tau = {m_t/1e3:.2f} GeV
    |
    |     Q = Q₀+h₁₁/K³ = 289/432 -> m_b = {m_b:.1f} MeV
    |     m_s = {m_s:.1f} MeV  (Q₀+h₁₀/K³ + √(d₁₀³Q₀²/d₁₁) bridge)
    |     lambda(M_Pl) = -delta_bridge -> m_H = {grav_data['mH_pred']:.1f} GeV
    |
    +-- SU(3)_3 / D(6) -> CKM
    |     lambda = tan(2/9),  A = sqrt(2/3),  etabar = pi/9,  rhobar = sqrt(2)/9
    |
    +-- SU(3)_3 / Z_C -> PMNS
    |     sin2_12 = {mixing_data['sin2_12']:.5f} ({p12:+.1f}sigma)
    |     sin2_23 = {mixing_data['sin2_23']:.5f} ({p23:+.1f}sigma)
    |     sin2_13 = {mixing_data['sin2_13']:.5f} ({p13:+.1f}sigma)
    |     delta_CP = {mixing_data['delta_CP_deg']:.1f} deg
    |
    +-- F4 singlets -> neutrino structure
    |     26 -> 8_v + 8_s + 8_c + 1 + 1  ({grav_data['n_RH']} RH neutrinos)
    |     Rank-{grav_data['n_RH']} seesaw -> m1 = 0  (normal ordering)
    |
    +-- Embedding index -> alpha_s
    |     alpha_G2(mu*) = pi/32,  mu* = M_Pl e^-(9pi^2/2-6) = 253.5 GeV
    |     M_V = {M_V:.0f} GeV  (derived)
    |     alpha_s(M_Z) = {a_MZ:.4f}  ({e_th:+.1f}%)
    |
    +-- G2 instanton + Fano orientation -> eta_B
    |     eta_B = 7 * J_lep * exp(-2pi²) = {grav_data['eta_B']:.4e}  ({b_pct:.1f}%)
    |
    +-- (7,26) bridge -> gravity
    |     182 modes, xi = 1/(48pi),  h_bridge = 1
    |     PvP = 0 (universality),  Pv²P = 1/2 P (memory)
    |     G_ind/G_N = {grav_data['G_ratio_UV']:.4f} (UV) - {grav_data['G_ratio_broken']:.4f} (broken)
    |
    \\-- Bridge self-interference -> dark sector
          Ω_DM/Ω_b = 2π−1 = {dark_data['DM_baryon_ratio']:.4f}  (Planck: {dark_data['DM_baryon_ratio_obs']:.3f}, {dark_data['pull_DM_baryon']:+.1f}σ)
          CC scale: ρ_Λ ~ M_Pl²H₀² (10^{dark_data['CC_suppression']:.0f} resolved)

  TOTAL: 9 masses + 15 CKM + 3 PMNS + 1 alpha_s + 1 Higgs + 1 G_N + 1 eta_B + 3 EW + 1 DM/baryon
       = {n_checks} numerical checks + {n_structural} structural (m1=0, ordering, CC scale)
       from the electron anchor alone (M_Pl derived), all couplings derived.
       Structural and polynomial identities are asserted in test_interference.py:
           pytest -q   ->   58 tests including the canonical freeze table.
  Mass: {n_sub1}/9 <= 1%, {n_sub2}/9 <= 2%, max {max_err:.1f}%
  CKM:  chi2/n = {chi2_ckm/n_ckm:.2f}  ({n_ckm} obs)
  PMNS: chi2/n = {chi2_pmns/3:.2f}  (3 obs)
  G_N:  G_ind/G_N = {grav_data['G_ratio_mid']:.4f}
"""

print(tree)
print(f"  {'='*78}")
print(f"  END-TO-END DERIVATION COMPLETE, ALL 248 E8 MODES ASSIGNED, {n_checks} CHECKS, {n_structural} STRUCTURAL")
print(f"  (14,1) gauge + (1,52) matter/Higgs + (7,26) gravity = 248")
print(f"  {'='*78}")
