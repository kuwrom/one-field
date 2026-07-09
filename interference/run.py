#!/usr/bin/env python3
"""
End-to-end runner: every prediction from (d₁₀, d₁₁, n₇, n₂₆) + pi + m_e
(the electron anchors the scale. M_Pl and G are outputs).

Calls root -> masses -> mixing -> couplings -> gravity (incl. Higgs,
baryogenesis, neutrinos) in dependency order, then prints a compact
scorecard.

Usage:
    python3 interference/run.py

Reference:
    Kahsay, Kibrom Kidane (2026). One-field repository.
    Papers and archived releases are listed in README.md.
"""

import os
import sys

# ``python3 interference/run.py`` starts with interference/ on sys.path.
# Add the repository root and import the canonical package modules.
# Script and ``python -m`` entry modes use the same module objects.
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

# ═══════════════════════════════════════════════════════════════════════
#  Step 0: The wiring, found rather than stated.  Zero input: the
#  lattice, sub-algebra and level menus are scanned for zero closure
#  residual, and (n₇, n₂₆, d₁₀, d₁₁) are read off the winner.  This is
#  where the four integers and normalized current response come from.
# ═══════════════════════════════════════════════════════════════════════
from interference import wiring
wiring_data = wiring.derive()

# ═══════════════════════════════════════════════════════════════════════
#  Step 1: Algebraic root
# ═══════════════════════════════════════════════════════════════════════
from interference import root
assert root.WIRING_DATA == wiring_data, \
    "root must consume the canonical result produced by wiring.py"
R = root.derive()
from interference.root import PDG_MASSES, pct  # noqa: E402

# ═══════════════════════════════════════════════════════════════════════
#  Step 2: Generation word lemma (base(n) = boundary-walk counts;
#  runs BEFORE masses, which imports these integers)
# ═══════════════════════════════════════════════════════════════════════
from interference import words
words_data = words.derive()   # printed verification. the dict feeds no later step (dependency is by import, not by this value)

# ═══════════════════════════════════════════════════════════════════════
#  Step 3: Masses (9 fermions)
# ═══════════════════════════════════════════════════════════════════════
from interference import masses
masses_data = masses.derive(R)

# ═══════════════════════════════════════════════════════════════════════
#  Step 4: Mixing (CKM + PMNS)
# ═══════════════════════════════════════════════════════════════════════
from interference import mixing
mixing_data = mixing.derive(R, masses_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 5: Couplings (alpha_s + alpha(0))
# ═══════════════════════════════════════════════════════════════════════
from interference import couplings
couplings_data = couplings.derive(R, masses_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 6: Gravity + Higgs + Baryogenesis + Neutrinos
# ═══════════════════════════════════════════════════════════════════════
from interference import gravity
grav_data = gravity.derive(R, masses_data, mixing_data, couplings_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 7: Octonionic CG verification
# ═══════════════════════════════════════════════════════════════════════
from interference.searchlight import octonions
octonion_data = octonions.derive()   # printed verification; the dict feeds no later step (dependency is by import, not by this value)

# ═══════════════════════════════════════════════════════════════════════
#  Step 8: Dark sector (DM/baryon ratio + cosmological densities)
# ═══════════════════════════════════════════════════════════════════════
from interference import dark_sector
dark_data = dark_sector.derive(R, grav_data)

# ═══════════════════════════════════════════════════════════════════════
#  Step 9: Structural uniqueness proofs
# ═══════════════════════════════════════════════════════════════════════
from interference.searchlight import embedding_uniqueness
embed_data = embedding_uniqueness.derive()   # printed verification; the dict feeds no later step (dependency is by import, not by this value)

from interference.searchlight import protected_forgetting
pf_data = protected_forgetting.derive()   # printed verification; the dict feeds no later step (dependency is by import, not by this value)

# The Z₃-NLS substrate dynamics is the parallel verification
# program, not the canon: no scorecard number depends on it, so it
# lives in tests/substrate/ and is exercised by the probes, not here.


# ═══════════════════════════════════════════════════════════════════════
#  GLOBAL SCORECARD
# ═══════════════════════════════════════════════════════════════════════
#
# The structural and polynomial identities (d_11^2-1 = 4 d_10, the C_2
# factorisations, h_bridge = 1, c_coset = 0, the Fraction ratios 38/9,
# 217/18, 1165/12, 289/432, 32/27, 649/972, sin^2 theta_W = 3/13,
# and the E8 current trace Tr_+(Q^2) = 30 are asserted in
# test_interference.py and by inline asserts. These are regression
# certificates for the implemented formulas.

n_reported_numeric = 35  # includes internal coordinates and correlated comparisons
n_structural = 1  # conditional minimal-type-I rank(M_nu) <= 2

print("\n" + "=" * 78)
print(f"  GLOBAL SCORECARD, {n_reported_numeric} REPORTED NUMERICAL COORDINATES")
print("  CKM/PMNS COORDINATES ARE NOT INDEPENDENT; COVARIANCES ARE NOT MODELED")
print("=" * 78)

# # Mass sector #
print(f"\n  --- Mass sector (1 anchor + 8 predictions) ---")

mass_map = {
    'e':   masses_data['m_e'],   'mu':  masses_data['m_mu'],  'tau': masses_data['m_tau'],
    'u':   masses_data['m_u'],   'd':   masses_data['m_d'],   's':   masses_data['m_s'],
    'c':   masses_data['m_c'],   'b':   masses_data['m_b'],   't':   masses_data['m_t'],
}

methods = {
    'e': 'Z₃ Koide (Q₀=d₁₀/d₁₁)', 'mu': 'Z₃ Koide (Q₀=d₁₀/d₁₁)', 'tau': 'Z₃ Koide (Q₀=d₁₀/d₁₁)',
    'u': '(d₁₀²+h₁₀) m_e', 'd': '(d₁₁²+h₁₀) m_e',
    's': 'Q₀+h₁₀/K³ + √(32/27·(1+h₁₁/K²)) bridge',
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

# # CKM sector #
print(f"\n  --- CKM sector (15 comparison coordinates; correlations unmodeled) ---")
pull_score_ckm = mixing_data['diagonal_pull_score_ckm']
n_ckm = len(mixing_data['ckm_obs'])
print(f"  diagonal pull score = {pull_score_ckm:.2f}; "
      f"score/coordinate = {pull_score_ckm/n_ckm:.2f} (not chi2)")
print(f"  Max marginal pull: {mixing_data['max_pull_ckm']:.2f}sigma")

# # PMNS sector #
print(f"\n  --- PMNS sector (3 correlated angles + convention-dependent U_e3 phase) ---")
pull_score_pmns = mixing_data['diagonal_pull_score_pmns']
print(f"  diagonal pull score = {pull_score_pmns:.2f}; "
      f"score/angle = {pull_score_pmns/3:.2f} (not chi2)")
print(f"  -arg(U_e3) = {mixing_data['ue3_phase_deg']:.1f} deg  "
      "(current phase convention; not physical delta_CP)")

# # Neutrino sector #
print(f"\n  --- Neutrino sector (structural rank bound) ---")
print(f"  Identified RH neutrinos: {grav_data['n_RH']} from the F4 singlets")
print(f"  Conditional ({grav_data['neutrino_mass_mechanism_assumption']}): "
      f"rank(M_nu) <= {grav_data['light_mass_rank_max']}; "
      f"at least {grav_data['massless_light_states_min']} massless state")
print(f"  Ordering:   {grav_data['ordering']}  (NO -> m1=0; IO -> m3=0)")

# # Strong coupling #
print(f"\n  --- Strong coupling ---")
print(f"  alpha_s(M_Z) = {couplings_data['alpha_s_MZ_thresh']:.4f}  ({couplings_data['err_thresh']:+.2f}%)")
print(f"  strong CP: bar(theta) = "
      f"{couplings_data['bar_theta']:.1f} mod 2pi")

# # Gravity sector #
print(f"\n  --- Gravity sector ---")
print(f"  Bridge: {grav_data['N_bridge']} modes, xi = 1/(48pi)")
print(f"  G_ind/G_N (UV)     = {grav_data['G_ratio_UV']:.6f}  ({grav_data['err_UV']:+.2f}%)")
print(f"  G_ind/G_N (broken) = {grav_data['G_ratio_broken']:.6f}  ({grav_data['err_broken']:+.2f}%)")
print(f"  G_ind/G_N INTERNAL    = {grav_data['G_ratio_canonical']:.9f}")

# # Higgs transport #
print(f"\n  --- Higgs running coordinate ---")
print(f"  lambda(M_Pl) = -{grav_data['delta_bridge']:.6f}  (candidate bridge boundary)")
print(f"  lambda(m_t; beta coordinate) = {grav_data['lambda_mt_transport']:.6f}")
print(f"  sqrt(2lambda)v = {grav_data['higgs_transport_coordinate_GeV']:.2f} GeV")

# # The electroweak scale #
# Derived in masses.py and consumed by the three predictions below.
# Not one of the counted reported coordinates: it is the scale on which the
# electroweak entries stand. V_EW_GF is a reference, like the PDG values
# quoted below.
V_EW_GF = 246.219651
v_EW = masses_data['v_EW_pred_GeV']
print(f"\n  --- Electroweak scale ---")
print(f"  v_EW = M_Pl·exp(-S_quark) = {v_EW:.6f} GeV"
      f"  (from G_F: {V_EW_GF}, {(v_EW/V_EW_GF - 1)*1e6:+.2f} ppm)")

# # Electroweak #
print(f"\n  --- Electroweak sector (3 predictions: sin²θ_W, M_Z, M_W) ---")
sin2W = grav_data['sin2W']
print(f"  sin²θ_W = 3/13 + Q₀²·α/(2π) = {sin2W:.6f}  (PDG 2024 fit: 0.23129(4), pull {(sin2W-0.23129)/4e-5:+.2f}σ)")
print(f"  M_Z = {grav_data['M_Z_derived']:.4f} GeV  (PDG: {grav_data['M_Z_PDG']}, {pct(grav_data['M_Z_derived'], grav_data['M_Z_PDG']):+.4f}%)")
print(f"  M_W = {grav_data['M_W_derived']:.4f} GeV  (PDG: {grav_data['M_W_PDG']}(133), pull {(grav_data['M_W_derived']-grav_data['M_W_PDG'])/0.0133:+.2f}σ)")

# # Baryogenesis #
print(f"\n  --- Baryogenesis ---")
b_pct = abs(1 - grav_data['ratio_B']) * 100
print(f"  eta_B = {grav_data['dim_G2_fund']} * J_lep * exp(-2pi²) = {grav_data['eta_B']:.4e}")
print(f"  Planck: {grav_data['eta_B_obs']:.2e}  ({b_pct:.2f}% agreement)")

# # Dark sector #
print(f"\n  --- Dark sector (1 algebraic output + density bookkeeping) ---")
DM_b = dark_data['DM_baryon_ratio']
DM_b_obs = dark_data['DM_baryon_ratio_obs']
err_DM_b = dark_data['err_DM_baryon']
pull_DM_b = dark_data['pull_DM_baryon']
print(f"  Ω_DM/Ω_b = 2π−1 = {DM_b:.5f}  (Planck: {DM_b_obs:.4f}, {err_DM_b:+.2f}%, {pull_DM_b:+.2f}σ)")
print(f"  ρ_crit(H₀) = {dark_data['rho_crit_mapped']:.2e} GeV⁴; "
      f"ρ_Λ(current) = {dark_data['rho_Lambda_current']:.2e} GeV⁴ "
      "[H₀ and Ω_Λ imports]")
print(f"  Critical-density ratio to M_Pl⁴: "
      f"10^{dark_data['critical_density_suppression']:.0f}")


# ═══════════════════════════════════════════════════════════════════════
#  THE DERIVATION GRAPH
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
M_V = couplings_data['M_V_GeV']

print("\n" + "=" * 78)
print("  THE DERIVATION GRAPH")
print("=" * 78)

tree = f"""
  E8(1) -> G2(1) x F4(1)    [conformal embedding, c=8]
    |
    +-- G2 sector -> leptons
    |     B/A = sqrt(2), theta = 2/9, Q = 2/3
    |     (e, mu, tau): 3 masses
    |
    +-- F4 sector -> quarks + v_EW + Higgs
    |     v_EW = M_Pl * exp(-(9pi²/2 - 6 + 15/512 - 16(a/2pi)²))
    |     E8 current trace Tr_+(Q²) = 30  [derived in wiring.py]
    |     sin2W = 3/13 + Q₀²·α/(2π) = {sin2W:.5f}
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
    |     candidate lambda(M_Pl) -> lambda(m_t; beta coordinate) = {grav_data['lambda_mt_transport']:.4f}
    |
    +-- SU(3)_3 / D(6) -> CKM
    |     lambda = h₁₀(1+h₁₁/K²) = 73/324,  A = sqrt(2/3),
    |     modular phase pi/9 -> etabar [IH],  rhobar = sqrt(2)/9 [IH]
    |
    +-- SU(3)_3 conjugation-orbit candidate -> PMNS
    |     sin2_12 = {mixing_data['sin2_12']:.5f} ({p12:+.1f}sigma)
    |     sin2_23 = {mixing_data['sin2_23']:.5f} ({p23:+.1f}sigma)
    |     sin2_13 = {mixing_data['sin2_13']:.5f} ({p13:+.1f}sigma)
    |     -arg(U_e3) = {mixing_data['ue3_phase_deg']:.1f} deg  [convention-dependent; not delta_CP]
    |
    +-- F4 singlets -> neutrino structure
    |     26 -> 8_v + 8_s + 8_c + 1 + 1  ({grav_data['n_RH']} RH neutrinos)
    |     minimal type-I only: {grav_data['n_RH']} identified RH neutrinos -> rank(M_nu) <= {grav_data['light_mass_rank_max']}
    |     conditionally at least one massless state; ordering not selected
    |
    +-- Embedding index -> alpha_s
    |     alpha_G2(mu*) = pi/32,  mu* = M_Pl e^-(9pi^2/2-6) = 253.5 GeV
    |     M_V = {M_V:.0f} GeV  (per-mode bridge reading)
    |     alpha_s(M_Z) = {a_MZ:.4f}  ({e_th:+.1f}%)
    |     strong CP: bar(theta) = {couplings_data['bar_theta']:.1f} mod 2pi from the phase action on A=1+B
    |
    +-- G2 winding action + Fano orientation -> eta_B
    |     eta_B = 7 * J_lep * exp(-2pi²) = {grav_data['eta_B']:.4e}  ({b_pct:.1f}%)  [G2 face action]
    |
    +-- (7,26) bridge -> gravity
    |     182 modes, xi = 1/(48pi),  h_bridge = 1
    |     PvP = 0 (universality),  per-mode E[v²] = 1/2, total Pv²P = P (memory)
    |     G_ind/G_N = {grav_data['G_ratio_canonical']:.9f} (canonical internal closure)
    |
    \\-- Bridge self-interference -> dark sector
          Ω_DM/Ω_b = 2π−1 = {dark_data['DM_baryon_ratio']:.4f}  (Planck: {dark_data['DM_baryon_ratio_obs']:.3f}, {dark_data['pull_DM_baryon']:+.1f}σ)
          Friedmann density: ρ_crit = (3/8π) M_Pl²H₀²; ρ_Λ(current) = Ω_Λ ρ_crit
          Critical-density ratio to M_Pl⁴: 10^{dark_data['critical_density_suppression']:.0f}

  TOTAL: 9 mass coordinates (1 anchor + 8 predictions) + 15 CKM + 3 PMNS + 1 alpha_s + 1 Higgs + 1 G_N + 1 eta_B + 3 EW + 1 DM/baryon
       = {n_reported_numeric} reported numerical coordinates
         + {n_structural} conditional structural result (minimal-type-I light-neutrino rank bound)
       from the electron anchor.  The mixing theorem gives no relative
       diagonalizer from closed frame-blind data.  The finite maps are
       identifications.  Under the Planck/healing-length identification,
       the gravity target is Sigma = 2pi*d11.  Two factorizations evaluate
       to this target.  The substrate suite contains substrate-only clock
       runs, confined-phase lattice probes and the Skyrmion-Dirac spectral-flow probe.
       Structural and polynomial identities are asserted in test_interference.py:
           pytest -q   ->   the fast suite (the substrate arc behind -m slow), including
                            the canonical freeze table.  Counts: README.
  Mass: {n_sub1}/9 <= 1%, {n_sub2}/9 <= 2%, max {max_err:.1f}%
  CKM:  diagonal pull score/coordinate = {pull_score_ckm/n_ckm:.2f}  ({n_ckm} comparison coordinates; covariance unmodeled, not chi2)
  PMNS: diagonal pull score/angle = {pull_score_pmns/3:.2f}  (3 correlated angles; not chi2)
  G_N:  G_ind/G_N = {grav_data['G_ratio_canonical']:.4f}  [internal]
"""

print(tree)
print(f"  {'='*78}")
print(f"  CANONICAL ALGEBRAIC CLOSURE COMPLETE: ALL 248 E8 MODES ASSIGNED")
print(f"  {n_reported_numeric} REPORTED NUMERICAL COORDINATES, "
      f"{n_structural} CONDITIONAL STRUCTURAL RESULT; "
      "KEY THEOREM/IH/IMPORT BOUNDARIES SHOWN")
print(f"  SUBSTRATE RESULTS: CLOCK RUNS, CONFINED-PHASE LATTICE, CHIRAL ZERO MODES")
print(f"  (14,1) gauge + (1,52) matter/Higgs + (7,26) gravity = 248")
print(f"  {'='*78}")
