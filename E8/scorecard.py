"""
Global Scorecard -- combined results and emergence tree.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Collects all predictions from every layer and prints the final
summary: mass sector, CKM, PMNS, α_s, gravity, and the emergence tree.
"""

import math
from .constants import PDG_MASSES
from .formatting import H, S, box, pct


def derive(lep: dict, quarks: dict, ckm: dict, higgs: dict, pmns: dict, alpha_s: dict, proofs=None, grav=None, nu=None, baryo=None):
    """
    Print the global scorecard and emergence tree.

    Parameters
    ----------
    lep     : dict from leptons.derive()
    quarks  : dict from quarks.derive()
    ckm     : dict from ckm.derive()
    higgs   : dict from higgs.derive()
    pmns    : dict from pmns.derive()
    alpha_s : dict from alpha_s.derive()
    proofs  : dict from proofs.derive() (optional)
    grav    : dict from gravity.derive() (optional)
    nu      : dict from neutrinos.derive() (optional)
    """

    n_checks = 30 if grav else 29
    # Neutrino sector adds: m₁=0 (structural) + ordering + Σmν
    if nu is not None:
        n_checks += 3
    # Baryogenesis adds: η_B
    if baryo is not None:
        n_checks += 1
    H(f"GLOBAL SCORECARD -- {n_checks} NUMERICAL CHECKS, ZERO FIT PARAMETERS")

    # ── Mass sector ──────────────────────────────────────────────────

    S("Mass sector (9 predictions)")

    mass_predictions = {
        'e':   (lep['lepton_pred']['e'],   'Brannen Z₃'),
        'mu':  (lep['lepton_pred']['mu'],  'Brannen Z₃'),
        'tau': (lep['lepton_pred']['tau'], 'Brannen Z₃'),
        'u':   (quarks['m_u'],             '(38/9) m_e'),
        'd':   (quarks['m_d'],             '(83/9) m_e'),
        's':   (quarks['m_s'],             'Rivero Q=2/3+h₁₀/K³ + √(32/27)'),
        'c':   (quarks['m_c'],             '(217/18) m_μ'),
        'b':   (quarks['m_b'],             'Q = 2/3 + h₁₁/K³'),
        't':   (quarks['m_t'],             '(1165/12) m_τ (emergence)'),
    }

    print(f"\n  {'Particle':>8s}  {'Predicted':>12s}  {'PDG':>12s}  {'Error':>8s}  Method")
    print(f"  {'─'*8}  {'─'*12}  {'─'*12}  {'─'*8}  {'─'*30}")

    n_sub1, n_sub2, max_err = 0, 0, 0.0
    for name in ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']:
        pred, method = mass_predictions[name]
        pdg = PDG_MASSES[name]
        err = pct(pred, pdg)
        if abs(err) <= 1: n_sub1 += 1
        if abs(err) <= 2: n_sub2 += 1
        max_err = max(max_err, abs(err))
        print(f"  {name:>8s}  {pred:12.4f}  {pdg:12.4f}  {err:+7.2f}%  {method}")

    print(f"\n  {n_sub1}/9 within 1%,  {n_sub2}/9 within 2%,  max error {max_err:.1f}%")

    # ── CKM sector ───────────────────────────────────────────────────

    S("CKM sector (15 observables)")

    chi2_ckm = ckm['chi2_ckm']
    n_ckm = len(ckm['ckm_obs'])
    print(f"  χ² = {chi2_ckm:.2f},  n = {n_ckm},  χ²/n = {chi2_ckm/n_ckm:.2f}")
    print(f"  Max pull: {ckm['max_pull_ckm']:.2f}σ")

    # ── PMNS sector ──────────────────────────────────────────────────

    S("PMNS sector (3 observables + δ_CP prediction)")

    chi2_pmns = pmns['chi2_pmns']
    print(f"  χ² = {chi2_pmns:.2f},  χ²/n = {chi2_pmns/3:.2f}")
    print(f"  δ_CP = {pmns['delta_CP_deg']:.1f}°  (testable by DUNE/Hyper-K)")

    # ── Neutrino sector ──────────────────────────────────────────────

    if nu is not None:
        S("Neutrino sector (3 predictions: m₁=0, ordering, Σmν)")

        print(f"  Structural: m₁ = 0  (rank-2 seesaw from 2 RH neutrinos in F₄)")
        print(f"  Ordering:   normal  (m₁ < m₂ < m₃)")
        print(f"  m₂ = {nu['m2_meV']:.2f} meV,  m₃ = {nu['m3_meV']:.2f} meV")
        print(f"  Σmν = {nu['sum_m_nu_meV']:.1f} ± {nu['sum_m_nu_err_meV']:.1f} meV  (testable by DESI/Euclid)")

    # ── α_s ──────────────────────────────────────────────────────────

    S("Strong coupling")

    print(f"  α_s(M_Z) = {alpha_s['alpha_s_MZ_thresh']:.4f}  ({alpha_s['err_thresh']:+.2f}%)")

    # ── Gravity sector ───────────────────────────────────────────────

    if grav is not None:
        S("Gravity sector (1 prediction: G_ind/G_N)")

        print(f"  Bridge sector: ({grav['N_bridge']} modes) × (1/6 − 1/(48π))")
        print(f"  ξ_bridge = 1/(48π) = {grav['xi_bridge']:.8f}")
        print(f"  G_ind/G_N (UV)     = {grav['G_ratio_UV']:.6f}  ({grav['err_UV']:+.2f}%)")
        print(f"  G_ind/G_N (broken) = {grav['G_ratio_broken']:.6f}  ({grav['err_broken']:+.2f}%)")
        print(f"  G_ind/G_N (mid)    = {grav['G_ratio_mid']:.6f}")

    # ── Baryogenesis sector ─────────────────────────────────────────

    if baryo is not None:
        S("Baryogenesis (1 prediction: η_B)")

        pull_pct = abs(1 - baryo['ratio']) * 100
        print(f"  η_B = dim(G₂_fund) × J_lep × exp(−2π²)")
        print(f"      = {baryo['dim_G2_fund']} × {baryo['J_lep']:.6f} × {baryo['tunnelling']:.4e}")
        print(f"      = {baryo['eta_B']:.4e}  (Planck: {baryo['eta_B_obs']:.2e}, {pull_pct:.2f}%)")

    if proofs is not None:
        S("Proof certificates")
        counts = proofs.get('status_counts', {})
        proved = counts.get('PROVED', 0)
        standard = counts.get('STANDARD_QFT', 0) + counts.get('STANDARD_RGE', 0)
        print(f"  Local mathematical certificates: {proved}")
        print(f"  Standard QFT/RGE certificates:   {standard}")

    # ── Emergence tree ───────────────────────────────────────────────

    m_c = quarks['m_c']
    m_t_em = quarks['m_t']  # emergence IS the primary top mass now
    m_b = quarks['m_b']
    m_s = quarks['m_s']
    p12 = pmns['pulls']['12']
    p23 = pmns['pulls']['23']
    p13 = pmns['pulls']['13']
    a_MZ = alpha_s['alpha_s_MZ_thresh']
    a_MZ_no = alpha_s['alpha_s_MZ_no_thresh']
    e_th = alpha_s['err_thresh']
    e_no = alpha_s['err_no_thresh']
    M_V = alpha_s['M_V_derived']

    H("THE EMERGENCE TREE")

    tree_lines = [
        "E₈(1) → G₂(1) × F₄(1)    [conformal embedding, c=8]",
        "  │",
        "  ├── G₂ sector → leptons",
        "  │     B/A = √2, θ = 2/9, Q = 2/3",
        "  │     (e, μ, τ): 3 masses",
        "  │",
        "  ├── F₄ sector → quarks + v_EW + Higgs",
        "  │     v_EW = M_Pl · exp(−(9π²/2 − 6 + 15/512))",
        "  │     30 modes = 26 (F₄ fund) + 4 (Higgs DOFs)",
        "  │     α(emergence) = π/512,  α(0) = 256(2π−1)/π² = 1/137.04",
        "  │     m_u = (38/9)m_e,  m_d = (83/9)m_e  (triality + h₁₀)",
        "  │",
        "  │     ── Emergence pattern ──",
        f"  │     Gen 1: m_u = (4+2/9) m_e",
        f"  │     Gen 2: m_c = (217/18) m_μ = {m_c:.1f} MeV",
        f"  │     Gen 3: m_t = (1165/12) m_τ = {m_t_em/1e3:.2f} GeV",
        "  │",
        f"  │     Yukawa cross-check: y_t = m_t/(v_EW/√2) ≈ 1 − 1/128",
        f"  │     Q = 2/3 + h₁₁/K³ = 289/432 → m_b = {m_b:.1f} MeV",
        f"  │     m_s = {m_s:.1f} MeV  (Rivero + √(32/27) bridge)",
        f"  │     λ(M_Pl)=−δ_bridge → m_H = {higgs['mH_pred']:.1f} GeV  (F₄ fusion + bridge threshold)",
        "  │",
        "  ├── SU(3)₃ / D⁽⁶⁾ → CKM",
        "  │     λ = tan(2/9),  A = √(2/3),  η̄ = π/9,  ρ̄ = √2/9",
        "  │",
        "  ├── SU(3)₃ / Z_C → PMNS",
        f"  │     sin²θ₁₂ = {pmns['sin2_12']:.5f} ({p12:+.1f}σ)",
        f"  │     sin²θ₂₃ = {pmns['sin2_23']:.5f} ({p23:+.1f}σ)",
        f"  │     sin²θ₁₃ = {pmns['sin2_13']:.5f} ({p13:+.1f}σ)",
        f"  │     δ_CP = {pmns['delta_CP_deg']:.1f}°",
        "  │",
    ]

    if nu is not None:
        tree_lines.extend([
            "  ├── F₄ singlets → neutrino masses",
            f"  │     26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ 1 ⊕ 1  ({nu['n_RH']} RH neutrinos)",
            f"  │     Rank-2 seesaw → m₁ = 0  (normal ordering)",
            f"  │     Σmν = {nu['sum_m_nu_meV']:.1f} meV  (testable by DESI/Euclid)",
            "  │",
        ])

    tree_lines.extend([
        "  ├── Embedding index → α_s",
        f"  │     α_G₂^WZW(v_EW) = π/32",
        f"  │     M_V = {M_V:.0f} GeV  (derived)",
        f"  │     α_s(M_Z) = {a_MZ:.4f}  ({e_th:+.1f}%)",
        "  │",
    ])

    if baryo is not None:
        b_pct = abs(1 - baryo['ratio']) * 100
        tree_lines.extend([
            "  ├── G₂ instanton + Fano orientation → η_B",
            f"  │     η_B = 7 × J_lep × exp(−2π²) = {baryo['eta_B']:.4e}  ({b_pct:.1f}%)",
            "  │",
        ])

    # Add gravity branch if available
    if grav is not None:
        tree_lines.extend([
            "  └── (7,26) bridge → gravity",
            f"        182 modes, ξ = 1/(48π),  h_bridge = 1",
            f"        PvP = 0 (universality),  Pv²P = ½P (memory)",
            f"        G_ind/G_N = {grav['G_ratio_UV']:.4f} (UV) – {grav['G_ratio_broken']:.4f} (broken)",
        ])
    else:
        # Close the tree without gravity
        tree_lines[-2] = tree_lines[-2].replace("├──", "└──")
        tree_lines[-1] = tree_lines[-1].replace("  │     ", "        ")

    tree_lines.extend([
        "",
        f"TOTAL: 9 masses + {'3 ν + ' if nu else ''}15 CKM + 3 PMNS + 1 α_s + 1 Higgs + {'1 G_N + ' if grav else ''}{'1 η_B' if baryo else ''}",
        f"     = {n_checks} numerical checks from M_Pl alone, zero free parameters",
        f"Mass: {n_sub1}/9 ≤ 1%, {n_sub2}/9 ≤ 2%, max {max_err:.1f}%",
        f"CKM:  χ²/n = {chi2_ckm/n_ckm:.2f}  ({n_ckm} obs)",
        f"PMNS: χ²/n = {chi2_pmns/3:.2f}  (3 obs)",
    ])

    if grav is not None:
        tree_lines.append(f"G_N:  G_ind/G_N = {grav['G_ratio_mid']:.4f}  (0.6–0.7%)")

    box(tree_lines)

    print(f"\n  {'═'*82}")
    if grav is not None:
        print(f"  END-TO-END DERIVATION COMPLETE -- ALL 248 E₈ MODES ASSIGNED")
        print(f"  (14,1) gauge + (1,52) matter/Higgs + (7,26) gravity = 248")
    else:
        print(f"  END-TO-END DERIVATION COMPLETE")
    print(f"  {'═'*82}\n")
