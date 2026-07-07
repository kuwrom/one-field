"""
The registry: every live bet, registered prediction, derivation
program, promotion, and rule, as code.

The registry lives in code so that the ledger cannot drift from the
engine: tests/test_registry.py asserts every registered value against
the solved web.  A registered number that stops matching canon fails
CI. The ledger is machine-checked, not remembered.

Registration discipline (unchanged):
  - A watch is registered before the deciding data exists, with its
    kill condition stated at registration time.
  - Watches are never silently edited; a supersession keeps the
    superseded value in its note, and git history keeps the rest.
  - DEPTH PRE-REGISTRATION RULE: before any new echo depth is
    evaluated in any channel, the complete grammar-permitted
    edge menu at that depth, with multiplicities, is committed first.
    A term identified after a residual is known, off-menu, is
    inadmissible.  First use: the m_H closure (see PROMOTIONS).

Render the ledger:  python interference/registry.py
"""

from fractions import Fraction

# ═══════════════════════════════════════════════════════════════════════
#  WATCHES: live bets ahead of the data, each with a kill condition
# ═══════════════════════════════════════════════════════════════════════

WATCHES = [
    dict(id=1, name="4/13 absence (alpha_s)",
         value=0.1184,
         bet="alpha_s(M_Z) = 0.1184 with NO anomalous gauge dressing; "
             "the G2-absorbed share of the bridge self-echo is invisible "
             "to gauge observables (no-self-dilution symmetry)",
         killed_by="alpha_s(M_Z) = 0.1177 or 0.1191 (a ±(4/13)/(2π) shift)",
         decided_by="lattice alpha_s at ~0.03%",
         note="the 4/13 back-reaction unit (0.0490) is nearly degenerate "
              "with the matching-scale unit (32/3)(15/512)/(2pi) = 0.0497 "
              "(1.4% apart); the scale is fixed first "
              "(mu* = M_Pl e^{-(9pi^2/2-6)}), then the dressing is tested "
              "on top; decidable at ~0.03% alpha_s"),
    dict(id=2, name="alpha(0) cesium commitment",
         value=137.035999050,
         bet="1/alpha(0) = 137.035999050 (depth-3 cubic, FINAL: no deeper "
             "term exists), +0.13 sigma of Berkeley Cs",
         killed_by="the Cs/Rb discrepancy (5.5 sigma) resolving on the Rb "
                   "side",
         decided_by="next-generation recoil measurement"),
    dict(id=3, name="m_s/m_ud ratio",
         value=27.318,
         bet="m_s/m_ud = 27.318 (RG-invariant, no coordinate)",
         killed_by="lattice average tightening at 27.30(8) and staying",
         decided_by="FLAG-class lattice ratio",
         note="supersedes 27.130"),
    dict(id=4, name="sin2(theta_W) emission weight",
         value=0.231285,
         bet="sin2(th_W)(M_Z) = 3/13 + Q0^2 alpha/(2pi) = 0.231285",
         killed_by="the 5/13 weight (0.231216) or any other family weight",
         decided_by="FCC-ee Z-pole (sigma ~ 9e-6)",
         note="supersedes 0.231234"),
    dict(id=5, name="m_b down-type vent",
         value=4193.8,
         bet="m_b = 4193.8 MeV with NO vent; if a vent exists it is "
             "±(9/4)(alpha/2pi) with epsilon-orientation signs (+d,-s,+b)",
         killed_by="persistence at >3 sigma without the registered "
                   "vent value",
         decided_by="PDG m_b error shrinking"),
    dict(id=6, name="Omega_DM/Omega_b",
         value=5.2832,
         bet="2pi - 1 = 5.2832 (sits between Planck -1.3 sigma and "
             "ACT +1.2 sigma; exact, no corrections possible: c_coset=0)",
         killed_by="both datasets converging away from 5.283",
         decided_by="CMB-S4 class"),
    dict(id=7, name="neutrino structural set",
         value=0.0,
         bet="m_nu1 = 0 exactly (rank-2 seesaw), normal ordering",
         killed_by="inverted ordering, or a measured m_nu1 > 0",
         decided_by="JUNO / cosmology"),
]

# ═══════════════════════════════════════════════════════════════════════
#  PREDICTIONS: registered exposures with kill conditions
# ═══════════════════════════════════════════════════════════════════════

PREDICTIONS = [
    dict(id="A", name="leptonic CP phase", value=76.9,
         claim="delta_CP = 76.9 deg",
         decided_by="DUNE / Hyper-K",
         kills="measured value incompatible with 76.9",
         note="STATUS: NuFit 6.0 NO global fit centers near "
              "177 deg and disfavors the first quadrant; the tension is "
              "present-tense, recorded so the exposure is not understated"),
    dict(id="B", name="fine structure", value=137.035999050,
         claim="1/alpha(0) = 137.035999050, cesium side, cubic FINAL",
         decided_by="Cs/Rb recoil dispute resolution",
         kills="dispute resolves toward Rb"),
    dict(id="C", name="weak angle form", value=0.2312854,
         claim="sin2 = 3/13 + Q0^2 alpha/(2pi) (depth-1 h7 + depth-2 "
               "h7/d11^2; supersedes 0.2312338) vs 5/13 "
               "alternative 0.2312161",
         decided_by="FCC-ee Z-pole (sigma ~ 1e-5)",
         kills="data selects the 5/13 form or neither"),
    dict(id="D", name="neutrino structure", value=0.0,
         claim="m1 = 0 exactly, normal ordering",
         decided_by="oscillation + cosmology",
         kills="inverted ordering, or m1 > 0 established"),
    dict(id="E", name="dark/baryon ratio", value=5.2832,
         claim="Omega_DM/Omega_b = 2pi - 1, exact",
         decided_by="CMB (Planck 5.36 ± ~0.05, ACT 5.22)",
         kills="converged CMB value excludes 5.2832"),
    dict(id="F", name="precision block",
         value=dict(G=6.674003e-11, m_mu=105.6583758, m_tau=1776.9092814,
                    v_EW=246.219645, M_Pl_GeV=1.2209171e19),
         claim="G (-2.0 sigma CODATA); m_mu (ppb, +0.12 sigma); m_tau "
               "(-0.23 sigma); v_EW vs G_F (-0.1 sigma at 1e-7); M_Pl "
               "(output)",
         decided_by="improved G; CODATA/PDG updates; Belle II (m_tau)",
         kills="improved G confirms CODATA central at <= 1e-5; any "
               "precision value moves outside its stated pull"),
    dict(id="G", name="RG-invariant light ratios",
         value=dict(mu_over_md=Fraction(38, 83), Q_ellipse=22.383,
                    ms_over_mud=27.318),
         claim="m_u/m_d = 38/83 = 0.45783 (-0.9 sigma); Q_ellipse = "
               "22.383 (+0.4 sigma dispersive, -1.7 sigma lattice, "
               "determinations disagree); m_s/m_ud carried as watch 3",
         decided_by="FLAG-class lattice ratios",
         kills="a converged ratio excludes the algebraic value"),
    dict(id=1, name="sum of neutrino masses", value=58.78,
         claim="Sigma m_nu = 58.78 meV (m1=0 + normal ordering + NuFit "
               "6.0 splittings dm21=7.49e-5, dm31=2.513e-3 eV^2, "
               "data-side inputs): m2 = 8.65, m3 = 50.13 meV",
         decided_by="DESI + CMB bound (~70 meV and tightening; the most "
                    "immediate exposure)",
         kills="a cosmological bound below ~58 meV, OR a "
               "detection above ~60 meV"),
    dict(id=2, name="KATRIN effective mass", value=8.81,
         claim="m_beta = 8.81 meV; direct kinematics sees nothing until "
               "sensitivity reaches ~9 meV",
         decided_by="KATRIN / Project 8",
         kills="any direct-kinematics detection above ~9 meV"),
    dict(id=3, name="neutrinoless double-beta", value=(1.51, 3.70),
         claim="m_bb in [1.51, 3.70] meV (Majorana phases only); every "
               "planned 0vbb experiment returns null",
         decided_by="LEGEND-1000, nEXO",
         kills="any 0vbb observation implying m_bb > 4 meV"),
    dict(id=4, name="dark matter is not a particle", value=None,
         claim="the dark sector is bridge coherence, coupling only "
               "through the common mode (PvP = 0); all direct/indirect/"
               "collider DM searches return null, indefinitely",
         decided_by="direct, indirect, collider searches",
         kills="any confirmed non-gravitational dark-matter detection"),
    dict(id=5, name="Fibonacci anyons at nu = 12/5", value=12/5,
         claim="k = 3 forces Read-Rezayi Z3 at nu = 2 + 2/(k+2) = 12/5: "
               "Fibonacci anyons; window: GaAs, B = 8.5 T, "
               "n = 4.9e11 cm^-2, T < 10 mK, mobility > 1e7, gap ~ 100 mK",
         decided_by="interferometry / thermal Hall at nu = 12/5",
         kills="established Abelian or non-Fibonacci order at 12/5 under "
               "the stated conditions"),
    dict(id=6, name="no fourth generation, ever", value=3,
         claim="exactly three windings close; a fourth is undefined, not "
               "suppressed",
         decided_by="colliders",
         kills="discovery of any fourth-generation fermion"),
    dict(id=7, name="equivalence principle exact", value=0.0,
         claim="PvP = 0 is an algebraic identity: null at EVERY future "
               "precision (MICROSCOPE eta < 1e-15 consistent)",
         decided_by="EP tests at any sensitivity",
         kills="any confirmed equivalence-principle violation"),
    dict(id=8, name="strong CP: theta_QCD = 0 exactly", value=0.0,
         claim="pi3(G2) = Z; the G2 -> SU(3)_c embedding has index 1 "
               "(instanton number preserved); E8(1) has a unique "
               "integrable vacuum module -> theta_G2 = 0; no axion, no "
               "new symmetry",
         decided_by="neutron EDM",
         kills="a measured neutron EDM attributable to theta_QCD != 0"),
    dict(id=9, name="Bell at exactly the Tsirelson bound", value=2**1.5,
         claim="local and NON-REALIST: the rejected premise is ±1 value "
               "assignment (Kochen-Specker), not locality; CHSH = "
               "2 sqrt(2) for ideal singlets at every distance, every "
               "knot species; no attenuation mechanism exists "
               "(tests/probes/bell.py). OPEN: derive E(a,b) = -a.b from "
               "two knots of one closure (Born conjecture -> theorem)",
         decided_by="loophole-free Bell tests",
         kills="any confirmed S > 2 sqrt(2), or a sub-quantum cap"),
]

SCOREBOARD_LOGIC = (
    "These sixteen entries (A-G, 1-9) are independent exposures. The "
    "framework survives only if ALL hold. A single clean failure (Rb "
    "wins, DUNE excludes 76.9, cosmology closes below 58 meV, a WIMP "
    "shows up, 0vbb at 10 meV) ends it, by its own rules. No new "
    "formulas were constructed for unmeasured constants by searching "
    "combinations of (2, 3, 7, 26, pi): in a four-integer grammar such "
    "searches always 'succeed', which is precisely why their successes "
    "would carry no weight. Every entry uses only edges already in the "
    "graph.")

# ═══════════════════════════════════════════════════════════════════════
#  PROMOTIONS: edges admitted into canon, with full provenance
# ═══════════════════════════════════════════════════════════════════════

PROMOTIONS = [dict(
    edge="lambda(M_Pl) = -N_bridge * alpha_G2^2 * E[v^2] * (1 - h10)",
    result="m_H = 125.2965 GeV (+0.88 sigma of 125.20(11))",
    lambda_MPl=-0.012450153,
    mH=125.2965,
    provenance=dict(
        menu="21 dictionary terms committed BEFORE evaluation (WZW "
             "weights h10/h11/h7/h26, Koide Q0/Q0^2/delta, altitudes "
             "1/K, 1/2K, h/K, dims 1/4, 1/9, face-splits 4/13, 9/13, "
             "mode-vents n/182, EM echoes), decision band fixed in "
             "advance at 125.20 ± 0.11 GeV",
        landing="EXACTLY ONE term landed: (1 - h10). Alternatives: "
                "1/4 -> +2.4 sigma, 1/6 -> -2.0, 26/182 -> -3.3, "
                "4/13 -> +5.4. RK4 certified at 2x steps.",
        scheme_fork="a scheme correction Delta shifts the menu common-"
                    "mode; landing windows: 2/9 iff Delta in "
                    "[-0.21, +0.01] GeV, 1/4 iff [-0.37, -0.15], "
                    "4/13 iff [-0.70, -0.48]. Raw MS-bar layer-3 import "
                    "(Delta = -0.40): NOTHING lands and canon degrades "
                    "(124.0 -> 123.6), rejected on the framework's own "
                    "scheme grounds: it mixes conventions and lands nothing.",
        convention="forced by the mass-coordinate doctrine: EW readings "
                   "at the derived M_Z, "
                   "tree relations; the extraction-scale test (m_H "
                   "moves -7.5 GeV if read at m_t) shows the doctrine "
                   "is load-bearing and single-valued, the same rule "
                   "every mass prediction uses",
        universal_residual="the scheme-universal layer-3 piece computed "
                           "exactly: beta2(g3, nf=6) = -65/2 shifts m_H "
                           "by +0.014 GeV, inside the window",
        truncation="layer-2 back-reaction is worth +7.9 GeV; quartic "
                   "layers self-certify at <= 35 MeV; iteration and "
                   "RK4 certified converged (< 1 MeV): layers, not "
                   "runtime, are the resolution limit"),
    open_theorem="WHY h10 is the vented share (derivation program 3)",
    kills="a derivation forcing a different factor; an improved m_H "
          "measurement excluding 125.30 at > 2 sigma; a derived "
          "in-convention layer-3 residual outside [-0.21, +0.01] GeV",
    isolation="sin2W, M_W, M_Z, G, v_EW unchanged to all displayed "
              "digits",
    frozen_in="tests/test_coverage.py::"
              "test_higgs_vent_closure",
)]

# ═══════════════════════════════════════════════════════════════════════
#  DERIVATION PROGRAMS: pre-committed proof obligations
#  (components named before the proofs are written; witnesses frozen
#   in tests/test_coverage.py; the proof cannot redefine its target)
# ═══════════════════════════════════════════════════════════════════════

DERIVATION_PROGRAMS = [
    dict(id=1, name="first-invariant-order theorem (why K^3)",
         claim="the Koide insertion is the Z3 sector-changer; the "
               "confined sector admits only Z3-invariant operators "
               "([Delta,S] = 0, selector theorem); tr S = tr S^2 = 0, "
               "tr S^3 = 3: orders 1-2 are FORBIDDEN, not small; first "
               "admissible order is cubic. COROLLARY: paired S S-dagger "
               "is invariant at order 2 -> bridge^2 gets K^2. One "
               "theorem, both altitudes.",
         witness="tests/test_coverage.py::test_first_invariant_order_theorem",
         obligation="promote the trace computation to the operator "
                    "statement in the MTC",
         kills="an admissible Z3-invariant first- or second-order "
               "unpaired insertion in the SU(3)_3 primary table"),
    dict(id=2, name="conversion lemma (integer ranks vs amplitudes)",
         claim="quantum dims are integers (1,2,3); the only "
               "irrationality source in the D6 walk algebra is the "
               "conversion vertex |C_3bar| = sqrt(2) (octonions.py, "
               "Fano-plane CG); every mixed terminus word carries an "
               "ODD number of channel switches -> amplitude = "
               "sqrt(2)^odd x rational = irrational -> cannot enter an "
               "integer rank. Pure towers: zero switches, admissible.",
         witness="tests/test_coverage.py::"
                 "test_conversion_vertex_irrationality_witness",
         obligation="parity bookkeeping: every switch inserts exactly "
                    "one C_3bar factor; no even-power completion inside "
                    "a length-4 word",
         kills="a mixed word whose full 6j evaluation is exactly "
               "rational"),
    dict(id=3, name="vent-share theorem (why h10, K^0, the scalar)",
         claim="(i) N-ality superselection (words.py, executable): a "
               "confining interface transmits only zero-triality; the "
               "7 of G2 decomposes 3+3bar+1, trialities (1,2,0): the "
               "fundamental share cannot reach the condensing channel "
               "[why the scalar]. (ii) emission-share grammar (sin2 "
               "th_W precedent): the blocked share in the weight metric "
               "is h10 [why the fundamental weight]. (iii) vents are "
               "selection rules (Casimir-vent precedent, kind='vent'): "
               "no altitude factor [why K^0].",
         witness="tests/test_coverage.py::test_vent_share_witnesses",
         obligation="the single-count rule: 3 and 3bar are two "
                    "orientations of one cycle; emission-share reads "
                    "the weight ONCE (sin2 precedent), the orientation "
                    "rule counts cycles TWICE (B = 2rho precedent); "
                    "derive which reading applies to a blocked share",
         kills="if the derivation forces the double count, the vent is "
               "(1 - 4/9) -> m_H ~ 126.6 GeV, excluded at > 2 sigma "
               "(frozen in the test suite): the wrong answer is "
               "already dead"),
]

# ═══════════════════════════════════════════════════════════════════════
#  OUT-OF-DICTIONARY CONSISTENCY CHECKS (data existed, never consulted)
# ═══════════════════════════════════════════════════════════════════════

CONSISTENCY_CHECKS = [
    dict(name="joint cosmology closure",
         claim="eta_B (G2 instanton) and Omega_DM/Omega_b (bridge "
               "venting) are independent; jointly: Omega_b h^2 = "
               "eta_10/273.9 = 0.02255 (+1.2 sigma) and x(2pi-1) -> "
               "Omega_DM h^2 = 0.1192 (-0.7 sigma vs Planck 0.1200(12))",
         kills="Planck-class data pinning both densities at values "
               "whose ratio and product cannot be met by "
               "(eta_B, 2pi - 1)"),
    dict(name="neutron-proton splitting, QCD part",
         claim="m_d - m_u = 5 m_e = 2.5550 MeV exactly; lattice (BMW "
               "2015): 2.52 ± 0.23 -> +0.15 sigma",
         kills="lattice QCD part excluding 5 m_e"),
    dict(name="constituent-quark scale (scale-level only)",
         claim="Lambda_conf = 314.2 MeV lands at the constituent scale "
               "(~1%); the 0.47% closeness to m_p/3 specifically is "
               "partly fortuitous and NOT claimed. FORCED qualitative: "
               "SU(3)_3 fusion gives h(3bar) = 2/9 < h(6) = 5/9 -> "
               "m_N < m_Delta (observed 939 < 1232). Once the composite "
               "Hamiltonian is derived, one measurement fixes the "
               "conversion and the baryon sector is over-constrained.",
         kills="derived composite spectrum contradicting the ordering"),
]

# ═══════════════════════════════════════════════════════════════════════
#  CLOSURE RECORD: the four structural questions, asked and closed
# ═══════════════════════════════════════════════════════════════════════

CLOSURE_RECORD = [
    dict(name="chirality",
         how="tests/probes/zero_mode.py: winding ±1 carry one true zero "
             "mode each in opposite chirality blocks; windingless "
             "gapped. Handedness = winding sign = the propagated Fano "
             "bit. V-A inherited, not imposed.",
         corroboration="Callan-Harvey 1985, Jackiw-Rossi 1981"),
    dict(name="fermionic statistics",
         how="tests/probes/bookkeeping.py + winding_texture.py: WZ "
             "coefficient = 3 microscopic Dirac species (the same count "
             "alpha_s uses); level 3 odd -> fermions; |B| = 0.91 ~ 1 "
             "with pi_3 protection. Generations and statistics share "
             "one origin: k = 3.",
         corroboration="Witten 1983; 't Hooft anomaly matching"),
    dict(name="saturation -> Koide bridge",
         how="classical 1D arc (knot/stationary/spectrum/universality "
             "probes): matter forms because of the G2 coupling; the "
             "circulant fringe is classically ABSENT; superselection "
             "places the three-line spectrum at the MTC layer "
             "(mtc_spectrum.py: Q = 2/3, mu/e = 206.77; theta_menu.py: "
             "theta = 2/9 unique).",
         corroboration="Cardy boundary-condition-changing operators"),
    dict(name="the selector",
         how="no selector, no wall: confinement = unbroken Z3 centre "
             "('t Hooft); admissible mass operators satisfy "
             "[Delta,S] = 0; diag(v) inadmissible by symmetry. "
             "Executable: tests/probes/selector_theorem.py (idempotent "
             "measures on Z3 = {Haar, delta_e}, the dichotomy).",
         corroboration="centre symmetry, standard confinement lore"),
]

# ═══════════════════════════════════════════════════════════════════════
#  EPISTEMIC RECORD: negative results and ordering, kept verbatim
# ═══════════════════════════════════════════════════════════════════════

EPISTEMIC_RECORD = [
    "The dynamical origin of B/A = sqrt(2) via instanton tunnelling "
    "was tested and REJECTED by rigorous computation.",
    "theta = 2/9 and r = sqrt(2) were first read from PDG masses; "
    "octonions.py derives both constructively.",
    "The 1D endpoint is a seed-dependent turbulent breather; standing "
    "character is universal, uniqueness of the attractor is not "
    "established in 1D (tests/probes/universality.py).",
    "Lines below omega = 0.11 were beneath the spectroscopy window's "
    "resolution (tests/probes/spectrum.py).",
]


# ═══════════════════════════════════════════════════════════════════════
#  PAPERS
# ═══════════════════════════════════════════════════════════════════════

PAPERS = [
    ("The Innocent Lepton", "10.5281/zenodo.19899091"),
    ("One Substrate, Three Generations", "10.5281/zenodo.20069456"),
    ("The Echo of Standing Waves", "10.5281/zenodo.20144381"),
    ("The Octavian Coherence Gate: Rigidity of the Four Irreducible "
     "Integers of the E8(1) > G2(1) x F4(1) Conformal Embedding",
     "10.5281/zenodo.20493955"),
]
CITE = ("Kahsay, Kibrom Kidane (2026). one-field. "
        "https://github.com/kuwrom/one-field (the repository "
        "supersedes the papers as the living canon).")


def render():
    """Human-readable ledger dump."""
    out = ["=" * 72, "  THE REGISTRY (machine-checked ledger)", "=" * 72]
    out.append("\n-- WATCHES (live bets, kill conditions) --")
    for w in WATCHES:
        out.append(f"  [{w['id']}] {w['name']}: {w['bet'][:64]}...")
        out.append(f"      killed by: {w['killed_by']}")
    out.append("\n-- PREDICTIONS (registered exposures) --")
    for p in PREDICTIONS:
        out.append(f"  [{p['id']}] {p['name']}: kills = {p['kills'][:60]}")
    out.append("\n-- PROMOTIONS --")
    for pr in PROMOTIONS:
        out.append(f"  {pr['edge']}")
        out.append(f"      -> {pr['result']}; open: {pr['open_theorem']}")
    out.append("\n-- DERIVATION PROGRAMS --")
    for d in DERIVATION_PROGRAMS:
        out.append(f"  [{d['id']}] {d['name']}")
        out.append(f"      witness: {d['witness']}")
        out.append(f"      kills:   {d['kills'][:64]}")
    out.append("\n-- SCOREBOARD --")
    out.append("  " + SCOREBOARD_LOGIC[:200] + "...")
    return "\n".join(out)


if __name__ == "__main__":
    print(render())
