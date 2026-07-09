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
  - Watches are never silently edited. A supersession keeps the
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
#  AXIOMS: the stated inputs, each with its layer and status
#  (everything else in the ledger is derived, cited as corroboration,
#   or registered as an open obligation. This section is the complete
#   list of what is ASSUMED)
# ═══════════════════════════════════════════════════════════════════════

AXIOMS = [
    dict(id="A1", name="self-presence (interference)",
         statement="to be is, minimally, to be self-present. What "
                   "exists carries non-zero interference weight",
         layer="founding", status="axiom (stated first, README)"),
    dict(id="A2", name="zero free parameters (R1)",
         statement="the closure must exist with zero dimensionless "
                   "free parameters. k = 3 (R3) is a LEMMA of this "
                   "demand, not a separate axiom "
                   "(tests/probes/wiring_scan.py)",
         layer="web", status="axiom. Its consequences are the "
                             "scanned closure conditions"),
    dict(id="A3", name="3+1 dimensions",
         statement="three spatial dimensions plus time for the "
                   "substrate (pi_3 classification and codim-4 node "
                   "events read off it)",
         layer="substrate", status="assumption (prose-motivated, "
                                   "not derived)"),
    dict(id="A4", name="NLS kinetic form + locality",
         statement="the substrate is a local Z3-coupled NLS medium "
                   "with quadratic-gradient kinetics "
                   "(interference/nls_soliton.py Eq. 2.1)",
         layer="substrate", status="assumption. Binds only the "
                                   "conjecture layer "
                                   "(SUBSTRATE_CONJECTURES)"),
    dict(id="A5", name="electron anchor",
         statement="m_e (CODATA) fixes the one dimensional scale. "
                   "H0 and the oscillation splittings enter the "
                   "Lambda/neutrino sectors as measured references",
         layer="web", status="dimensional input, not a coupling"),
]

# ═══════════════════════════════════════════════════════════════════════
#  WATCHES: live bets ahead of the data, each with a kill condition
# ═══════════════════════════════════════════════════════════════════════

WATCHES = [
    dict(id=1, name="4/13 absence (alpha_s)",
         value=0.1184,
         bet="alpha_s(M_Z) = 0.1184 with NO anomalous gauge dressing. "
             "the G2-absorbed share of the bridge self-echo is invisible "
             "to gauge observables (no-self-dilution symmetry)",
         killed_by="alpha_s(M_Z) = 0.1177 or 0.1191 (a ±(4/13)/(2π) shift)",
         decided_by="lattice alpha_s at ~0.03%",
         note="the 4/13 back-reaction unit (0.0490) is nearly degenerate "
              "with the matching-scale unit (32/3)(15/512)/(2pi) = 0.0497 "
              "(1.4% apart). The scale is fixed first "
              "(mu* = M_Pl e^{-(9pi^2/2-6)}), then the dressing is tested "
              "on top. Decidable at ~0.03% alpha_s"),
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
         bet="m_b = 4193.8 MeV with NO vent. If a vent exists it is "
             "±(9/4)(alpha/2pi) with epsilon-orientation signs (+d,-s,+b)",
         killed_by="persistence at >3 sigma without the registered "
                   "vent value",
         decided_by="PDG m_b error shrinking"),
    dict(id=6, name="Omega_DM/Omega_b",
         value=5.2832,
         bet="2pi - 1 = 5.2832 (sits between Planck -1.3 sigma and "
             "ACT +1.2 sigma. Exact, no corrections possible: c_coset=0)",
         killed_by="both datasets converging away from 5.283",
         decided_by="CMB-S4 class"),
    dict(id=7, name="neutrino structural set",
         value=0.0,
         bet="m_nu1 = 0 exactly (rank-2 seesaw), normal ordering",
         killed_by="inverted ordering, or a measured m_nu1 > 0",
         decided_by="JUNO / cosmology"),
    dict(id=8, name="gravitational-wave speed",
         value=0.0,
         bet="|c_gw/c - 1| = 0 identically: tensor modes and light "
             "ride one emergent metric (gravity.py is web-layer, so "
             "this binds regardless of the architecture reading)",
         killed_by="any confirmed |c_gw/c - 1| > 1e-15",
         decided_by="multi-messenger events (GW170817-class bound "
                    "~1e-15 already consistent)"),
    dict(id=9, name="PPN gamma",
         value=1.0,
         bet="gamma = 1 exactly: no scalar admixture in the emergent "
             "metric read by light bending and time delay",
         killed_by="a confirmed |gamma - 1| > 2e-5",
         decided_by="Cassini-class time delay and successors "
                    "(current bound 2.3e-5)"),
    # The third gravity-sector exposure, the equivalence principle,
    # is PREDICTIONS #7 (PvP = 0, kill condition stated there). It is
    # not duplicated here.
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
              "177 deg and disfavors the first quadrant. The tension is "
              "present-tense, recorded so the exposure is not understated"),
    dict(id="B", name="fine structure", value=137.035999050,
         claim="1/alpha(0) = 137.035999050, cesium side, cubic FINAL",
         decided_by="Cs/Rb recoil dispute resolution",
         kills="dispute resolves toward Rb"),
    dict(id="C", name="weak angle form", value=0.2312854,
         claim="sin2 = 3/13 + Q0^2 alpha/(2pi) (depth-1 h7 + depth-2 "
               "h7/d11^2. Supersedes 0.2312338) vs 5/13 "
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
         claim="G (-2.0 sigma CODATA). m_mu (ppb, +0.12 sigma). m_tau "
               "(-0.23 sigma). v_EW vs G_F (-0.1 sigma at 1e-7). M_Pl "
               "(output)",
         decided_by="improved G. CODATA/PDG updates. Belle II (m_tau)",
         kills="improved G confirms CODATA central at <= 1e-5. Any "
               "precision value moves outside its stated pull"),
    dict(id="G", name="RG-invariant light ratios",
         value=dict(mu_over_md=Fraction(38, 83), Q_ellipse=22.383,
                    ms_over_mud=27.318),
         claim="m_u/m_d = 38/83 = 0.45783 (-0.9 sigma). Q_ellipse = "
               "22.383 (+0.4 sigma dispersive, -1.7 sigma lattice, "
               "determinations disagree). m_s/m_ud carried as watch 3",
         decided_by="FLAG-class lattice ratios",
         kills="a converged ratio excludes the algebraic value"),
    dict(id=1, name="sum of neutrino masses", value=58.78,
         claim="Sigma m_nu = 58.78 meV (m1=0 + normal ordering + NuFit "
               "6.0 splittings dm21=7.49e-5, dm31=2.513e-3 eV^2, "
               "data-side inputs): m2 = 8.65, m3 = 50.13 meV",
         decided_by="DESI + CMB bound (~70 meV and tightening. The most "
                    "immediate exposure)",
         kills="a cosmological bound below ~58 meV, OR a "
               "detection above ~60 meV"),
    dict(id=2, name="KATRIN effective mass", value=8.81,
         claim="m_beta = 8.81 meV. Direct kinematics sees nothing until "
               "sensitivity reaches ~9 meV",
         decided_by="KATRIN / Project 8",
         kills="any direct-kinematics detection above ~9 meV"),
    dict(id=3, name="neutrinoless double-beta", value=(1.51, 3.70),
         claim="m_bb in [1.51, 3.70] meV (Majorana phases only). Every "
               "planned 0vbb experiment returns null",
         decided_by="LEGEND-1000, nEXO",
         kills="any 0vbb observation implying m_bb > 4 meV"),
    dict(id=4, name="dark matter is not a particle", value=None,
         claim="the dark sector is bridge coherence, coupling only "
               "through the common mode (PvP = 0). All direct/indirect/"
               "collider DM searches return null, indefinitely",
         decided_by="direct, indirect, collider searches",
         kills="any confirmed non-gravitational dark-matter detection"),
    dict(id=5, name="Fibonacci anyons at nu = 12/5", value=12/5,
         claim="k = 3 forces Read-Rezayi Z3 at nu = 2 + 2/(k+2) = 12/5: "
               "Fibonacci anyons. Window: GaAs, B = 8.5 T, "
               "n = 4.9e11 cm^-2, T < 10 mK, mobility > 1e7, gap ~ 100 mK",
         decided_by="interferometry / thermal Hall at nu = 12/5",
         kills="established Abelian or non-Fibonacci order at 12/5 under "
               "the stated conditions"),
    dict(id=6, name="no fourth generation, ever", value=3,
         claim="exactly three windings close. A fourth is undefined, not "
               "suppressed",
         decided_by="colliders",
         kills="discovery of any fourth-generation fermion"),
    dict(id=7, name="equivalence principle exact", value=0.0,
         claim="PvP = 0 is an algebraic identity: null at EVERY future "
               "precision (MICROSCOPE eta < 1e-15 consistent)",
         decided_by="EP tests at any sensitivity",
         kills="any confirmed equivalence-principle violation"),
    dict(id=8, name="strong CP: theta_QCD = 0 exactly", value=0.0,
         claim="pi3(G2) = Z. The G2 -> SU(3)_c embedding has index 1 "
               "(instanton number preserved). E8(1) has a unique "
               "integrable vacuum module -> theta_G2 = 0. No axion, no "
               "new symmetry",
         decided_by="neutron EDM",
         kills="a measured neutron EDM attributable to theta_QCD != 0"),
    dict(id=9, name="Bell at exactly the Tsirelson bound", value=2**1.5,
         claim="local and NON-REALIST: the rejected premise is ±1 value "
               "assignment (Kochen-Specker), not locality. CHSH = "
               "2 sqrt(2) for ideal singlets at every distance, every "
               "knot species. No attenuation mechanism exists "
               "(tests/probes/bell.py). OPEN: derive E(a,b) = -a.b from "
               "two knots of one closure (Born conjecture -> theorem)",
         decided_by="loophole-free Bell tests",
         kills="any confirmed S > 2 sqrt(2), or a sub-quantum cap"),
    dict(id=10, name="proton stability (web bookkeeping)", value=None,
         claim="the conserved quantity is the fermion-line terminus: "
               "every edge in the web's taxonomy (ratio/echo/vent, "
               "frozen in root.EchoTerm) composes amplitudes AT a "
               "species node and no edge kind changes a knot's "
               "species, so baryon and lepton number are conserved "
               "node labels of the graph. A B-violating channel would "
               "require an identity-changing edge kind that does not "
               "exist and cannot be admitted without a new theorem "
               "under the depth pre-registration rule. Lowest "
               "admissible B-violating edge: NONE at any depth <= 3. "
               "Registered against tau_p > 1.6e34 yr (Super-K, "
               "consistent: the web predicts no decay at any "
               "sensitivity).",
         decided_by="proton-decay searches (Super-K, Hyper-K, DUNE)",
         kills="any confirmed proton decay (an identity-changing "
               "process the edge taxonomy cannot host)"),
]

SCOREBOARD_LOGIC = (
    "These seventeen entries (A-G, 1-10) are independent exposures. The "
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
                    "mode. Landing windows: 2/9 iff Delta in "
                    "[-0.21, +0.01] GeV, 1/4 iff [-0.37, -0.15], "
                    "4/13 iff [-0.70, -0.48]. Raw MS-bar layer-3 import "
                    "(Delta = -0.40): NOTHING lands and canon degrades "
                    "(124.0 -> 123.6), rejected on the framework's own "
                    "scheme grounds: it mixes conventions and lands nothing.",
        convention="forced by the mass-coordinate doctrine: EW readings "
                   "at the derived M_Z, "
                   "tree relations. The extraction-scale test (m_H "
                   "moves -7.5 GeV if read at m_t) shows the doctrine "
                   "is load-bearing and single-valued, the same rule "
                   "every mass prediction uses",
        universal_residual="the scheme-universal layer-3 piece computed "
                           "exactly: beta2(g3, nf=6) = -65/2 shifts m_H "
                           "by +0.014 GeV, inside the window",
        truncation="layer-2 back-reaction is worth +7.9 GeV. Quartic "
                   "layers self-certify at <= 35 MeV. Iteration and "
                   "RK4 certified converged (< 1 MeV): layers, not "
                   "runtime, are the resolution limit"),
    open_theorem="WHY h10 is the vented share (derivation program 3)",
    kills="a derivation forcing a different factor. An improved m_H "
          "measurement excluding 125.30 at > 2 sigma. A derived "
          "in-convention layer-3 residual outside [-0.21, +0.01] GeV",
    isolation="sin2W, M_W, M_Z, G, v_EW unchanged to all displayed "
              "digits",
    frozen_in="tests/test_coverage.py::"
              "test_higgs_vent_closure",
)]

# ═══════════════════════════════════════════════════════════════════════
#  DERIVATION PROGRAMS: pre-committed proof obligations
#  (components named before the proofs are written. Witnesses frozen
#   in tests/test_coverage.py. The proof cannot redefine its target)
#
#  GRADE (2026-07 adjudication): programs 1-3 select their structures
#  by menu uniqueness, the same grade as every accepted edge in the
#  web, and that standard is now statistically validated (null model
#  + blind holdout, test_coverage.py).  Their open halves are
#  symbolic strengthenings, not claim-blocking gaps.  Program 4 is
#  the substrate keystone and remains a genuine open front, together
#  with SUBSTRATE_CONJECTURES #1.
# ═══════════════════════════════════════════════════════════════════════

DERIVATION_PROGRAMS = [
    dict(id=1, name="first-invariant-order theorem (why K^3)",
         claim="the Koide insertion is the Z3 sector-changer. The "
               "confined sector admits only Z3-invariant operators "
               "([Delta,S] = 0, selector theorem). Tr S = tr S^2 = 0, "
               "tr S^3 = 3: orders 1-2 are FORBIDDEN, not small. First "
               "admissible order is cubic. COROLLARY: paired S S-dagger "
               "is invariant at order 2 -> bridge^2 gets K^2. One "
               "theorem, both altitudes.",
         witness="tests/test_coverage.py::test_first_invariant_order_theorem",
         obligation="promote the trace computation to the operator "
                    "statement in the MTC",
         kills="an admissible Z3-invariant first- or second-order "
               "unpaired insertion in the SU(3)_3 primary table"),
    dict(id=2, name="conversion lemma (integer ranks vs amplitudes)",
         claim="quantum dims are integers (1,2,3). The only "
               "irrationality source in the D6 walk algebra is the "
               "conversion vertex |C_3bar| = sqrt(2) (octonions.py, "
               "Fano-plane CG). Every mixed terminus word carries an "
               "ODD number of channel switches -> amplitude = "
               "sqrt(2)^odd x rational = irrational -> cannot enter an "
               "integer rank. Pure towers: zero switches, admissible.",
         witness="tests/test_coverage.py::"
                 "test_conversion_vertex_irrationality_witness",
         obligation="parity bookkeeping: every switch inserts exactly "
                    "one C_3bar factor. No even-power completion inside "
                    "a length-4 word",
         kills="a mixed word whose full 6j evaluation is exactly "
               "rational"),
    dict(id=3, name="vent-share theorem (why h10, K^0, the scalar)",
         claim="(i) N-ality superselection (words.py, executable): a "
               "confining interface transmits only zero-triality. The "
               "7 of G2 decomposes 3+3bar+1, trialities (1,2,0): the "
               "fundamental share cannot reach the condensing channel "
               "[why the scalar]. (ii) emission-share grammar (sin2 "
               "th_W precedent): the blocked share in the weight metric "
               "is h10 [why the fundamental weight]. (iii) vents are "
               "selection rules (Casimir-vent precedent, kind='vent'): "
               "no altitude factor [why K^0].",
         witness="tests/test_coverage.py::test_vent_share_witnesses",
         obligation="the single-count rule: 3 and 3bar are two "
                    "orientations of one cycle. Emission-share reads "
                    "the weight ONCE (sin2 precedent), the orientation "
                    "rule counts cycles TWICE (B = 2rho precedent). "
                    "derive which reading applies to a blocked share",
         kills="if the derivation forces the double count, the vent is "
               "(1 - 4/9) -> m_H ~ 126.6 GeV, excluded at > 2 sigma "
               "(frozen in the test suite): the wrong answer is "
               "already dead"),
    dict(id=4, name="texture stabilization in the confined phase",
         claim="BRANCH (i) CLOSED constructively (2026-07, "
               "tests/probes/effective_action.py): the fixed-norm "
               "sector's effective action contains the induced "
               "Faddeev-Skyrme quartic, coefficient rho/(8 m^2) with "
               "m the MEASURED charged-sector gap (2.749 = 0.74 mu, "
               "knot_spectrum.py), via the exact CP^1 split and the "
               "Mermin-Ho relation (both verified on the texture). "
               "Derrick minimum at R* = 1.0 healing lengths, local "
               "limit self-consistent (R* m = 1.4). The bare Z3-NLS "
               "does NOT protect pi_3 winding (measured: ballistic "
               "collapse, tau ~ R/c_s, exponent 1.29 over R in "
               "[1.8, 3.2] in the review runs. Unwinding through "
               "codim-4 amplitude-node events. Rotating textures "
               "have no conserved backing charge, commutation "
               "residual 8e-3 vs 5e-15 control), and the derivation "
               "explains why: the amplitude node EXITS the fixed-"
               "norm manifold, the one channel the quartic cannot "
               "see. Corroboration: Babaev-Faddeev-Niemi 2002.",
         witness="tests/test_coverage.py::"
                 "test_texture_stabilization_witness",
         obligation="the single remaining premise: the node channel "
                    "is suppressed in the saturated quantum regime "
                    "(= SUBSTRATE_CONJECTURES #1, sharpened there). "
                    "full records: tests/probes/knot_charge.py, "
                    "tests/probes/effective_action.py",
         kills="if the amplitude/node channel is not suppressed in "
               "the saturated quantum regime (the gap closes beyond "
               "Gaussian order, or the substrate diluteness "
               "rho xi^3 cannot be large), unit winding has no "
               "dynamical origin and the statistics/chirality chain "
               "breaks"),
]

# ═══════════════════════════════════════════════════════════════════════
#  OUT-OF-DICTIONARY CONSISTENCY CHECKS (data existed, never consulted)
# ═══════════════════════════════════════════════════════════════════════

CONSISTENCY_CHECKS = [
    dict(name="joint cosmology closure",
         claim="eta_B (G2 instanton) and Omega_DM/Omega_b (bridge "
               "venting) are independent. Jointly: Omega_b h^2 = "
               "eta_10/273.9 = 0.02255 (+1.2 sigma) and x(2pi-1) -> "
               "Omega_DM h^2 = 0.1192 (-0.7 sigma vs Planck 0.1200(12))",
         kills="Planck-class data pinning both densities at values "
               "whose ratio and product cannot be met by "
               "(eta_B, 2pi - 1)"),
    dict(name="neutron-proton splitting, QCD part",
         claim="m_d - m_u = 5 m_e = 2.5550 MeV exactly. Lattice (BMW "
               "2015): 2.52 ± 0.23 -> +0.15 sigma",
         kills="lattice QCD part excluding 5 m_e"),
    dict(name="constituent-quark scale (scale-level only)",
         claim="Lambda_conf = 314.2 MeV lands at the constituent scale "
               "(~1%). The 0.47% closeness to m_p/3 specifically is "
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
             "mode each in opposite chirality blocks. Windingless "
             "gapped. Handedness = winding sign = the propagated Fano "
             "bit. V-A inherited, not imposed. Labeled remainder: the "
             "Lorentz-chiral (V-A) form of the 4d vertex itself "
             "(tests/probes/orientation_bit.py).",
         corroboration="Callan-Harvey 1985, Jackiw-Rossi 1981"),
    dict(name="fermionic statistics",
         how="tests/probes/bookkeeping.py + winding_texture.py: WZ "
             "coefficient = 3 microscopic Dirac species (the same count "
             "alpha_s uses). Level 3 odd -> fermions. The static unit "
             "texture carries |B| -> 1 (estimator-converged). "
             "Protection of that winding is NOT substrate-level: the "
             "bare Z3-NLS unwinds it ballistically through codim-4 "
             "amplitude-node events, and no Noether charge backs "
             "internal rotation (tests/probes/knot_charge.py). "
             "Protection is an obligation of the confined-phase "
             "effective theory, where the level-3 WZ structure already "
             "lives (DERIVATION_PROGRAMS #4). Generations and "
             "statistics share one origin: k = 3.",
         corroboration="Witten 1983. 't Hooft anomaly matching"),
    dict(name="saturation -> Koide bridge",
         how="classical 1D arc (knot/stationary/spectrum/universality "
             "probes): matter forms because of the G2 coupling. The "
             "circulant fringe is classically ABSENT. Superselection "
             "places the three-line spectrum at the MTC layer "
             "(mtc_spectrum.py: Q = 2/3, mu/e = 206.77. theta_menu.py: "
             "theta = 2/9 unique).",
         corroboration="Cardy boundary-condition-changing operators"),
    dict(name="the selector",
         how="no selector, no wall: confinement = unbroken Z3 centre "
             "('t Hooft). Admissible mass operators satisfy "
             "[Delta,S] = 0. diag(v) inadmissible by symmetry. "
             "Executable: tests/probes/selector_theorem.py (idempotent "
             "measures on Z3 = {Haar, delta_e}, the dichotomy).",
         corroboration="centre symmetry, standard confinement lore"),
]

# ═══════════════════════════════════════════════════════════════════════
#  SUBSTRATE-CONJECTURE LEDGER: obligations of the physical layer
#
#  Two-layer architecture, stated as a rule: the recursive web is the
#  COMPUTATION (every scorecard number, watch, and prediction attaches
#  to it). The substrate is the PHYSICAL CONJECTURE.  The entries below
#  are what the substrate story owes before web-layer claims can be
#  read as statements about a physical medium.  The first is the
#  entry test: everything else waits on it.
# ═══════════════════════════════════════════════════════════════════════

SUBSTRATE_CONJECTURES = [
    dict(id=1, name="confined topological phase exists (ENTRY TEST)",
         claim="the Z3-coupled substrate at knot saturation enters a "
               "confined phase with unbroken Z3 centre: the premise "
               "under the WZ citation, winding quantization, and the "
               "substrate-web bridge. SHARPENED (2026-07) to one "
               "quantitative condition: the amplitude/node channel "
               "of the charged sector is suppressed in the saturated "
               "quantum regime. Everything else that protection "
               "needs is now derived or measured: the stabilizing "
               "quartic exists given the gap "
               "(effective_action.py), the gap is measured at "
               "Gaussian order (omega = 2.749 = 0.74 mu, "
               "knot_spectrum.py. Bogoliubov linear response = the "
               "quantum quasiparticle spectrum at that order), and "
               "the classical endpoint is defect-free "
               "(confined_phase.py). Beyond Gaussian order a node "
               "event is a finite-action excursion, suppressed as "
               "exp(-c rho xi^3). The diluteness rho xi^3 is a "
               "substrate-depth parameter the web does not fix. "
               "READ AS A PREDICTION: stable fermionic matter "
               "exists, so through the closed implication chain the "
               "framework requires the substrate to be a deep "
               "condensate (rho xi^3 >> 1). Any independent handle "
               "on substrate depth must land there.",
         status="OPEN on the single premise. Killed by: a derivation "
                "forcing rho xi^3 ~ 1, or the charged-sector gap "
                "closing beyond Gaussian order in the saturated "
                "regime",
         depends_on=None),
    dict(id=2, name="Lorentz universality",
         claim="one emergent light cone for all excitations of the "
               "substrate (the common mode's c_s is universal at long "
               "wavelength). Carries the 4d Lorentz-chiral (V-A) "
               "vertex form as a dependent item (mechanism cited: "
               "Callan-Harvey inflow. tests/probes/orientation_bit.py)",
         status="OPEN. Numeric kill: any confirmed photon-sector "
                "Lorentz violation (current bounds: linear-dispersion "
                "scale E_QG > 1e19 GeV from GRB time-of-flight. "
                "substrate corrections enter at (k xi_Pl)^2, far "
                "below all bounds if the healing length is Planckian)",
         depends_on=1),
    dict(id=3, name="gauge fields as substrate excitations",
         claim="the gauge sector of the web is carried by collective "
               "excitations of the confined phase, not postulated "
               "fields. Includes reproducing the perturbative gauge "
               "series on composite knots (g-2, scattering "
               "amplitudes, decay widths) and the multi-knot Fock/"
               "Pauli structure beyond two-body exchange, the "
               "quantum-observables half of the conjecture",
         status="OPEN",
         depends_on=1),
    dict(id=4, name="quantized knots (unit winding)",
         claim="unit pi_3 winding is created and protected by the "
               "confined-phase effective theory. The bare substrate "
               "does NOT protect it (tests/probes/knot_charge.py, "
               "measured). This is DERIVATION_PROGRAMS #4.",
         status="OPEN. Substrate-level alternatives exhausted",
         depends_on=1),
]

# ═══════════════════════════════════════════════════════════════════════
#  EPISTEMIC RECORD: negative results and ordering, kept verbatim
# ═══════════════════════════════════════════════════════════════════════

EPISTEMIC_RECORD = [
    "The dynamical origin of B/A = sqrt(2) via instanton tunnelling "
    "was tested and REJECTED by rigorous computation.",
    "theta = 2/9 and r = sqrt(2) were first read from PDG masses. "
    "octonions.py derives both constructively (postdicted, then "
    "derived. Discrete menu values, no continuous freedom).",
    "The 1D endpoint is a seed-dependent turbulent breather. Standing "
    "character is universal, uniqueness of the attractor is not "
    "established in 1D (tests/probes/universality.py).",
    "Lines below omega = 0.11 were beneath the spectroscopy window's "
    "resolution (tests/probes/spectrum.py).",
    "W = 1/d_lambda was carried as one of the gate paper's three "
    "stated conventions. The projector lemma (octonions.py, "
    "2026-07) upgraded it to forced: idempotency fixes the "
    "coefficient at 1/(loop value), and the loop value is the "
    "quantum dimension by definition of the quantum trace.",
]


# ═══════════════════════════════════════════════════════════════════════
#  PAPERS
# ═══════════════════════════════════════════════════════════════════════

PAPERS = [
    ("The Innocent Lepton", "10.5281/zenodo.19899091"),
    ("One Substrate, Three Generations", "10.5281/zenodo.20069456"),
    ("The Echo of Standing Waves", "10.5281/zenodo.20144381"),
    ("The Octavian Coherence Gate: The Four Irreducible Integers of "
     "the E8(1) > G2(1) x F4(1) Conformal Embedding",
     "10.5281/zenodo.20493955"),
]
CITE = ("Kahsay, Kibrom Kidane (2026). One-field. "
        "https://github.com/kuwrom/one-field (the repository "
        "supersedes the papers as the living canon).")


def render():
    """Human-readable ledger dump."""
    out = ["=" * 72, "  THE REGISTRY (machine-checked ledger)", "=" * 72]
    out.append("\n-- AXIOMS (the complete list of what is assumed) --")
    for a in AXIOMS:
        out.append(f"  [{a['id']}] {a['name']} ({a['layer']}): "
                   f"{a['status'][:52]}")
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
        out.append(f"      -> {pr['result']}. Open: {pr['open_theorem']}")
    out.append("\n-- DERIVATION PROGRAMS --")
    for d in DERIVATION_PROGRAMS:
        out.append(f"  [{d['id']}] {d['name']}")
        out.append(f"      witness: {d['witness']}")
        out.append(f"      kills:   {d['kills'][:64]}")
    out.append("\n-- SUBSTRATE CONJECTURES (physical-layer ledger) --")
    for s in SUBSTRATE_CONJECTURES:
        out.append(f"  [{s['id']}] {s['name']}: {s['status'][:60]}")
    out.append("\n-- CONSISTENCY CHECKS (out-of-dictionary) --")
    for c in CONSISTENCY_CHECKS:
        out.append(f"  {c['name']}: kills = {c['kills'][:60]}")
    out.append("\n-- CLOSURE RECORD --")
    for c in CLOSURE_RECORD:
        out.append(f"  {c['name']}: {c['how'][:64]}...")
    out.append("\n-- EPISTEMIC RECORD --")
    for e in EPISTEMIC_RECORD:
        out.append(f"  {e[:70]}")
    out.append("\n-- PAPERS --")
    for title, doi in PAPERS:
        out.append(f"  {title[:56]} ({doi})")
    out.append(f"  CITE: {CITE[:70]}")
    out.append("\n-- SCOREBOARD --")
    out.append("  " + SCOREBOARD_LOGIC[:200] + "...")
    return "\n".join(out)


if __name__ == "__main__":
    print(render())
