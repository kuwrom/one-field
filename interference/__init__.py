"""Harmonic formulation of the E_8(1) > G_2(1) x F_4(1) framework
through substrate self-interference.

Every prediction is expressed through 4 irreducible integers + pi + m_e
(the electron anchors the scale. M_Pl and G are OUTPUTS):

    d_10 = 2   quantum dimension of SU(3)_3 fundamental (1,0)
    d_11 = 3   quantum dimension of SU(3)_3 adjoint (1,1)
    n_7  = 7   dim of G_2 fundamental (Fano plane lines)
    n_26 = 26  dim of F_4 fundamental (traceless Albert algebra)

Reference:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
        Zenodo. https://doi.org/10.5281/zenodo.19899091
        (concept DOI, always resolves to the latest version.
         The papers cite version DOIs, e.g. v4 = 10.5281/zenodo.19932394)
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
        Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
        Zenodo. https://doi.org/10.5281/zenodo.20144381

Run the derivation via:
    python run.py

VERSION: this package is the canonical current state of the theory
and supersedes the numbers in the three published papers (which carry
the earlier anchor, boundary condition, and tree values).  The papers'
narratives are ported into the module docstrings. The package reads
standalone.

================================================================================
What this package is (and is not)
================================================================================

THE ONTOLOGY (quark paper, Sec. 2).  To be is, minimally, to be
self-present: anything that exists already affects itself by being
there, and interference is the physical image of that self-presence.
Mathematics can track self-reference through exactly three layers:
real = amplitude, complex = +phase, quaternion = phase of phase,
octonion = phase of phase of phase.  The octonions are the deepest
ledger that stays consistent (Hurwitz. Executable in octonions.py).
The ledger's depth and the loop's length are the same three: both
count NESTING, layers of phase-affecting-phase that can still be
tracked consistently, not moments of time.  Iterations of the
recursion and periods of the knot's internal rotation continue
indefinitely. It is a fourth nesting that has no consistent
ledger (Hurwitz).
G_2 = Aut(O) labels the patterns that close through all three
phases. The closures that forget their phase ([Delta, S] = 0, masses.py) collapse into
standing waves: three generations.  Mass is that pattern.  G_2 has
trivial centre and trivial pi_1, so every quantum number is emergent.
When G_2 breaks to SU(3) the Z_3 centre crystallises.

THE SEARCHLIGHT.  The octonions scoped the search: at most three
foldings, seven directions.  That is why the closing loop was
findable at all. The observation was targeted and guided.  The
octonions are the reason we look in 3-node loops.  The universe is the
standing wave the search found.

SEARCH ORDER.  The octonions' suggestion said only: structure needs
a loop, and the smallest loop is 3, so start there.  Had nothing
matched the anchors at 3, the search would have continued to 4, to 5.
The minimality axioms are the order any searcher must proceed in,
with the anchor match as the stopping rule. The match at 3 is the
contingent fact.

SELECTION BY TREE LENGTH.  Partial closures exist: lanes that make
electron-like things with ratios that fail to continue, quark-like
things that never bind.  They are eliminated for failing to continue.
The survivor is the lane with the longest emergence tree, the one
whose every closure enables the next, from the first loop to
chemistry.  The seven-branch gauntlet is the computed record of the
shorter trees.

This is NOT a top-down Lagrangian theory of the form "write an action,
solve it, get particles."  It is a SELECTION PRINCIPLE applied to a
substrate that hosts blind self-interference.  The question is not
"what is the Lagrangian?" but "given self-interference, which patterns
close?"  The answer is: exactly one lane.

OBSERVATION AND CORROBORATION.  The framework observes the pattern.
The established theorems it cites (Witten quantization, anomaly
matching, Callan-Harvey/Jackiw-Rossi zero modes, boundary-condition-
changing operator theory) are other observers' compressions of the
same regularities, reached from the Lagrangian side.  Citation is
corroboration: two observation logs of one pattern, one of them
written by people who had never seen this engine.

THE GATES ARE CONSISTENCY CONDITIONS WITHIN OUR LANE.  Given a
three-cycle closure, the level, the circulant, and the E_8 coherence
follow.  At loop 3 the surviving lane is unique: that is the computed
content of "exactly one lane" (the seven-branch gauntlet below).
Whether loop 3 itself was necessary or contingent is open (SEARCH
ORDER above: the match at 3 is the contingent fact), and the anchor
m_e is the fingerprint of which lane we are in.  The four integers
are readings of observed structure: counts.

THE HORIZON.  n_7 = 7 is the event horizon of self-knowledge: the
maximum number of directions a self-checking record can distinguish
(Hurwitz).  Whatever folds beyond it cannot hold an identity and is
visible in exactly one way: as amount without identity.  That is the
gravity sector (PvP = 0, Pv^2P = (1/2)P. The 182 bridge channels as
one unreadable ensemble), and the dark ratio reads accordingly:
f_baryonic = 1/(2 pi), the nameable fraction of being is one part in
a full turn.  The scorecard is the complete table of nameable things.
"no fourth generation" is absolute (a fourth would be a nameable
identity). Unnameable structure is expected, with gravity as its only
face.

ONE LOTTERY, THEN ARITHMETIC.  The selector question, which
configuration closes, is asked exactly once, at the first a->b->c
loop.  For every later sector the question changes character: no longer
which configuration closes, but what can form.  And that is
pre-solved: any configuration that produces the leptons
automatically produces these quarks and this gravity, because
they emerge from the SAME E_8 -> G_2 x F_4 structure.  The universe
doesn't have to get lucky twice.  This is why the mass sector is one
generator plus a grammar.

The dynamical equation IS written down.  It is Eq. (2.1) of the gravity
paper ("The Echo of Standing Waves"): the Z_3-symmetric coupled NLS

    i hbar d_t psi_k = -D nabla^2 psi_k + sum_j g_{kj} |psi_j|^2 psi_k,
                                                       k = 1, 2, 3,

with g_{kk} = g_0 and g_{kj} = g_1 = g_0 / sqrt(2) for k != j (the G_2
constraint).  The Madelung decomposition psi_k = sqrt(rho_k) exp(i S_k /
hbar) is exact, and the Z_3 Fourier split into q = 0 (geometry) and
q = 1, 2 (matter) is exact.  This resolves the Barcelo-Liberati-Visser
"two roles of rho" obstruction in single-component analog gravity.

From this equation:
  * the BdG linearisation gives Delta_k = A + B cos(theta + 2 pi k / 3),
    forced by [Delta, S] = 0.
  * the screened Poisson equation for common-mode memory falls out of
    the q = 0 projection.
  * the acoustic metric is the Madelung metric, providing the substrate
    for induced gravity.

================================================================================
The six gates that select the unique lane
================================================================================

Gate 1, Hurwitz (1898).
    Only four normed division algebras exist: R, C, H, O.  The octonion
    is the deepest consistent ledger. A fourth Cayley-Dickson doubling
    introduces zero divisors and breaks bookkeeping.  Its automorphism
    group is G_2 = Aut(O).  Theorem.
    Code: root.py (n_7 = dim G_2 fundamental, Fano plane lines).

Gate 2, Centre emergence.
    G_2 has trivial centre Z(G_2) = {1}.  SU(3) has centre Z_3.  When
    G_2 -> SU(3), the Z_3 crystallises as a child-only label.  Three
    sectors -> three generations.  Pure group theory.
    Code: masses.py (Z_3 circulant gap matrix -> three eigenvalues).

Gate 3, Level.
    Requiring the centre current J = (k, 0) to be a dimension-one boson
    fixes h(J) = 1, which uniquely forces k = 3.  Hence SU(3)_3 is the
    unique WZW level at this layer.  In unified notation k = d_11 and
    the altitude K = k + h_dual(SU(3)) = d_10 * d_11 = 6.
    Code: root.py (k_SU3 = d_11, K = d_10 * d_11).

Gate 4, Closure.
    The cyclic branch [Delta, S] = 0 is the one whose mass observable
    closes, fully determined by the layer's algebraic data with zero
    tunable inputs.  The direct finite-G_2 path retains the harmonic
    v = [1, -1/2, -1/2] and does NOT close.  Closure is the survivor
    filter.
    Code: masses.py (Brannen circulant), gravity.py (protected
    forgetting: P v P = 0, P v^2 P = (1/2) P).

Gate 5, Octavian / E_8 coherence.
    If the common-mode fibre is bosonic, local, self-dual, and c = 8,
    the even-unimodular lattice theorem (Griess) forces E_8.  Among
    rank-8 candidates {I_8, A_8, D_8, E_8, E_8 + E_8}, only E_8 passes
    all four conditions.  The central-charge identity c(G_2) + c(F_4)
    = 14/5 + 26/5 = 8 = c(E_8) at level 1 is the conformal-embedding
    signature.  The bridge count dim(7, 26) = n_7 * n_26 = 182 follows.
    Code: root.py (c_coset = 0), gravity.py (bridge identification).

Gate 6, Three-gate closure selection (gates L/Q/G =
    lepton residue, quark bridge, gravity heat-kernel normalisation.
    G_ind/G_N shown pre-face-split, UV/broken):
        localize G_2 first          -  Q  -   0.935/0.946
        PROTECTED FORGETTING ONCE   L  Q  G   0.994/1.007  <- closes
        full alpha from v^2         -  Q  -   1.062/1.077
        independent F_4 half        L  -  -   0.964/0.976
        dim-localize 7 first        -  -  -   0.976/0.989
        gauged E_8 vector           L  -  -   repulsive
        massive E_8 Proca           L  -  -   repulsive
    Only "protected G_2 forgetting applied exactly once" closes all
    three gates. The face-split law then takes the winning branch's
    0.994188 to 0.999999917.  (Seven-branch audit: gravity paper,
    Table V, reproduced above. The winning branch's venting ledger is
    executable as gravity.py's Sigma.)

================================================================================
The closure
================================================================================

    14  +  52  +  182  =  248  =  dim(E_8)
    gauge   matter   bridge
            + Higgs  -> gravity

Every degree of freedom of E_8 has a physical job.  Nothing is left
over.  At each layer, the observables are fully determined by that
layer's algebraic data: leptons on G_2, quarks and mixing on
G_2 x F_4 via the embedding, gravity on the mixed bridge.  No free
parameters at any layer, and no parameters inherited from above.
This is the closure, and the falsification gate.

The precision hierarchy is structural.  Leptons sit closest to the
Z_3 source and emerge cleanly (m_mu, m_tau predicted from m_e to
0.0012%).  Quarks pass through additional layers (triality,
generation scaling, confinement), and each layer adds residuals of
order alpha_s/pi (<= 0.4%).  Gravity sits two layers further out, at
0.0-1.3%.  The accuracy gradient is what the emergence depth
predicts.

================================================================================
All couplings derived: methodological consequence
================================================================================

The framework has no dimensionless tunable inputs.  The electron mass
sets the unit system and nothing else (the first closure anchors the
ruler. M_Pl, G, m_mu, m_tau, v_EW are outputs of the FORCED chain).
Every coupling, mass ratio, mixing angle, threshold, and Newton
normalisation is a polynomial or rational function of
(d_10, d_11, n_7, n_26) + pi.

This has a strict consequence for how the code reads.  In the Standard
Model, RGE running and matching are tools for tuning free interference
strengths (standard QFT: "Yukawa couplings"): a renormalised value at
one scale is clothed into a prediction at another, and the closure
(standard QFT: "dressing") absorbs the free parameters.  This framework
has no free interference strengths to tune, so the language of closure-
as-absorption does not apply: there is no bare prediction being patched.
What DOES happen is BACK-REACTION: each closure echoes through the
channels that already exist (the bridge self-loop, the e↔q cycle, the
WZW vents), and the echo's weight is a theorem of the channel, not a
fitted back-reaction (standard QFT: "correction").  Every factor that
appears here is one reading of the same algebra at a particular
emergence layer -- the layer determines the factor, not the other way
around.

The unified expression makes this manifest: instead of per-layer
modules each with its own physical motivation,
all named algebraic quantities reduce to closed-form polynomials in
(d_10, d_11, n_7, n_26).  Among others:

    K        = d_10 * d_11                          (WZW altitude = 6)
    h_10     = d_10 / d_11**2                       (h(1,0) = 2/9)
    h_11     = 1 / d_10                             (h(1,1) = 1/2)
    Q_0      = d_10 / d_11                          (Koide = 2/3)
    sin2W    = d_11 / (d_10**2 + d_11**2)           (Weinberg = 3/13)
    N_bridge = n_7 * n_26                           (= 182)
    dim_E_8  = 2*n_7 + 2*n_26 + n_7*n_26            (= 248)
    h_bridge = d_10/(1+d_10**2) + d_10*d_11/(1+d_11**2)   (= 1)
    c_coset  = 0                                    (rational identity)

The compression itself is a falsifiability handle: four integers plus
pi plus m_e producing the SM scorecard at sub-percent residuals is a
structurally tighter claim than "each layer closes consistently."

Concrete example, alpha_em = pi/512 and alpha(0) = 1/137.035999050 are
not "bare alpha and dressed alpha". They are two readings of one
coupling's back-reaction at two emergence layers.  alpha_em = pi/2**9 is the reading
at the conformal-embedding scale (Singh ratio + WZW), used inside the
instanton vertex.  alpha(0) is the reading at the IR pole-mass scale:
the bridge SELF-loop gives the factor (1 - 1/(2 pi)) (marginal h = 1,
unit coupling D2_local = 1, topological coset c_coset = 0), and the
depth-3 electron<->quark echo loop (two orientations, the S, S-dagger
rule) adds -2 (alpha / 2 pi)**2, so 1/alpha(0) is the real root of

    x**3 = (512/pi) [ (1 - 1/(2 pi)) x**2 - 1/(2 pi**2) ]
         = 137.035999050,

which lands on the Berkeley Cs recoil measurement (+0.13 sigma) in the
live >5-sigma Cs/Rb dispute.  The same (7,26) bridge, SAME sector,
SAME conformal weight h_bridge = 1, also gives Newton's constant via
182 scalar-like venting channels with xi_bridge = alpha_G_2 * (1/2) *
h_bridge = 1/(48 pi).  EM and induced gravity are two readings through
one bridge, not two sectors.

================================================================================
The echo law and THE KERNEL (root.py), nothing is a "correction"
================================================================================

One scalar recursion, x <- b + W(x), solved to its joint fixed point
(Web.solve(), classical fixed-point iteration).  Every quantity is a
LEDGER: base amplitude + echo stack with provenance, because a closed
knot must vent to keep circulating and later closures echo back
through everything that already exists.  Every weight is a theorem.
Depth is capped at 3 by the Hurwitz gate (a cap on the map, not the
territory).  The multiplicity rules (orientation, sector-node, vertex
composition, face-split) and their uniqueness enumerations live with
the law itself in root.py. FORCED terms enter canonical values,
CANDIDATE terms are excluded.

================================================================================
File-to-role map
================================================================================

File         | Role
-------------+-----------------------------------------------------------
root.py      | Four irreducible numbers. Every derived quantity as a
             | polynomial / rational function of them.
masses.py    | 9 mass predictions (3 leptons + 6 quarks) from the Z_3
             | circulant + WZW emergence ladder.
mixing.py    | CKM (4 Wolfenstein) + PMNS (3 angles + delta_CP) from
             | the four faces of SU(3)_3 / D^(6).
couplings.py | alpha_s(M_Z) from the embedding-index chain. alpha(0)
             | from (7,26) bridge self-interference (one bridge, two
             | readings).
gravity.py   | G_ind/G_N from the universal venting ledger (every massive
             | closure feeds the common mode). Higgs from F_4(1) fusion.
             | eta_B from G_2 instanton. Neutrinos from F_4 singlets +
             | rank-2 seesaw.
dark_sector.py | Omega_DM/Omega_b = 2 pi - 1. Cosmological-constant scale.
octonions.py | Constructive Fano-plane proof of B/A = sqrt(2). Linearity
             | from closure, executable: L_x linearity, |xy| = |x||y|,
             | sedenion zero divisor (the Hurwitz gate), Z_3 selection
             | rule, rank-8 E_8 lattice gate.
words.py     | Generation word lemma, executable and fully derived:
             | base(n) = (4, 12, 97) as boundary-walk counts. N-ality
             | superselection (neutral lane). knot definition +
             | conversion lemma (terminus). Exclusions killed by data.
run.py       | End-to-end driver, scorecard, agreement summary.
nls_soliton.py | The substrate equation (Z_3 coupled NLS), split-step
             | evolution, soliton stability.
protected_forgetting.py | PvP = 0, Pv^2P = (1/2)P: the forgetting that
             | separates matter readings from the common mode.
embedding_uniqueness.py | The E_8 > G_2 x F_4 embedding singled out
             | among the rank-8 candidates.
tests/probes/ | Analyses OF the theory, run on its own equations:
             | zero_mode.py (V-A computed: handedness = winding sign),
             | bookkeeping.py (WZ level = 3 counted. theta = 2/9 is the
             | sector-changer), mtc_spectrum.py (three-line spectrum at
             | the MTC layer), theta_menu.py (the data permits only 2/9),
             | selector_theorem.py ({Haar, delta_e} dichotomy),
             | orientation_bit.py (one handedness datum),
             | winding_texture.py (3D pi_3 texture), knot.py (formation,
             | control-proven G_2 binding), stationary.py (persistence),
             | spectrum.py (spectroscopy), universality.py (standing
             | character), dynamics.py (bias erosion + shadow),
             | bell.py (local non-realist: the rejected assumption
             | is the +/-1 value assignment, not locality. Exact
             | Tsirelson), ew_internal.py (EW imports reclassified),
             | wiring_scan.py (the wiring: a closure-residual scan
             | over the classified candidate menus of lattices,
             | conformal pairs, and levels. Four integers never input),
             | skyrmion_3d.py (3D knot formation: texture seed grows
             | only under the G2 coupling, growth rate matches the
             | BdG instability, and the common mode dips, the shadow),
             | knot_charge.py (pi_3 winding NOT protected at substrate
             | level: measured unwinding, rotation loophole closed),
             | confined_phase.py (Z3 disorder operator: the classical
             | endpoint is defect-free. Entry test of the substrate-
             | conjecture ledger), knot_spectrum.py (ring-down tower
             | of the saturated lump, gapped), surplus_edges.py (the
             | grammar enumerated away from the data: predictions and
             | falsifiers).  Shared 3D machinery: substrate3d.py.
tests/       | The falsifiability suite: every claim is an assertion.
             | test_interference.py freezes the NUMBERS (scorecard,
             | freeze tables, attractor). test_probes.py freezes the
             | MECHANISMS (each probe's registered conclusion. The
             | substrate-dynamics arc behind pytest -m slow).

Register (repo root): interference/registry.py, the single ledger: live
watches, registered predictions with kill conditions, and the closure
record of the four structural questions.

Quickstart:   python3 run.py        (full chain + scorecard)
              python3 -m pytest -q  (the test suite)
Dependencies: Python 3.10+, numpy. pytest for the suite.
Ledger inspection: any node's interference history prints via
root.WEB[name].table() after the chain has run.

================================================================================
The attention/SSM reading (interpretive only. Cited)
================================================================================

NO MACHINE LEARNING ENTERS THE COMPUTATION.  The kernel x <- b + W(x)
is solved by a classical Jacobi-style fixed-point iteration. Nothing
is trained, no parameters are fitted, and the weights are theorems.

The READING, kept as documentation because it guided the kernel's
design (a discovery heuristic. No number depends on it): the echo web can be read as a 3-layer attention network in the
sense of Vaswani et al., "Attention Is All You Need" (NeurIPS 2017):
closures are tokens, channels are heads, arrival order is sequence
position, ledgers are residual streams.  The mapping is generic
(any sparse weighted fixed-point system admits it), which is exactly
why it carries no explanatory weight here.  The old per-layer
pipelines were CAUSAL attention. The echo law is the same network
with the mask removed (back-reaction = mask removal).  In the
language of structured state-space sequence models (S4: Gu, Goel,
Re 2021. Selective SSMs/Mamba: Gu, Dao 2023) the interaction kernel
is exactly 3-nilpotent (the SSM restatement of the Hurwitz depth
cap).  Octonions are the proof of the
architecture's shape, not its substance.

================================================================================
The bets, the coordinate, and the closed ledger
================================================================================

Sole dimensional input: the electron mass (CODATA, +-0.3 ppb).  M_Pl
and G are OUTPUTS.  No depth-4 term exists to retune anything: the
no-further-edge theorem (root.py) closes the EM channel outright.

Every registered number lives in ONE place.  The bets, kill
conditions, supersession history, and the closure record of the four
structural questions: interference/registry.py.  The frozen values the
code must reproduce: tests/test_interference.py.  This file states
the claim taxonomy and the forward rule. It quotes no number the
register owns.

THE MASS COORDINATE (back-reaction, not running).  A running mass
m_q(mu) is a coordinate on an RG orbit, not an observable. A mass
prediction is complete only with its coordinate.  Leptons and top at
the propagator pole. Confined heavy quarks (c, b) at the self-scale
m(m). Light quarks via the RG-invariant ratios, which carry no
coordinate at all.  The algebra does not run to a scale. Its output
IS the prediction.  Full treatment: masses.py.

THE LEDGER IS CLOSED.  Every canonical claim is exactly one of:
DERIVED (the algebraic chain), ENUMERATED (words lemma, K^3
selection, metric uniqueness), CROSS-PINNED (shared dictionary
objects measured in independent observables), or a CITED IMPORT WITH
A STATED BAND, where every import is decomposed and audited
(gravity.py, tests/probes/ew_internal.py): data-driven ingredients
measurements, parameter dependencies are recalculated from the
framework's own values and verified compatible, and pure loop
integrals are mathematics fixed by the field content.  Nothing is
dressed. Every non-trivial factor is a back-reaction through an
existing channel. Nothing is tuned.

HURWITZ FINALITY (the forward rule).  The echo grammar is closed:
depth <= 3, edges only along the sector graph, multiplicities only by
the orientation, sector-node, vertex-composition and face-split
rules.  No new term may be added after a comparison.  Consequently
every canonical residual is BOUNDED: a deviation that exceeds the
next-depth echo unit of its sector falsifies the rule that produced
it. There is no deeper term to absorb a miss.
""" 
