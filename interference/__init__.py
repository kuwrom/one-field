"""Harmonic formulation of the E_8(1) > G_2(1) x F_4(1) framework
through substrate self-interference.

THE CLAIM: hand an interference solver a substrate
(three components in cyclic Z_3 exchange, self-reference trackable) and exactly one family of standing waves closes; its
fixed point is the Standard Model, gravity, and the dark sector.

Four layers, named consistently across the repository:
    ENGINE      the solver x = b + W(x) plus the solved web; canon,
                every scorecard number attaches here.  The derivation
                is flat, so interference/ itself carries it
    SEARCHLIGHT the octonionic algebra that located the coordinates
                (searchlight/, the one subdirectory, set aside
                because it is logically dispensable: the engine never
                imports it); it bounds what is trackable, not what
                exists
    SUBSTRATE   the medium the claim is about; its physical
                realisation is a verification program, so it sits
                with the scrutiny (tests/substrate/), probed by
                tests/probes/
    AUDIT       the instruments that hold the engine to its own rules
                (tests/audit/): the committed edge menus, the forward
                search, the blindness freeze
Checking whether the theory is flawless is testing: the freezes and
the probes live in tests/.  The layer names are not the directory
names, and deliberately so: what a thing IS and where it is READ from
are different questions.

Every prediction is expressed through 4 integers (depth-3 searchlight
vocabulary) + pi + m_e (the electron anchors the scale. M_Pl and G
are OUTPUTS):

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
    python3 interference/run.py

This package is the canonical statement of the theory: where a number
here and a number in the published papers differ, this one is the
canon.  The papers' narratives are carried in the module docstrings,
and the package reads standalone.

================================================================================
What this package is (and is not)
================================================================================

THE ONTOLOGY (quark paper, Sec. 2).  To be is, minimally, to be
self-present: anything that exists already affects itself by being
there, and interference is the physical image of that self-presence.
Closure takes a record to be a finite-dimensional unital real composition
algebra: composition is bilinear and its positive-definite norm, the
ledger's probability weight, is multiplicative. Hurwitz then classifies
the admissible algebras as R, C, H, and O. Taking R as the amplitude
baseline, the three non-real extensions C, H, O are the three searchlight
depths; the next Cayley-Dickson algebra fails by zero divisors. This is a
complete theorem from the framework's ledger criterion, not
a claim that Hurwitz alone speaks about self-reference.

Records, rather than Hilbert rays, are primitive.  A physical phase must
therefore be carried by a relation between closed records; permutations,
braids, and interference paths provide such relations, while a scalar
acting identically on the complete vacuum does not add a physical line.
For gauge winding, the complete matter-coupled character is bar(theta),
and its action factors through the closed vacuum-record algebra.  This is
part of the ontology, not a consequence of local MTC data.  Applied to
A=1+B in couplings.py, it leaves only bar(theta)=0.

The recursive Web, implemented and ancestry-checked in root.py, may have
greater depth. The octonionic searchlight reveals the a->b->c->a loop at
depth three; the Web refines that pattern at depth four and beyond.
G_2 = Aut(O) labels what the
searchlight can track; it does not limit what reality can do. The
protected-forgetting branch ([Delta, S] = 0, masses.py)
is a three-eigenvalue circulant identified with the three generations,
and mass is the squared gap of that pattern. G_2 has
trivial centre and trivial pi_1, so every quantum number is emergent.
When G_2 breaks to SU(3) the Z_3 centre crystallises.

THE SEARCHLIGHT. The octonions reduced blind interference fluctuations
to seven directions and, through the unique (7,26) bridge, a finite
7x26 search. That is why the closing loop was findable at all. We got
lucky in a precise sense: the smallest structure becomes observable at
the same depth three that the searchlight can resolve. The Chladni-plate
analogy is literal at the level of method: the form first becomes
visible at finite resolution and is refined by later recursive readings.
The engine's arithmetic never evaluates an octonion or a sedenion.
The algebra's entire output is a handful of discrete selections,
each with a closure status in the engine:

  selection              | engine constant | closure
  -----------------------+-----------------+--------------------------
  normed-ledger gate     | searchlight = 3 | DERIVED given closure's
  + Hurwitz              |                 | unital positive composition-
                         |                 | algebra criterion; not Web cap
  quantum dims of SU(3)_3| d_10=2, d_11=3  | the four integers
  Fano CG |C_3bar|/|C_1| | B/A = sqrt(2)   | COMPUTED + IH: Fano CG,
                         |                 | positive spherical trace and
                         |                 | orientation count enter the
                         |                 | channel-to-circulant map
  Sugawara weight h_10   | theta = 2/9     | COMPUTED + IH: lightest
                         |                 | triality changer and theta=h
                         |                 | mass-phase map; sole positive
                         |                 | angle-menu survivor
  Q_0 = 1/3 + d_10/6     | Q_0 = 2/3       | DERIVED downstream from B/A
                         |                 | by the circulant identity

The last row is the model for the rest: the program
(executable as test_coverage.py's test_searchlight_absorption) is
that every searchlight selection is eventually read back out of the
engine the way Q_0 already is, at which point the algebra has fully
retired into motivation.  The strong reading, tested by that
program, is that the octonions will turn out to be not a separate
searchlight shone on the territory but the territory's own
accounting seen from the bookkeeping side. Until each selection is
read back out, they remain the searchlight. Deeper self-interference is represented by the
recursive web, not by extending the Cayley-Dickson tower. The Chladni
image is exact: the pattern becomes visible at depth three and gains
finer structure as later backreactions are evaluated.

DEPTH VOCABULARY. ``searchlight depth`` is the octonionic tracking
resolution (at most three). ``web depth`` is the recursive-interference
depth stored by Ledger and EchoTerm. ``beta order`` is the imported
perturbative expansion used for transport and is not web depth. A module's
``layer`` is only its place in the dependency or physical-regime narrative.
No value acquires a second physical identity by changing calculation or
comparison coordinates.

SEARCH ORDER. The framework chose the three-label cyclic closure as its
first search target. This is not a graph-theoretic claim that every loop
needs three nodes: the Web itself contains self-loops and a two-cycle.
Had nothing matched the anchors at 3, the search would have continued
to 4, then 5. The match at 3 is a contingent result of that search rule.

THE GENERATIVE DIRECTION.  The engine derives the scorecard.  The
substrate probes test whether the same structure emerges dynamically:
knot formation, circulation, common-mode sourcing, confinement and
clock locking.  Each probe states its own verdict and kill
condition.  The circulant ratios are closed algebraically.
emergent_clock.py reports broadband clocks in substrate-only 1D and 3D
runs.  A confined-phase 3D run that locks to the algebraic circulant
remains an outstanding dynamical verification.

WHERE THE QUESTIONS END.  The engine takes one number, and even
that number is a unit name rather than physics: m_e/M_Pl is derived
outright and reproduces nature at two sigma of G's 22 ppm
uncertainty (asserted at import, root.py), so the dimensionless
output takes no input.  How strongly the substrate interferes is not
a second input; it is the thing being derived, alpha_G_2(M_Pl),
alpha(0), alpha_s, one channel at a time.  What no theory can answer
is why the amplitude is nonzero at all.  Self-presence is
UNIVERSAL, not conditional: the empty world is present to itself too,
and precisely
because it is nothing the self-touching generates no amplitude.  It
is the engine's zero orbit of x = b + W(x), zero in and zero out,
not an exception to the axiom but its trivial fixed point.  Nothing
does not escape the axiom; it just does not close down the road.
And because self-interference cannot be reduced further, anything
beneath it is undecidable from inside, not merely unanswered.

THE GAUGE-FIELDS VERIFICATION PROGRAM.  The gauge sector of the web
is carried by collective excitations of the confined phase, not
postulated fields.  The algebraic structure (representations,
couplings, confinement) is derived from the E_8/MTC data.  The
computational demonstration program includes reproducing the
perturbative gauge series on composite knots (g-2, scattering
amplitudes, decay widths) and the multi-knot Fock/Pauli structure
beyond two-body exchange.
ENDGAME (the grounded-strength reading): closure makes the emergent
coupling a function of substrate parameters, and the derived
alpha_G2(M_Pl) = 1/(24 pi) then INVERTS into a derived depth
rho xi^3: the medium's interference weight becomes inherent to the
minimal form, and the contingency surface collapses to existence
itself plus unit names.

SELECTION BY DERIVATION-PATH LENGTH. Partial closures exist: lanes that make
electron-like things with ratios that fail to continue, quark-like
things that never bind.  They are eliminated for failing to continue.
The survivor is the lane with the longest derivation path, the one
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
matching, Skyrmion-Dirac spectral flow, boundary-condition-
changing operator theory) are other observers' compressions of the
same regularities, reached from the Lagrangian side.  Citation is
corroboration: two observation logs of one pattern, one of them
written by people who had never seen this engine.

THE GATES ARE CONSISTENCY CONDITIONS WITHIN OUR LANE.  Given a
three-cycle closure, the level, the circulant, and the E_8 coherence
follow.  At loop 3 the surviving lane is unique: that is the computed
content of "exactly one lane" (the seven-branch gauntlet below).
Whether loop 3 itself was necessary or contingent is answered by the
SEARCH ORDER above: the match at 3 is the contingent fact.  The anchor
m_e is the fingerprint of which lane we are in.  The four integers
are readings of observed structure: counts.

THE HORIZON. Given the norm-multiplicative ledger gate, n_7 = 7 is the
maximum number of imaginary directions this octonionic searchlight
distinguishes. Whatever folds beyond it is read as amount without identity. That is the
gravity sector (PvP = 0, E[v^2] = 1/2 per mode. The 182 bridge channels as
one unreadable ensemble), and the dark ratio reads accordingly:
f_baryonic = 1/(2 pi), the nameable fraction of being is one part in
a full turn.  The scorecard is the complete table of nameable things.
"no fourth generation" is absolute (a fourth would be a nameable
identity). Unnameable structure is expected, with gravity as its only
face.  This is also why the bridge reads as a scalar: tracking its
182 channels as one coherent multi-layer object would need a normed
division algebra that DOES NOT EXIST, so the ensemble is read as
real scalars, the vanishing coset charge 8 - 14/5 - 26/5 = 0 makes
that decomposition exact, and each channel's heat-kernel coefficient
+1/6 - xi is positive and makes gravity ATTRACTIVE.  Reading the same
channels as vectors would give -1/6, reverse the sign, and matter
would never bind.

ONE LOTTERY, THEN ARITHMETIC.  The closure question, which
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

with g_{kk} = g_0 and g_{kj} = g_1 = g_0 / sqrt(2) for k != j.
The Fano-plane cross product gives
|C_3bar|/|C_1| = sqrt(2) (searchlight/octonions.py).  Gate 3 selects
the SU(3)_3 modular data.  Their canonical positive spherical trace gives
the channel weight 1/d(fund) = 1/2.  The channel-to-circulant
identification combines that trace with two cycle orientations and reads
B/A = 2 sqrt(2)/2 = sqrt(2).  A separate identification transplants this
amplitude ratio into the substrate coupling g_1/g_0 = 1/(B/A) in
tests/substrate/nls_soliton.py and tests/probes/circulant_normalization.py.
Those probes test consequences of the transplant; they do not yet
reconstruct the algebraic circulant from confined dynamics, so the
algebra-to-dynamics lock remains open.  The Madelung decomposition psi_k = sqrt(rho_k) exp(i S_k /
hbar) is exact, and the Z_3 Fourier split into q = 0 (geometry) and
q = 1, 2 (matter) is exact.  This resolves the Barcelo-Liberati-Visser
"two roles of rho" obstruction in single-component analog gravity.
Standing of the equation itself: the NLS form is
DERIVED up to locality as the universal envelope normal form of weak
local self-interference (tests/probes/envelope_universality.py); the
locality is an axiom of self-presence: a nonlocal self-term would be
interference with a presence elsewhere, i.e. with a second thing,
which is coupling, not self-interference.

From this equation:
  * the BdG linearisation gives Delta_k = A + B cos(theta + 2 pi k / 3),
    forced by [Delta, S] = 0.
  * the screened Poisson equation for common-mode memory falls out of
    the q = 0 projection.
  * the q = 0 phonon dispersion and sound speed are computed.  The
    framework identifies the Madelung common mode with acoustic geometry,
    which supplies the geometry channel for the induced-gravity ledger.

================================================================================
The six gates that select the unique lane
================================================================================

Gate 1, the normed-ledger criterion plus Hurwitz classification.
    Closure requires a finite-dimensional unital real composition algebra:
    bilinear composition and a positive-definite multiplicative norm.
    Hurwitz therefore leaves R, C, H, O, while the next Cayley-Dickson
    algebra has zero divisors. Relative to the real amplitude baseline,
    C, H, O are exactly three added record depths. G_2 = Aut(O).
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
    forgetting: P v P = 0; per-mode E[v^2] = 1/2, total
    P v^2 P = P; the two contractions are distinct and easily
    conflated).

Gate 5, Octavian / E_8 coherence.
    The common-mode fibre is bosonic (integer-spin venting), local
    (self-interference), and self-dual (the substrate is its own
    observer).  The central charge c = 8 is not an additional
    assumption: it is DERIVED from the conformal embedding selected by
    the preceding gates, c(G_2) + c(F_4) = 14/5 + 26/5 = 8.  Given
    these, the even-unimodular lattice theorem (Griess) forces E_8.
    Among rank-8 candidates {I_8, A_8, D_8, E_8, E_8 + E_8}, only E_8
    passes all conditions (wiring.py stage W1 computes the residuals).
    The bridge count dim(7, 26) = n_7 * n_26 = 182 follows.
    Code: wiring.py (lattice scan), root.py (c_coset = 0),
    gravity.py (bridge identification).

Gate 6, Three-gate closure selection (gates L/Q/G =
    lepton residue, quark bridge, gravity heat-kernel normalisation.
    G_ind/G_N shown before apportionment, UV/broken):
        localize G_2 first          -  Q  -   0.935/0.946
        PROTECTED FORGETTING ONCE   L  Q  G   0.994/1.007  <- closes
        full alpha from v^2         -  Q  -   1.062/1.077
        independent F_4 half        L  -  -   0.964/0.976
        dim-localize 7 first        -  -  -   0.976/0.989
        gauged E_8 vector           L  -  -   repulsive
        massive E_8 Proca           L  -  -   repulsive
    Only "protected G_2 forgetting applied exactly once" closes all
    three gates. The apportionment law then takes the winning branch's
    0.994188 to 0.999999917.  (Seven-branch audit: gravity paper,
    Table V, reproduced above. The winning branch's venting ledger is
    executable as gravity.py's Sigma.)

WHAT THE GATES LEAVE.  One structural lane closes through the gates;
the closure search holds it fixed and runs the committed mass grammar
as a generator. It prints the look-elsewhere denominator: canon is
one of 504 mass grammar hits at the scorecard's 1% band out of
2,370,816 eligible realizations, the data resolves the word counts
uniquely but not the sub-percent corrections, and the joint 0.1% tier is empty, canon
included, because the u and b references are coarser than that.
A candidate that does not close describes nothing that exists:
what exists is what carries non-zero interference weight.  Code: tests/audit/closure_search.py, which prints the full
accounting and states this in its docstring.

================================================================================
The closure
================================================================================

    14  +  52  +  182  =  248  =  dim(E_8)
    gauge   matter   bridge
            + Higgs  -> gravity

Every degree of freedom of E_8 has a physical job.  Nothing is left
over.  The implemented scorecard is numerically fixed by the layer data;
what is identified rather than derived is the mixing sector.

The precision hierarchy follows the dependency graph. Leptons sit
closest to the Z_3 source. Quarks cross additional triality,
generation-scaling, and confinement edges. Gravity first appears at
web depth four. These graph distances are not perturbative beta orders.

================================================================================
Coupling chain: methodological consequence
================================================================================

The framework has no dimensionless tunable inputs.  The electron mass
sets the unit system and nothing else (the first closure anchors the
ruler. M_Pl, G, m_mu, m_tau, v_EW are outputs of the FORCED chain).
Every coupling, mass ratio, mixing angle, threshold, and Newton
normalisation is a polynomial or rational function of
(d_10, d_11, n_7, n_26) + pi.

This has a strict consequence for how the code reads. The framework has
no free interference strengths to tune and no sequence of preliminary
and corrected physical masses. Each closure echoes through channels that already
exist (the bridge self-loop, the e↔q cycle, the WZW vents), and the echo's
weight is a theorem of the channel. A beta-function or matching formula
may translate a completed reading into the coordinate required by an
external calculation; it does not create another framework state.

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
two coordinates of one coupling. alpha_em = pi/2**9 is the reading at
the conformal-embedding scale (Singh ratio + WZW), used inside the
instanton vertex. alpha(0) is the zero-momentum reading:
the bridge SELF-loop gives the factor (1 - 1/(2 pi)) (h = 1,
unit coupling D2_local = 1, c_coset = 0), and the
depth-4 recursive reading of the electron<->quark carrier (two
orientations, the S, S-dagger rule) adds -2 (alpha / 2 pi)**2, so
1/alpha(0) is the real root of

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
through everything that already exists.  Every weight has an explicit
status: derived weights carry a witness, while the links
identified rather than derived are stated beside the physics that
produces them, with their exact consumers.
Searchlight resolution is 3 by the norm-multiplicative ledger criterion
composed with Hurwitz classification. Recursive-web depth is not capped.
Every later edge must name and read its
depth-(n-1) parent and be derived and registered before a
comparison. The multiplicity rules (orientation, sector-node, vertex
composition, apportionment) and their uniqueness enumerations live with
the law itself in root.py. FORCED terms enter canonical values,
CANDIDATE terms are excluded.

================================================================================
File-to-role map
================================================================================

The directory groups are the four layers.  Each derivation module
exposes one derive(...) -> dict, and run.py calls them in dependency
order.

File                    | Role
------------------------+--------------------------------------------
run.py                  | End-to-end driver, scorecard, agreement
                        | summary.
root.py                 | Depth-3 searchlight vocabulary.  Every derived
                        | quantity as a polynomial / rational function
                        | of them.  The echo web and its depth gate.
masses.py               | 9 mass predictions (3 leptons + 6 quarks)
                        | from the Z_3 circulant + WZW emergence
                        | ladder.
mixing.py               | WZW/D^(6)* connection and mixing candidates;
                        | the sector-operator lift is depth-4
                        | (mixing_z3_dichotomy.py, mixing_operator_search.py).
couplings.py            | alpha_s(M_Z) from the embedding-index chain.
                        | alpha(0) from (7,26) bridge self-interference
                        | (one bridge, two readings).
gravity.py              | G_ind/G_N from the universal venting ledger
                        | (every massive closure feeds the common
                        | mode).  eta_B from G_2 instanton.  Neutrinos
                        | from the F_4-singlet RH-neutrino identification:
                        | conditional minimal-type-I rank at most 2;
                        | ordering not selected.
higgs_ew.py             | SM RGE transport, Higgs quartic boundary
                        | (F_4(1) fusion + bridge), M_W, M_Z, import
                        | audit.  Engine import inventory.
dark_sector.py          | Omega_DM/Omega_b = 2 pi - 1. Friedmann
                        | critical density versus current dark-energy
                        | density; CKN status explicit.
words.py                | Generation word lemma, executable and fully
                        | derived: base(n) = (4, 12, 97) as boundary-
                        | walk counts.  N-ality superselection (neutral
                        | lane).  Knot definition + conversion lemma
                        | (terminus).  Exclusions killed by data.
wiring.py               | The wiring: a closure-residual scan over the
                        | classified candidate menus of lattices,
                        | conformal pairs, and levels.  The four
                        | integers are never an input.
reference.py            | Comparison data, single home, quarantined by
                        | tests/test_blindness.py.
searchlight/octonions.py| Constructive Fano-plane proof of
                        | B/A = sqrt(2).  Linearity from closure,
                        | executable: L_x linearity, |xy| = |x||y|,
                        | sedenion zero divisor (the Hurwitz gate),
                        | Z_3 selection rule, rank-8 E_8 lattice gate.
searchlight/embedding_uniqueness.py
                        | The E_8 > G_2 x F_4 embedding singled out
                        | among the rank-8 candidates, with its
                        | near miss frozen by assertion.
searchlight/protected_forgetting.py
                        | PvP = 0, E[v^2] = 1/2 per mode: the
                        | forgetting that separates matter readings
                        | from the common mode.  Weyl-orbit uniqueness.

The SUBSTRATE and the AUDIT instruments read from the scrutiny
side, so they live with the tests:

tests/substrate/nls_soliton.py
                        | The substrate equation (Z_3 coupled NLS),
                        | split-step evolution, BdG, soliton stability,
                        | Madelung sourcing.
tests/substrate/substrate3d.py
                        | The same equation in 3D, shared by the
                        | dynamics probes.
tests/audit/edge_grammar.py
                        | The committed grammar: menus, corrections,
                        | the null model's alphabet.
tests/audit/closure_search.py
                        | The grammar as a generator: blind gates
                        | first, data contact once, look-elsewhere
                        | denominator printed.

Checking whether the theory is flawless is testing, and the suite is
at the repo root in tests/.  test_interference.py freezes the NUMBERS
(scorecard, freeze tables, attractor).  test_probes.py freezes the
MECHANISMS, one conclusion per probe, with the substrate-
dynamics arc behind pytest -m slow.  test_coverage.py freezes the
EPISTEMICS (menu sparsity, cross-bracing, citation resolution, the
searchlight absorption scoreboard).  test_blindness.py is the
reference-data firewall.  The analyses themselves live in
tests/probes/, each carrying its own docstring and its own
verdict; that inventory is not restated here, because a list of
probes is the first thing to rot.

Quickstart:   python3 interference/run.py   (full chain + scorecard)
              python3 -m pytest -q          (the test suite)
Dependencies: Python 3.10+, numpy, scipy. pytest for the suite.
Ledger inspection: any node's interference history prints via
root.WEB[name].table() after the chain has run.

================================================================================
The attention/SSM reading (interpretive only. Cited)
================================================================================

NO MACHINE LEARNING ENTERS THE COMPUTATION.  The kernel x <- b + W(x)
is solved by a classical Jacobi-style fixed-point iteration. Nothing
is trained, no parameters are fitted, and the weights are theorems.

The READING, kept as documentation because it guided the kernel's
design (a discovery heuristic. No number depends on it): the echo web
can be read as a sparse attention network in the sense of Vaswani et
al., "Attention Is All You Need" (NeurIPS 2017):
closures are tokens, channels are heads, arrival order is sequence
position, ledgers are residual streams.  The mapping is generic
(any sparse weighted fixed-point system admits it), which is exactly
why it carries no explanatory weight here.  A per-layer pipeline is
CAUSAL attention; the echo law is the same network with the mask
removed (backreaction = mask removal).  In the language of
structured state-space sequence models (S4: Gu, Goel,
Re 2021. Selective SSMs/Mamba: Gu, Dao 2023), the finite committed
edge graph is the state transition. Octonions located its smallest
observable cycle; they do not cap its recursive depth.

================================================================================
The bets, the coordinate, and the closed chain
================================================================================

Sole dimensional input: the electron mass (CODATA, +-0.3 ppb). M_Pl
and G are OUTPUTS. The depth-4 response refines the ruler.
The EM sector graph admits no additional independent edge after it:
repeated traversals of that carrier are already resummed by the fixed
point. This is channel finality, not a global depth cap.

Every number lives in ONE place: the module that computes it, with
its kill condition stated beside it.  The frozen values the code must
reproduce: tests/test_interference.py.  mixing_z3_dichotomy.py proves that
closed frame-blind data contain no relative sector diagonalizer.  The finite
mixing maps have identification status.  sakharov_normalisation.py computes
the target 2*pi*d11 from two equal factorizations under the
Planck/healing-length identification.  The substrate suite contains clock
runs, confined-phase lattice probes, common-mode response calculations, and
the Skyrmion-Dirac spectral-flow probe.

THE MASS READING. The algebra produces one completed Web reading per
fermion. Standard experimental tables attach conventional comparison
coordinates to their entries; those labels are dictionary metadata, not
additional stages of the derivation. Full treatment: masses.py.

THE CHAIN IS CLOSED.  Every canonical claim is exactly one of:
DERIVED (the algebraic chain), ENUMERATED (words lemma, K^3
selection, metric uniqueness), TRIANGULATED (shared dictionary
objects measured in independent observables), or a CITED IMPORT WITH
A STATED BAND, where every import is decomposed and audited
(gravity.py, tests/probes/ew_internal.py): data-driven ingredients
measurements, parameter dependencies are recalculated from the
framework's own values and verified compatible, and pure loop
integrals are mathematics fixed by the field content.  Nothing is
split into preliminary and subsequently corrected physical values. Every
non-trivial factor is a backreaction through an existing channel. Nothing
is tuned.

FORWARD RULE. The searchlight closes at depth 3. The recursive web is
open to later resolution but closed against hindsight: edges only along
the sector graph, with multiplicities only by the orientation,
sector-node, vertex-composition and apportionment rules, and every new depth
registered with its derivation before comparison. Consequently
every canonical residual is BOUNDED: a deviation that exceeds the
next admissible echo unit of its sector falsifies the rule that produced
it. A deeper term off the committed menu cannot be invented to absorb
a miss.

HISTORICAL STATUS.  The Koide relation (1982) and Brannen
parameterization (2006) predate this implementation.  The lepton mass
scorecard is postdictive.

COMPUTATION ORDER.  The six gates select SU(3)_3.  The lightest
triality-changing primary has h = 2/9, and the identified mass-phase map
sets theta = h.  The identified channel-to-circulant map applies the Fano
CG, canonical positive spherical trace and orientation count to obtain B/A.
Positivity and nondegeneracy leave theta = 2/9 as the sole value in the
committed structural menu.  A separate identified amplitude-to-coupling
transplant remains dynamically unclosed.  masses.py then compares the
calculated mass ratios with the PDG entries.

WHAT WAS NOT DONE. No formula was constructed for an unmeasured
constant by searching combinations of (2, 3, 7, 26, pi). In a
four-integer grammar such searches always succeed, which is precisely
why their successes would carry no weight. Every exposure in this
repository uses only edges already in the graph, and the framework
survives only if ALL of them hold: one clean failure ends it, by its
own rules. The look-elsewhere denominator is printed, not asserted
(tests/audit/closure_search.py).
"""

# The derivation is a regular Python package. Sibling modules use
# package-relative imports so stateful modules such as ``root`` load only
# as ``interference.root`` and the solved WEB singleton remains unique.
