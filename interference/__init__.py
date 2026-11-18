"""Harmonic formulation of the E_8(1) > G_2(1) x F_4(1) framework
through substrate self-interference.

THE CLAIM: give an interference solver a substrate (three components
in cyclic Z_3 exchange, with trackable self-reference). Exactly one
family of standing waves closes. Its fixed point is the Standard Model,
gravity, and the dark sector.

Four layers, named consistently across the repository:
    ENGINE      the solver x = b + W(x) plus the solved web. This is
                the canonical layer: every scorecard number comes from
                it.  The derivation is flat, so it lives directly in
                interference/
    SEARCHLIGHT the octonionic algebra that located the coordinates
                (searchlight/, the one subdirectory). It is kept apart
                because it is logically dispensable: the engine never
                imports it. It limits what can be tracked, not what
                exists
    SUBSTRATE   the medium the claim is about. Its physical
                realisation is the realization program, so it sits
                with the scrutiny (tests/substrate/), probed by
                tests/probes/
    AUDIT       the instruments that hold the engine to its own rules
                (tests/audit/): the committed edge menus, the forward
                search, the blindness freeze
Checking the theory is testing, so the freezes and the probes live in
tests/.  The layer names are deliberately not the directory names.
What a thing IS and where it is READ from are different questions.

Every prediction is expressed through 4 integers (the depth 3
searchlight vocabulary), pi, and m_e. The electron sets the scale.
M_Pl and G are OUTPUTS.  The integers are:

    d_10 = 2   quantum dimension of SU(3)_3 fundamental (1,0)
    d_11 = 3   quantum dimension of SU(3)_3 adjoint (1,1)
    n_7  = 7   dim of G_2 fundamental (Fano plane lines)
    n_26 = 26  dim of F_4 fundamental (traceless Albert algebra)

Reference:
    Kahsay, Kibrom Kidane (2026). The Octavian Coherence Gate: The Four
        Irreducible Integers of the E₈(1) ⊃ G₂(1) × F₄(1) Conformal Embedding.
        Zenodo. https://doi.org/10.5281/zenodo.20493955
    The chain from the four integers to the scorecard is this repository.

Run the derivation via:
    python3 interference/run.py

This package is the canonical statement of the theory.  Where a
number here differs from a number in the published papers, the number
here is canonical.  The module docstrings carry the papers' arguments,
so the package can be read on its own.

================================================================================
What this package is (and is not)
================================================================================

THE ONTOLOGY.  To be is, minimally, to be self-present.  Anything that
exists already affects itself by being there.  Interference is the
physical image of that self-presence.  Closure takes a record to be a
finite dimensional unital real composition algebra.  Composition is
bilinear, and its positive definite norm (the ledger's probability
weight) is multiplicative.  Hurwitz's theorem then leaves exactly R, C,
H, and O.  With R as the amplitude baseline, the three nonreal
extensions C, H, O are the three searchlight depths.  The next
Cayley-Dickson algebra fails because it has zero divisors.  This is a
complete theorem given the framework's ledger criterion.  It is not a
claim that Hurwitz's theorem by itself says anything about
self-reference.

Records, rather than Hilbert rays, are primitive.  A physical phase must
therefore be carried by a relation between closed records. Permutations,
braids, and interference paths provide such relations, while a scalar
acting identically on the complete vacuum does not add a physical line.
For gauge winding, the complete matter coupled character is bar(theta),
and its action factors through the closed vacuum record algebra.  This is
part of the ontology, not a consequence of local MTC data.  Applied to
A=1+B in couplings.py, it leaves only bar(theta)=0.

The recursive Web (implemented, with ancestry checks, in root.py) can
go deeper than the searchlight.  The octonionic searchlight finds the
a->b->c->a loop at depth three.  The Web refines that pattern at depth
four and beyond.  G_2 = Aut(O) labels what the searchlight can track.
It does not limit what reality can do.  The protected forgetting branch
([Delta, S] = 0, masses.py) is a three eigenvalue circulant, identified
with the three generations.  Mass is the squared gap of that pattern.
G_2 has trivial centre and trivial pi_1, so every quantum number is
emergent.  When G_2 breaks to SU(3), the Z_3 centre appears.

THE SEARCHLIGHT.  The octonions reduced blind interference fluctuations
to seven directions.  Through the unique (7,26) bridge they reduced the
search to a finite 7x26 one.  That is why the closing loop could be
found at all.  The smallest structure becomes observable at depth three,
the same depth the searchlight can resolve.  The Chladni plate analogy
holds as a method: the pattern first becomes visible at finite
resolution, and later recursive readings refine it.  The engine's
arithmetic never evaluates an octonion or a sedenion.  The algebra's
entire output is a handful of discrete selections, each with a closure
status in the engine:

  selection              | engine constant | closure
  -----------------------+-----------------+--------------------------
  normed ledger gate     | searchlight = 3 | DERIVED given closure's
  + Hurwitz              |                 | unital positive composition-
                         |                 | algebra criterion. Not Web cap
  quantum dims of SU(3)_3| d_10=2, d_11=3  | the four integers
  Fano CG |C_3bar|/|C_1| | B/A = sqrt(2)   | IDENTIFIED: the Fano CG
                         |                 | ratio is computed, and the
                         |                 | triality amplitude rule is
                         |                 | identified (the classical
                         |                 | trace gives 5/12, excluded,
                         |                 | tests/probes/amplitude_gate.py)
  Sugawara weight h_10   | theta = 2/9     | IDENTIFIED: lightest
                         |                 | triality changer and theta=h
                         |                 | mass phase map. Sole positive
                         |                 | nondegenerate angle menu survivor
  Q_0 = 1/3 + d_10/6     | Q_0 = 2/3       | DERIVED downstream from B/A
                         |                 | by the circulant identity

The last row is the model for the rest.  The program (executable as
test_searchlight_absorption in test_coverage.py) is to derive every
searchlight selection from the engine, as Q_0 already is.  At that
point the algebra serves only as motivation.  The strong hypothesis
that the program tests is this: the octonions are not a separate
searchlight shone on the territory, but the territory's own accounting
seen from the ledger's side.  The absorption program is to derive each
selection this way.  Deeper self-interference is carried
by the recursive web, not by extending the Cayley-Dickson tower.  The
pattern becomes visible at depth three and gains finer structure as
later backreactions are evaluated.

DEPTH VOCABULARY. ``searchlight depth`` is the octonionic tracking
resolution (at most three). ``web depth`` is the recursive interference
depth stored by Ledger and EchoTerm. ``beta order`` is the imported
perturbative expansion used for transport and is not web depth. A module's
``layer`` is only its place in the dependency or physical regime narrative.
No value acquires a second physical identity by changing calculation or
comparison coordinates.

GATE VOCABULARY. Three sets of six carry the word gate. The lane's Gates
1 to 6 (below) are the conditions the closing lane meets. The coherence
gates (a) to (f) scan the E_8 embedding on Lie table invariants
(searchlight/embedding_uniqueness.py). The gravity lane's gates S, U, M,
C, T, P (sign, universality, magnitude, one cone, tensor carrier, no
preferred location) test gravity on that lane (gravity.py).

WHY THREE.  G_2 -> SU(3) supplies a Z_3 centre (Gate 2).  Given it, a
single sector changing insertion S first returns to its sector at
third order: tr S = tr S^2 = 0 and tr S^3 = 3 (masses.py, FIRST
INVARIANT ORDER THEOREM).  Orders one and two are forbidden, not
small, so within a sector the first closure is the three cycle
a->b->c->a.  The Web's bridge self-loop and the e<->q two cycle
(root.py, eq_loop) join the G_2 and F_4 sectors.  That bound is the
territory's.  The map's bound is Hurwitz: three is the deepest the
octonionic searchlight resolves.  The two bounds are independent, and
that they meet at three is why the closing loop could be found at all
(THE SEARCHLIGHT).

THE GENERATIVE DIRECTION.  The engine derives the scorecard.  The
substrate probes test whether the same structure emerges dynamically:
the self-trapped lump, circulation, common mode sourcing, confinement and
clock locking.  Each probe states its own verdict and kill
condition.  The circulant ratios are closed algebraically.
emergent_clock.py reports broadband clocks in substrate only 1D and 3D
runs.  A confined phase 3D run that locks to the algebraic circulant
belongs to the realization program.

WHERE THE QUESTIONS END.  The engine takes one number, and that
number is a choice of unit, not physics.  The ratio m_e/M_Pl is
derived outright.  G, read through the common mode, sits within one
standard deviation of the G measurements under an added spread or the
unweighted median, and
2.4 below CODATA's multiplicative combination (root.py,
tests/probes/g_measurements.py G6).  That
reading is an identification: the boundary is read as the full Planck
mass, and the common mode as the count of the lepton circulant.  So
the dimensionless output takes no input.  How strongly the substrate
interferes is not a second input either.  It is what is being derived,
one channel at a time: alpha_G_2(M_Pl), alpha(0), alpha_s.  What no
theory can answer is why the amplitude is nonzero at all.
Self-presence is UNIVERSAL, not conditional.  The empty world is
present to itself too, and because it is nothing, its self-interaction
generates no amplitude.  This is the zero orbit of x = b + W(x): zero
in, zero out.  It is not an exception to the axiom but its trivial
fixed point.  Nothing does not escape the axiom.  It simply never
closes into anything further.  Because self-interference cannot be
reduced further, anything beneath it is undecidable from inside, not
merely unanswered.

THE GAUGE FIELDS.  The gauge sector of the web is carried by collective
excitations of the confined phase, not postulated fields.  The algebraic
structure (representations, couplings, confinement) is derived from the
E_8/MTC data.  The realization program is to reproduce the perturbative
gauge series on composite knots (g-2, scattering amplitudes, decay
widths) and the multiknot Fock/Pauli structure beyond two body exchange.

THE SUBSTRATE STRENGTH.  Closure makes the emergent coupling a
function of the substrate parameters.  The derived alpha_G2(M_Pl) =
1/(24 pi), read as the medium's exchange coupling, places the depth at
rho xi^3 = 3 or 4.25 in the range convention
(tests/audit/parameter_ledger.py P9).  With the
depth read this way, the only contingent inputs are existence itself
and the choice of units.

CLOSURE BY DERIVATION PATH LENGTH. Partial closures exist.  Some lanes
make electron like objects whose ratios do not carry on to the next
sector.  Others make quark like objects that never bind.  These lanes
are eliminated because they stop.  The survivor is the lane with the
longest derivation path: each closure enables the next, from the first
loop to chemistry.  The seven branch gauntlet is the computed record
of the shorter paths.

This is NOT a top down Lagrangian theory of the form "write an action,
solve it, get particles."  It is a CLOSURE PRINCIPLE applied to a
substrate with blind self-interference.  The question is not "what is
the Lagrangian?" but "given self-interference, which patterns close?"
The answer is exactly one lane.

OBSERVATION AND CORROBORATION.  The framework observes the pattern.
The established theorems it cites (Witten quantization, anomaly
matching, Skyrmion-Dirac spectral flow, boundary condition changing
operator theory) describe the same regularities, reached from the
Lagrangian side.  A citation is therefore corroboration: two
independent records of one pattern, one of them written by people who
had never seen this engine.

THE GATES ARE CONSISTENCY CONDITIONS WITHIN OUR LANE.  Given a three
cycle closure, the level, the circulant, and the E_8 coherence follow.
At loop 3 the surviving lane is unique.  That is the computed content
of "exactly one lane" (the seven branch gauntlet below).  Loop 3
itself follows from the Z_3 of Gate 2 and the first invariant order
(WHY THREE above).  The anchor m_e
identifies which lane we are in.  The four integers are readings of
observed structure: they are counts.

THE HORIZON. Given the norm multiplicative ledger gate, n_7 = 7 is the
maximum number of imaginary directions the octonionic searchlight can
distinguish.  Anything beyond it is read as an amount with no
identity.  That is the gravity sector (PvP = 0, E[v^2] = 1/2 per mode,
with the 182 bridge channels as one ensemble whose members cannot be
told apart).  The dark ratio follows: f_baryonic = 1/(2 pi), so the
fraction that can be named is one part in a full turn.  The scorecard
is the complete table of things that can be named.  "No fourth
generation" is absolute, because a fourth would be a nameable
identity.  Unnamed structure is expected, and gravity is its only
visible effect.  This is also why the bridge reads as a scalar.
Tracking its 182 channels as one coherent multilayer object would need
a normed division algebra beyond O, and none EXISTS.  So the ensemble
is read as real scalars.  The vanishing coset charge
8 - 14/5 - 26/5 = 0 makes that decomposition exact.  Each channel's
heat kernel coefficient +1/6 is positive and makes gravity ATTRACTIVE
(boundary.py).  Reading the same channels as vectors would give -1/3
per state, reverse the sign, and matter would never bind.

ONE LOTTERY, THEN ARITHMETIC.  The closure question (which
configuration closes) is asked exactly once, at the first a->b->c
loop.  For every later sector the question is different: not which
configuration closes, but what can form.  That is already settled.
Any configuration that produces the leptons also produces these quarks
and this gravity, because they come from the SAME E_8 -> G_2 x F_4
structure.  The universe does not have to get lucky twice.  This is
why the mass sector is one generator plus a grammar.

The dynamical equation is the hermitian exchange form of
tests/substrate/nls_soliton.py, the Z_3-symmetric NLS of three
components,

    i hbar d_t psi_k = -D nabla^2 psi_k + dV/d psi*_k,   k = 1, 2, 3,
    V = (g_0/2) R^2 + g_1 sum_k rho_k Re(psi*_{k+1} psi_{k+2}) + (kappa/3) R^3,

R the total density, with |g_1/g_0| = 1/sqrt(2) and the self-bound
sign g_0 < 0 < g_1 (tests/probes/self_bound_vacuum.py).  The common
term binds, and the g_1 term locks the three phases to the clock
(tests/probes/lump_stability.py).  In density density notation the
same couplings are

    i hbar d_t psi_k = -D nabla^2 psi_k + sum_j g_{kj} |psi_j|^2 psi_k,

with g_{kk} = g_0 and g_{kj} = g_1 for k != j.  The phase lock comes
from the exchange term.

The Fano plane cross product gives |C_3bar|/|C_1| = sqrt(2)
(searchlight/octonions.py).  Gate 3 selects the SU(3)_3 modular data.
Their canonical positive spherical trace gives the channel weight
1/d(fund) = 1/2.  The channel to circulant identification combines that
trace with two cycle orientations and reads B/A = 2 sqrt(2)/2 = sqrt(2).
The substrate realizes the same ratio as its coupling, |g_1/g_0| =
1/(B/A) (tests/substrate/nls_soliton.py,
tests/probes/circulant_normalization.py).  The substrate is a
realization, and no engine number depends on it.  A classical Q = 2/3 on
the knot forming branch would falsify that separation.  The Madelung
decomposition psi_k = sqrt(rho_k) exp(i S_k / hbar) is exact, and the
Z_3 Fourier split into q = 0 (geometry) and q = 1, 2 (matter) is exact.
This resolves the Barcelo-Liberati-Visser "two roles of rho" obstruction
in single component analog gravity.  Status of the equation itself: the
NLS form is DERIVED, up to locality, as the universal envelope normal
form of weak local self-interference
(tests/probes/envelope_universality.py).  Locality is an axiom of
self-presence.  A nonlocal self term would be interference with
something elsewhere: a second thing.  That is coupling, not
self-interference.

From this equation:
  * the BdG linearisation gives Delta_k = A + B cos(theta + 2 pi k / 3),
    forced by [Delta, S] = 0.
  * in the empty field the equation holds self-trapped lumps with the
    clock inside, whose weight is their energy
    (tests/probes/lump_in_empty_field.py).
  * the gradient bilinear summed over the three components is the
    metric of the carrier written without a background, which supplies
    the geometry channel for the induced gravity ledger
    (tests/probes/composite_metric.py).

================================================================================
The six gates that select the unique lane
================================================================================

Gate 1, the normed ledger criterion plus Hurwitz classification.
    Closure requires a finite dimensional unital real composition algebra:
    bilinear composition and a positive definite multiplicative norm.
    Hurwitz therefore leaves R, C, H, O, while the next Cayley-Dickson
    algebra has zero divisors. Relative to the real amplitude baseline,
    C, H, O are exactly three added record depths. G_2 = Aut(O).
    Code: root.py (n_7 = dim G_2 fundamental, Fano plane lines).

Gate 2, Centre emergence.
    G_2 has trivial centre Z(G_2) = {1}.  SU(3) has centre Z_3.  When
    G_2 -> SU(3), the Z_3 crystallises as a child only label.  Three
    sectors -> three generations.  Pure group theory.
    Code: masses.py (Z_3 circulant gap matrix -> three eigenvalues).

Gate 3, Level.
    Requiring the centre current J = (k, 0) to be a dimension one boson
    fixes h(J) = 1, which uniquely forces k = 3.  Hence SU(3)_3 is the
    unique WZW level at this layer.  In unified notation k = d_11 and
    the altitude K = k + h_dual(SU(3)) = d_10 * d_11 = 6.
    Code: root.py (k_SU3 = d_11, K = d_10 * d_11).

Gate 4, Closure.
    The cyclic branch [Delta, S] = 0 is the branch whose mass observable
    closes, fully determined by the layer's algebraic data with zero
    tunable inputs.  The direct finite-G_2 path retains the harmonic
    v = [1, -1/2, -1/2] and does NOT close.  Closure is the survivor
    filter.
    Code: masses.py (Brannen circulant), gravity.py (protected
    forgetting: P v P = 0; per mode E[v^2] = 1/2, total
    P v^2 P = P.  The two contractions are distinct and easily
    conflated).

Gate 5, Octavian / E_8 coherence.
    The common mode fibre is bosonic (integer spin venting), local
    (self-interference), and self-dual (the substrate is its own
    observer).  The central charge c = 8 is not an additional
    assumption: it is DERIVED from the conformal embedding selected by
    the preceding gates, c(G_2) + c(F_4) = 14/5 + 26/5 = 8.  Given
    these, the even unimodular lattice theorem (Griess) forces E_8.
    Among the rank-8 candidates {I_8, A_8, D_8, E_8} only E_8 passes
    all conditions. The rank-16 even unimodular lattices E_8 + E_8 and
    D_16+ enter through the minimal rank axiom (wiring.py stage W1
    computes the residuals).
    The bridge count dim(7, 26) = n_7 * n_26 = 182 follows.
    Code: wiring.py (lattice scan), root.py (c_coset = 0),
    gravity.py (bridge identification).

Gate 6, Three gate closure.
    Gates L/Q/G are the lepton residue, the quark bridge, and the
    gravity sign.  The last column is G_ind/G_N, UV/broken, as in
    Table V of the gravity record (The Echo of Standing Waves):
        localize G_2 first          -  Q  -   0.935/0.946
        PROTECTED FORGETTING ONCE   L  Q  G   0.994/1.007  <- closes
        full alpha from v^2         -  Q  -   1.062/1.077
        independent F_4 half        L  -  -   0.964/0.976
        dim localize 7 first        -  -  -   0.976/0.989
        gauged E_8 vector           L  -  -   repulsive
        massive E_8 Proca           L  -  -   repulsive
    Only "protected G_2 forgetting applied exactly once" closes L and
    Q, and its ledger passes the G sign gate.  (Seven branch audit:
    the gravity record, Table V, reproduced above.  The winning
    branch's venting ledger is executable as Sigma in gravity.py.)

WHAT THE GATES LEAVE.  One structural lane closes through the gates.
The closure search holds it fixed and runs the committed mass grammar
as a generator.  It prints the look elsewhere denominator.  The
canonical grammar is one of 504 mass grammar hits within the
scorecard's 1% band, out of 2,370,816 eligible realizations.  The data
fix the word counts uniquely but not the subpercent echo forms.  The
joint 0.1% tier is empty, the canonical grammar included, because the
u and b references are coarser than that.  A candidate that does not
close describes nothing that exists: what exists is what carries
non-zero interference weight.  Code: tests/audit/closure_search.py,
which prints the full accounting and states this in its docstring.

================================================================================
The closure
================================================================================

    14  +  52  +  182  =  248  =  dim(E_8)
    gauge   matter   bridge
            + Higgs  -> gravity

Every degree of freedom of E_8 has a physical job.  Nothing is left
over.  The layer data fix every number in the implemented scorecard.
The mixing sector is identified rather than derived.

The precision hierarchy follows the dependency graph. Leptons sit
closest to the Z_3 source. Quarks cross additional triality,
generation scaling, and confinement edges. Gravity first appears at
web depth four. These graph distances are not perturbative beta orders.

================================================================================
Coupling chain: methodological consequence
================================================================================

The framework has no dimensionless tunable inputs.  The electron mass
sets only the unit system: the first closure sets the
units.  M_Pl, G, m_mu, m_tau, v_EW are outputs of the FORCED
chain.
Every coupling, mass ratio, mixing angle, threshold, and Newton
normalisation is a closed form expression in
(d_10, d_11, n_7, n_26) + pi.

This has a strict consequence for how the code reads. The framework has
no free interference strengths to tune.  Each closure echoes through channels
that already exist (the bridge self-loop, the e<->q cycle, the WZW
vents), and each echo's weight is a theorem about its channel.  A beta
function or matching formula may convert a completed reading into the
coordinate an external calculation needs.  It does not create another
framework state.

The unified expression makes this explicit.  There are no per layer
modules with separate physical motivations.  All named algebraic
quantities reduce to closed form polynomials in (d_10, d_11, n_7, n_26).
Among others:

    K        = d_10 * d_11                          (WZW altitude = 6)
    h_10     = d_10 / d_11**2                       (h(1,0) = 2/9)
    h_11     = 1 / d_10                             (h(1,1) = 1/2)
    Q_0      = d_10 / d_11                          (Koide = 2/3)
    sin2W    = d_11 / (d_10**2 + d_11**2)           (Weinberg = 3/13)
    N_bridge = n_7 * n_26                           (= 182)
    dim_E_8  = 2*n_7 + 2*n_26 + n_7*n_26            (= 248)
    h_bridge = d_10/(1+d_10**2) + d_10*d_11/(1+d_11**2)   (= 1)
    c_coset  = 0                                    (rational identity)

The compression itself makes the claim easier to falsify.  Four
integers plus pi plus m_e producing the SM scorecard at subpercent
residuals is a much tighter claim than "each layer closes
consistently."

Concrete example, alpha_em = pi/512 and alpha(0) = 1/137.035999050 are
two coordinates of one coupling. alpha_em = pi/2**9 is the zero edge
base of the embedding layer (alpha_G2(mu*) over the Singh ratio), used
inside the instanton vertex. alpha(0) is the zero momentum reading:
the bridge SELF-loop gives the factor (1 - 1/(2 pi)) (h = 1,
unit coupling D2_local = 1, c_coset = 0), and the
depth 4 recursive reading of the electron<->quark carrier (two
orientations, the S, S-dagger rule) adds -2 (alpha / 2 pi)**2, so
1/alpha(0) is the real root of

    x**3 = (512/pi) [ (1 - 1/(2 pi)) x**2 - 1/(2 pi**2) ]
         = 137.035999050,

which matches the Berkeley Cs recoil measurement (+0.13 sigma) in the
open >5-sigma Cs/Rb dispute.  The same (7,26) bridge, in the SAME
sector with the SAME conformal weight h_bridge = 1, also enters the
induced stiffness through 182 minimally coupled scalar venting
channels (boundary.py).  EM and induced gravity are two readings of
one bridge, not two sectors.

================================================================================
The echo law and THE KERNEL (root.py)
================================================================================

One scalar recursion, x <- b + W(x), solved to its joint fixed point
(Web.solve(), classical fixed point iteration).  Every quantity is a
LEDGER: a base amplitude plus a stack of echoes, each with its source.
The reason is that whatever does not cancel locally vents into the
common mode, and later closures echo back through everything that
already exists.  Every weight names its justification.  Every
identified link is stated next to the physics that produces it,
together with the exact quantities that use it.  The searchlight
resolution is 3, by the norm multiplicative ledger criterion combined
with the Hurwitz classification.  Recursive web depth is not capped.
Every later edge must name and read its depth-(n-1) parent, and must be
derived and registered before any comparison with data.  The
multiplicity rules (orientation, sector node, vertex composition,
apportionment) and their uniqueness enumerations live with the law
itself in root.py.  FORCED terms enter canonical values.  CANDIDATE
terms are excluded.

================================================================================
File to role map
================================================================================

The directory groups are the four layers.  Each derivation module
exposes one derive(...) -> dict, and run.py calls them in dependency
order.

File                    | Role
------------------------+--------------------------------------------
run.py                  | End to end driver, scorecard, agreement
                        | summary.
root.py                 | Depth-3 searchlight vocabulary.  Every derived
                        | quantity as a polynomial / rational function
                        | of them.  The echo web and its depth gate.
masses.py               | 1 anchor + 8 mass predictions (2 leptons + 6 quarks)
                        | from the Z_3 circulant + WZW emergence
                        | ladder.
mixing.py               | WZW/D^(6)* connection and mixing candidates;
                        | the sector operator lift is depth 4
                        | (mixing_z3_dichotomy.py, mixing_operator_search.py).
couplings.py            | alpha_s(M_Z) from the embedding index chain.
                        | alpha(0) from (7,26) bridge self-interference
                        | (one bridge, two readings).
boundary.py             | The one list read by the running of alpha_G_2
                        | and by the stiffness.
gravity.py              | The gravity lane's gates on the universal venting
                        | ledger (every massive closure feeds the common
                        | mode). G is the anchor inversion.
                        | eta_B from G_2 instanton.  Neutrinos
                        | from the F_4-singlet RH neutrino identification:
                        | conditional minimal-type-I rank at most 2;
                        | ordering not selected.
higgs_ew.py             | SM RGE transport, Higgs quartic boundary
                        | (F_4(1) fusion + bridge), M_W, M_Z, import
                        | audit.  Engine import inventory.
dark_sector.py          | Omega_DM/Omega_b = 2 pi - 1. Friedmann
                        | critical density versus current dark energy
                        | density, CKN status explicit.
words.py                | Generation word lemma, executable and fully
                        | derived: base(n) = (4, 12, 97) as boundary-
                        | walk counts.  N-ality superselection (neutral
                        | lane).  Knot definition + conversion lemma
                        | (terminus).  Exclusions are algebraic (the
                        | grading and superselection lemmas). The data
                        | confirms them separately.
wiring.py               | The wiring: a closure residual scan over the
                        | classified candidate menus of lattices,
                        | conformal pairs, and levels.  The four
                        | integers are never an input.
reference.py            | Comparison data, single home, quarantined by
                        | tests/audit/test_blindness.py.
searchlight/octonions.py| Constructive Fano plane proof of
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
                        | from the common mode.  Weyl orbit uniqueness.

The SUBSTRATE and the AUDIT instruments read from the scrutiny
side, so they live with the tests:

tests/substrate/nls_soliton.py
                        | The substrate equation (Z_3 coupled NLS),
                        | split step evolution, BdG, soliton stability,
                        | Madelung sourcing.
tests/substrate/substrate3d.py
                        | The same equation in 3D, shared by the
                        | dynamics probes.
tests/audit/edge_grammar.py
                        | The committed grammar: menus, echo terms,
                        | the null model's alphabet.
tests/audit/closure_search.py
                        | The grammar as a generator: blind gates
                        | first, data contact once, look elsewhere
                        | denominator printed.

Checking whether the theory is flawless is testing, and the suite is
at the repo root in tests/, one folder per branch.  The
test_<branch>_scorecard.py files freeze the NUMBERS (scorecard, freeze
tables, attractor).  The test_<branch>_probes.py files freeze the
MECHANISMS, one conclusion per probe, with the substrate-
dynamics arc behind pytest -m slow.  test_coverage.py freezes the
EPISTEMICS (menu sparsity, cross bracing, citation resolution, the
searchlight absorption scoreboard).  test_blindness.py is the
reference data firewall.  The analyses themselves live in
tests/probes/, each carrying its own docstring and its own
verdict. That inventory is not restated here, because a list of
probes is the first thing to rot.

Quickstart:   python3 interference/run.py   (full chain + scorecard)
              python3 -m pytest -q          (the test suite)
Dependencies: Python 3.10+, numpy, scipy. pytest for the suite.
Ledger inspection: any node's interference history prints via
root.WEB[name].table() after the chain has run.

================================================================================
The bets, the coordinate, and the closed chain
================================================================================

Sole dimensional input: the electron mass (CODATA, +-0.3 ppb).  M_Pl
and G are OUTPUTS.  The depth 4 response refines the unit.  The EM
sector graph admits no further independent edge after it, because the
fixed point already resums repeated traversals of that carrier.  This
is channel finality, not a global depth cap.

Every number lives in ONE place: the module that computes it, with its
kill condition stated next to it.  The frozen values the code must
reproduce are the branch scorecards,
tests/<branch>/test_<branch>_scorecard.py.  Mixing is four readings of
the solved D(6)* boundary, one per face (mixing_z3_dichotomy.py).  The
substrate suite contains the lump in the empty field, clock runs,
confined phase lattice probes, the composite metric, and the
Skyrmion-Dirac spectral flow probe.

THE MASS READING. The algebra produces one completed Web reading per
fermion.  Standard experimental tables attach conventional comparison
coordinates to their entries.  Those labels are metadata for the
comparison, not extra stages of the derivation.  Full treatment:
masses.py.

THE CHAIN IS CLOSED.  The framework has no Lagrangian.  It observes the
interference, seeded by the anchor at depth three, so every reading is
an IDENTIFICATION: a chain object observed as a physical one, with the
exhaustion of its committed menu and the reason only this lane closes,
each separately falsifiable.  Within that, many links are DERIVED
outright (the algebraic chain) or ENUMERATED (the words lemma, K^3
closure, metric uniqueness).  The anchor sets the unit.  Each CITED
IMPORT WITH A STATED BAND is transport: it carries an engine reading to
the coordinate a measurement is quoted in, and no engine constant reads
it (tests/audit/test_blindness.py).  Every import is decomposed and
audited (higgs_ew.py, tests/probes/ew_internal.py).  Data driven
ingredients are measurements.  Parameter dependencies are recalculated
from the framework's own values and checked for compatibility.  Pure
loop integrals are mathematics fixed by the field content.  Every
nontrivial factor is a backreaction through an existing channel.
Nothing is tuned.

FORWARD RULE. The searchlight closes at depth 3.  The recursive web is
open to later resolution but closed against hindsight.  Edges run only
along the sector graph.  Multiplicities come only from the orientation,
sector node, vertex composition and apportionment rules.  Every new
depth is registered with its derivation before comparison.  As a
result every canonical residual is BOUNDED: a deviation larger than the
next admissible echo unit of its sector falsifies the rule that
produced it.  A deeper term off the committed menu cannot be invented
to absorb a miss.

COMPUTATION ORDER.  The six gates select SU(3)_3.  The lightest
triality changing primary has h = 2/9, and the identified mass phase map
sets theta = h.  The identified channel to circulant map applies the Fano
CG, canonical positive spherical trace and orientation count to obtain B/A.
Positivity and nondegeneracy leave theta = 2/9 as the sole value in the
committed structural menu.  The substrate realizes the amplitude ratio
as its coupling.  masses.py then compares the
calculated mass ratios with the PDG entries.

NO SEARCHED FORMULAS. No formula was constructed for an unmeasured
constant by searching combinations of (2, 3, 7, 26, pi).  In a four
integer grammar such searches always succeed, which is exactly why
their successes would carry no weight.  Every testable prediction in
this repository uses only edges already in the graph.  The framework
survives only if ALL of them hold: by its own rules, one clean failure
ends it. The look elsewhere denominator is printed, not asserted
(tests/audit/closure_search.py).
"""

# The derivation is a regular Python package. Sibling modules use
# package relative imports so stateful modules such as ``root`` load only
# as ``interference.root`` and the solved WEB singleton remains unique.
