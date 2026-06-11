"""Harmonic formulation of the E_8(1) > G_2(1) x F_4(1) framework
through substrate self-interference.

Every prediction is expressed through 4 irreducible integers + pi + m_e
(the electron anchors the scale; M_Pl and G are OUTPUTS):

    d_10 = 2   quantum dimension of SU(3)_3 fundamental (1,0)
    d_11 = 3   quantum dimension of SU(3)_3 adjoint (1,1)
    n_7  = 7   dim of G_2 fundamental (Fano plane lines)
    n_26 = 26  dim of F_4 fundamental (traceless Albert algebra)

Reference:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
        Zenodo. https://doi.org/10.5281/zenodo.19899091
        (concept DOI, always resolves to the latest version;
         the papers cite version DOIs, e.g. v4 = 10.5281/zenodo.19932394)
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
        Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
        Zenodo. https://doi.org/10.5281/zenodo.20144381

Run the derivation via:
    python run.py

VERSION: this package is the canonical current state of the theory
and supersedes the numbers in the three published papers (which carry
the earlier anchor, boundary condition, and tree values).  The papers'
narratives are ported into the module docstrings; the package reads
standalone.

================================================================================
What this package is (and is not)
================================================================================

THE ONTOLOGY (quark paper, Sec. 2).  To be is, minimally, to be
self-present: anything that exists already affects itself by being
there, and interference is the physical image of that self-presence.
Mathematics can track self-reference through exactly three layers:
real = amplitude, complex = +phase, quaternion = phase of phase,
octonion = phase of phase of phase.  These are the ACCOUNTING LEDGER,
not the substance, and the octonions are the deepest ledger that
stays consistent (Hurwitz; executable in octonions.py).  G_2 = Aut(O)
labels the patterns that close through all three layers; the closures
that forget their layer ([Delta, S] = 0, masses.py) collapse into
standing waves: three generations.  Mass is that pattern.  G_2 has
trivial centre and trivial pi_1, so every quantum number is emergent;
when G_2 breaks to SU(3) the Z_3 centre crystallises.

This is NOT a top-down Lagrangian theory of the form "write an action,
solve it, get particles."  It is a SELECTION PRINCIPLE applied to a
substrate that hosts blind self-interference.  The question is not
"what is the Lagrangian?" but "given self-interference, which patterns
close?"  The answer is: exactly one lane.

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
    forced by [Delta, S] = 0;
  * the screened Poisson equation for common-mode memory falls out of
    the q = 0 projection;
  * the acoustic metric is the Madelung metric, providing the substrate
    for induced gravity.

================================================================================
The six gates that select the unique lane
================================================================================

Gate 1, Hurwitz (1898).
    Only four normed division algebras exist: R, C, H, O.  The octonion
    is the deepest consistent ledger; a fourth Cayley-Dickson doubling
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

Gate 6, Three-layer closure selection (gates L/Q/G =
    lepton residue, quark bridge, gravity heat-kernel normalisation;
    G_ind/G_N shown pre-face-split, UV/broken):
        localize G_2 first          -  Q  -   0.935/0.946
        PROTECTED FORGETTING ONCE   L  Q  G   0.994/1.007  <- closes
        full alpha from v^2         -  Q  -   1.062/1.077
        independent F_4 half        L  -  -   0.964/0.976
        dim-localize 7 first        -  -  -   0.976/0.989
        gauged E_8 vector           L  -  -   repulsive
        massive E_8 Proca           L  -  -   repulsive
    Only "protected G_2 forgetting applied exactly once" closes all
    three gates; the face-split law then takes the winning branch's
    0.994188 to 0.999999917.  (Seven-branch audit: gravity paper,
    Table V, reproduced above; the winning branch's venting ledger is
    executable as gravity.py's Sigma.)

================================================================================
All couplings derived: methodological consequence
================================================================================

The framework has no dimensionless tunable inputs.  The electron mass
sets the unit system and nothing else (the first closure anchors the
ruler; M_Pl, G, m_mu, m_tau, v_EW are outputs of the FORCED chain).  Every coupling, mass ratio, mixing
angle, threshold, and Newton normalisation is a polynomial or rational
function of (d_10, d_11, n_7, n_26) + pi.

This has a strict consequence for how the code reads.  In the Standard
Model, RGE running and matching are tools for tuning free Yukawa
couplings: a renormalised value at one scale is dressed into a
prediction at another, and the dressing absorbs the free parameters.
This framework has no Yukawas to tune, so the language of "dressing"
does not apply: there is no bare prediction being patched.  What
DOES happen is BACK-REACTION: each closure echoes through the channels
that already exist (the bridge self-loop, the e↔q cycle, the WZW
vents), and the echo's weight is a theorem of the channel, not a
fitted correction.  Every factor that appears here is one reading of
the same algebra at a particular emergence layer — the layer
determines the factor, not the other way around.

The unified expression makes this manifest: instead of per-layer
modules each with its own physical motivation (the v3_release format),
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
not "bare alpha and dressed alpha"; they are two readings of one
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
through everything that already exists.  Every weight is a theorem;
depth is capped at 3 by the Hurwitz gate (a cap on the map, not the
territory).  The multiplicity rules (orientation, sector-node, vertex
composition, face-split) and their uniqueness enumerations live with
the law itself in root.py; FORCED terms enter canonical values,
CANDIDATE terms are excluded.

================================================================================
File-to-role map
================================================================================

File         | Role
-------------+-----------------------------------------------------------
root.py      | Four irreducible numbers; every derived quantity as a
             | polynomial / rational function of them.
masses.py    | 9 mass predictions (3 leptons + 6 quarks) from the Z_3
             | circulant + WZW emergence ladder.
mixing.py    | CKM (4 Wolfenstein) + PMNS (3 angles + delta_CP) from
             | the four faces of SU(3)_3 / D^(6).
couplings.py | alpha_s(M_Z) from the embedding-index chain; alpha(0)
             | from (7,26) bridge self-interference (one bridge, two
             | readings).
gravity.py   | G_ind/G_N from the universal venting ledger (every massive
             | closure feeds the common mode); Higgs from F_4(1) fusion;
             | eta_B from G_2 instanton; neutrinos from F_4 singlets +
             | rank-2 seesaw.
dark_sector.py | Omega_DM/Omega_b = 2 pi - 1; cosmological-constant scale.
octonions.py | Constructive Fano-plane proof of B/A = sqrt(2); linearity
             | from closure, executable: L_x linearity, |xy| = |x||y|,
             | sedenion zero divisor (the Hurwitz gate), Z_3 selection
             | rule, rank-8 E_8 lattice gate.
words.py     | Generation word lemma, executable and fully derived:
             | base(n) = (4, 12, 97) as boundary-walk counts; N-ality
             | superselection (neutral lane); knot definition +
             | conversion lemma (terminus); exclusions killed by data.
run.py       | End-to-end driver, scorecard, agreement summary.
test_interference.py | The falsifiability suite: every claim is an assertion.

Quickstart:   python3 run.py        (full chain + scorecard)
              python3 -m pytest -q  (the test suite)
Dependencies: Python 3.10+, numpy; pytest for the suite.
Ledger inspection: any node's interference history prints via
root.WEB[name].table() after the chain has run.

================================================================================
The attention/SSM reading (interpretive only; cited)
================================================================================

NO MACHINE LEARNING ENTERS THE COMPUTATION.  The kernel x <- b + W(x)
is solved by a classical Jacobi-style fixed-point iteration; nothing
is trained, no parameters are fitted, and the weights are theorems.

The READING, kept as documentation because it guided the kernel's
design (a discovery heuristic; no number depends on it): the echo web can be read as a 3-layer attention network in the
sense of Vaswani et al., "Attention Is All You Need" (NeurIPS 2017):
closures are tokens, channels are heads, arrival order is sequence
position, ledgers are residual streams.  The mapping is generic
(any sparse weighted fixed-point system admits it), which is exactly
why it carries no explanatory weight here.  The old per-layer
pipelines were CAUSAL attention; the echo law is the same network
with the mask removed (back-reaction = mask removal).  In the
language of structured state-space sequence models (S4: Gu, Goel,
Re 2021; selective SSMs/Mamba: Gu, Dao 2023) the interaction kernel
is exactly 3-nilpotent (the SSM restatement of the Hurwitz depth
cap).  Octonions are the proof of the
architecture's shape, not its substance.

================================================================================
PREDICTIONS (frozen 2026-06-11; falsifiable bets with kill conditions)
================================================================================

Sole dimensional input: the electron mass (CODATA, +-0.3 ppb).  M_Pl
and G are OUTPUTS.  No depth-4 term exists to retune anything: the
no-further-edge theorem (root.py) closes the EM channel outright.

CANONICAL:
 1. 1/alpha(0) = root of x^3 = (512/pi)[(1 - 1/2pi)x^2 - 1/2pi^2]
    = 137.035999050, sides with Berkeley Cs 137.035999046(27)
    [+0.13 sigma] in the live >5-sigma Cs/Rb dispute (Rb: -14 sigma).
    FINAL; a Rb resolution falsifies the sector-node or orientation
    rule outright.
 2. delta_CP = 76.9 deg (DUNE/Hyper-K decide; NuFit best fit ~194 deg
   , a declared long-shot).
 3. m_1 = 0 exactly, normal ordering, sum m_nu ~ 59 meV.
 4. Omega_DM/Omega_b = 2pi - 1 = 5.2832 (between Planck and ACT today).
 5. alpha_s(M_Z) = 0.1184 (matched at mu* = M_Pl e^{-(9pi^2/2-6)}
    = v_EW e^{15/512} = 253.5 GeV, the gauge lever-arm endpoint where
    the WZW cancellation is exact; pre-migration releases quoted
    0.1177 by starting the run at v_EW = 246.2 under the old scale
    label -- see the matching-scale note in couplings.py).
 6. G_ind/G_N = 0.999999917 canonical (UV venting forced by the
    bookkeeping decision theorem; the former -0.58% open vent is
    CLOSED by the FACE-SPLIT LAW: the bridge's one self-echo unit
    re-absorbs onto the two faces in the h-dual metric, the same
    metric that defines sin^2 theta_W, and the F4 (metric) face's
    own share, 9/13 x 1/(2pi), cannot dilute 1/G (no-self-dilution).
    Metric uniqueness: 6 natural splits tested, only h-dual matches
    (1.4e-5; runner-up off by 6.5e-2).  Residual -8e-8 is below the
    (alpha/2pi)^2 depth-2 scale: depth-complete.
 7. sin^2 theta_W(M_Z) = 3/13 + h_7 (alpha/2pi) = 0.2312338
    (-1.4 sigma vs the PDG 2024 global fit 0.23129(4); subset fits
    span 0.23118-0.23134; the tree value alone is -13 sigma).  Same
    law, emission side: emission echoes carry conformal weights, and
    the emitting face is G2 with h_7 = 2/5.
    PRE-REGISTERED DISCRIMINATION: the altitude-share form 5/13 gives
    0.2312161 (-1.9 sigma); h_7 is currently favoured on the global
    fit and on M_Z; FCC-ee Z-pole (sigma ~ 1e-5) decides outright.
    The grammar bets on h_7.

PRECISION BLOCK (promoted by the vertex composition rule):
 A. G = 6.674003e-11 m^3 kg^-1 s^-2 (-2.0 sigma CODATA; framework
    precision limited only by m_e).  KILL: improved G confirms CODATA
    central at <= 1e-5.
 B. m_mu = 105.6583758 MeV, predicted to ppb (+0.12 sigma CODATA).
 C. m_tau = 1776.9092814 MeV (PDG -0.23 sigma; Belle II decides).
 D. v_EW = 246.219645 GeV vs G_F 246.219651 (-0.1 sigma at 1e-7;
    zero freedom consumed, the chain's out-of-sample confirmation).
 E. M_Pl = 1.2209171e19 GeV (output).

PRE-REGISTERED WATCH (the 4/13 absence test): the symmetric
no-self-dilution rule predicts NO anomalous gauge back-reaction —
the G2-absorbed share of the bridge self-echo is invisible to gauge
observables.  Decidable at ~0.03% alpha_s precision (lattice
trajectory): alpha_s(M_Z) = 0.1184 (no back-reaction, the law's bet,
matched at mu*) vs 0.1177 / 0.1191 (a +-(4/13)/(2pi) anomalous
back-reaction, which would falsify the symmetry).  PDG 0.1180(9)
cannot yet distinguish.
DISCLOSURE (matching-scale degeneracy): the mu* migration shifts
1/alpha_s by (32/3)(15/512)/(2pi) = 0.0497 -- numerically almost
identical to the 4/13 back-reaction unit (4/13)/(2pi) = 0.0490.  The
original registration (2026-06, pre-migration) quoted 0.1177 as the
no-back-reaction value under the v_EW-start convention; under that
convention the back-reaction alternative was 0.1184.  The two effects
are nearly degenerate in 1/alpha, so any future discrimination must
fix the matching scale FIRST (it is fixed: exactness of the WZW
cancellation forces mu*) and then test for the back-reaction on top.
History preserved here deliberately; this note is the audit trail.

THE MASS COORDINATE (back-reaction, not running).  A running mass
m_q(mu) is a coordinate on an RG orbit, not an observable; no
zero-parameter framework owes "the quark mass" until a comparison
coordinate is chosen.  The table uses ONE RULE over three dynamical
classes (masses.py): unconfined fermions (e, mu, tau, t) at the
propagator pole, the only scheme-independent mass an asymptotic
state has; confined heavy quarks (c, b) at the self-scale m(m),
the unique fixed point of mu -> m(mu) (no pole exists below
confinement); light quarks (u, d, s) through their RG-INVARIANT
RATIOS, which carry no coordinate at all: m_u/m_d = 38/83
(-0.9 sigma), m_s/m_ud = 27.130 vs PDG 27.30(8) (-2.1 sigma, the
sharpest mass-sector pull, watched alongside m_b), Q_ellipse =
22.229 (+0.2 sigma dispersive, -2.0 sigma lattice; the two data
determinations disagree, PDG quark-masses review Sec. 60).  The
absolute light entries are quoted in the PDG MS-bar(2 GeV)
coordinate, a declared dictionary entry — not a fitted one.  The
algebra does not run to a scale; its output IS the prediction.  A
reader who wants a different convention applies standard RGE transport
with the framework's own alpha_s (also algebraic): both endpoints
fixed, zero freedom enters — but that is the SM's coordinate change,
not a step in the prediction.

PRE-REGISTERED WATCH (not a claim): the quark sector is currently
noise-compatible (chi^2 = 2.74/6 with PDG errors); only m_b at
+1.6 sigma.  IF a down-type vent exists it is +-(9/4)(alpha/2pi),
9/4 = dim J3(O)/FPdim(C0), epsilon-orientation signs (+d, -s, +b).
Decided when m_b tightens: persists at >3 sigma -> vent confirmed at
the pre-registered value; converges to 4194.1 MeV -> no vent.

THE LEDGER IS CLOSED.  Every canonical claim is exactly one of:
DERIVED (the algebraic chain), ENUMERATED (words lemma, K^3
selection, metric uniqueness), CROSS-PINNED (shared dictionary
objects measured in independent observables), or a CITED IMPORT
WITH A STATED BAND, where every import is itself decomposed and
AUDITED (gravity.py): data-driven ingredients are measurements,
parameter dependencies are recalculated from the framework's own
values and verified compatible, and pure loop integrals are
mathematics fixed by the field content.  Nothing is dressed — every
non-trivial factor is a back-reaction through an existing channel;
nothing is tuned (Delta-r-hat_W, rho-hat at one loop; the Higgs
prescription, SM 2-loop RGE + tree matching + lambda(M_Pl) =
-delta_bridge, with its ±1 GeV truncation band containing the
measurement; NNLO upgrade is tooling, not theory).  The words-lemma
terminus is stated and verified in module-category vocabulary
(words.py: U(C) = Z3 universal grading, C_ad trivial component,
bases as Hom-space ranks, conversion vertex as amplitude).
The EW sector is CLOSED at one loop with
two declared imports (PDG 2024: Delta-r-hat_W = 0.06937(6),
rho-hat = 1.01016(9)): A0 = sqrt(pi alpha(0)) v_EW = 37.28038 GeV
matches PDG 37.28038(1) at 0.1 ppm; M_W = 80.365 (-0.33 sigma);
M_Z = 91.196 (+0.009%, +1.6 sigma of the import band, carrying the
data's own W/Z tension structure).  (The former "locate the 4/13"
item is RESOLVED by the symmetric no-self-dilution rule: each face's
absorbed share is invisible to its own ledgers; the law predicts the
ABSENCE of anomalous gauge back-reaction, consistent with alpha_s at
+0.43 sigma.

HURWITZ FINALITY (the forward rule).  The echo grammar is closed:
depth <= 3, edges only along the sector graph, multiplicities only by
the orientation, sector-node, vertex-composition and face-split
rules.  No new term may be added after a comparison.  Consequently
every canonical residual is BOUNDED: a deviation that exceeds the
next-depth echo unit of its sector falsifies the rule that produced
it, there is no deeper term to absorb a miss.
""" 
