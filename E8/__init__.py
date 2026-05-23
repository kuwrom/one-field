"""E8 -- modular E_8(1) > G_2(1) x F_4(1) derivation chain.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
        Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
        Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
        Zenodo. https://doi.org/10.5281/zenodo.20144381

Run the full derivation via:
    python -m E8

================================================================================
What this package is (and is not)
================================================================================

This is NOT a top-down Lagrangian theory of the form "write an action,
solve it, get particles."  It is a SELECTION PRINCIPLE applied to a
substrate that hosts blind self-interference.  The question is not
"what is the Lagrangian?" but "given self-interference, which patterns
close?"  The answer is: exactly one lane.

The dynamical equation IS written down.  It is Eq. (2.1) of the
gravity paper ("The Echo of Standing Waves"): the Z_3-symmetric
coupled NLS

    i hbar d_t psi_k = -D nabla^2 psi_k + sum_j g_{kj} |psi_j|^2 psi_k,
                                                       k = 1, 2, 3,

with g_{kk} = g_0 and g_{kj} = g_1 = g_0 / sqrt(2) for k != j (the G_2
constraint).  The Madelung decomposition psi_k = sqrt(rho_k) exp(i S_k
/ hbar) is exact, and the Z_3 Fourier split into q = 0 (geometry) and
q = 1, 2 (matter) is exact.  This resolves the Barcelo-Liberati-Visser
"two roles of rho" obstruction in single-component analog gravity.

From this equation:
  * the BdG linearisation gives Delta_k = A + B cos(theta + 2 pi k / 3)
    -- forced by [Delta, S] = 0;
  * the screened Poisson equation for common-mode memory falls out of
    the q = 0 projection;
  * the acoustic metric is the Madelung metric, providing the substrate
    for induced gravity.

================================================================================
The six gates that select the unique lane
================================================================================

Gate 1 -- Hurwitz (1898).
    Only four normed division algebras exist: R, C, H, O.  The octonion
    is the deepest consistent ledger; a fourth Cayley-Dickson doubling
    introduces zero divisors and breaks bookkeeping.  Its automorphism
    group is G_2 = Aut(O).  Theorem, not postulate.
    Code: algebra.py, octonions.py.

Gate 2 -- Centre emergence.
    G_2 has trivial centre Z(G_2) = {1}.  SU(3) has centre Z_3.  When
    G_2 -> SU(3), the Z_3 crystallises as a child-only label.  Three
    sectors -> three generations.  Pure group theory.
    Code: algebra.py (Sec. 0.4), leptons.py.

Gate 3 -- Level.
    Requiring the centre current J = (k, 0) to be a dimension-one boson
    fixes h(J) = 1, which uniquely forces k = 3.  Hence SU(3)_3 is the
    unique WZW level at this layer.  All conformal weights, quantum
    dimensions, and altitudes used downstream are fixed by k = 3 alone.
    Code: leptons.py (theta = h(3) = 2/9), quarks.py (WZW emergence).

Gate 4 -- Closure.
    The cyclic branch [Delta, S] = 0 is the one whose mass observable
    closes -- fully determined by the layer's algebraic data with zero
    tunable inputs.  The direct finite-G_2 path retains the harmonic
    v = [1, -1/2, -1/2] and does NOT close.  Closure is the survivor
    filter, not an assumption.
    Code: leptons.py (Brannen formula), gravity.py (protected forgetting
    P v P = 0, P v^2 P = (1/2) P).

Gate 5 -- Octavian / E_8 coherence.
    If the common-mode fibre is bosonic, local, self-dual, and c = 8,
    the even-unimodular lattice theorem (Griess) forces E_8.  Among
    rank-8 candidates {I_8, A_8, D_8, E_8, E_8 + E_8}, only E_8 passes
    all four conditions.  The bridge count dim(7, 26) = 182 is then
    fixed by the unique exceptional conformal embedding E_8(1) ⊃
    G_2(1) × F_4(1).
    Code: algebra.py (central charge sum rule c(G_2) + c(F_4) = 8),
    gravity.py (bridge identification and five-argument scalar reading).

Gate 6 -- Three-paper closure audit.
    Seven competing branches are tested across the lepton, quark, and
    gravity layers.  Only "protected G_2 forgetting applied exactly
    once" closes all three.  Gauging E_8 makes gravity repulsive;
    localising G_2 first leaves finite bias; full alpha overshoots;
    independent F_4 half undershoots.  The other five fail specific
    gates with specific numbers.  See gravity paper Table V.
    Code: gravity.py (Sec. 10.7), scorecard.py.

================================================================================
All couplings derived: methodological consequence
================================================================================

The framework has no dimensionless tunable inputs.  M_Pl sets the
unit system and nothing else.  Every coupling, mass ratio, mixing
angle, threshold, and Newton normalisation is an algebraic output of
the six gates above.

This has a strict consequence for how the code reads.  In the
Standard Model, RGE running and matching are tools for tuning free
Yukawa couplings against measurements -- a renormalised value at one
scale is dressed into a prediction at another.  This framework has no
Yukawas to tune, so the language of "correction" and "dressing" does
not apply: there is no bare prediction being patched.  Every factor
that appears in this code (the QED vertex factor (1 - alpha(0) / (2 pi))
in leptons.py, the 30-mode instanton vertex in scale.py, the Singh
ratio 16 in constants.py, the bridge factor sqrt(32/27) in quarks.py,
the heat-kernel coefficient 1/6 - xi_bridge in gravity.py) is one
reading of the E_8(1) ⊃ G_2(1) × F_4(1) algebra at a particular
emergence layer.  Different layers expose different readings of the
same structure; none is more fundamental than another.  Either a
reading is kept (and the residual is what the algebra delivers) or
the previous reading is reported (with whatever residual that layer
left).

Concrete example -- alpha(0) and alpha_em = pi/512.  These are not
"physical alpha and bare alpha"; they are two readings of one
coupling at two emergence layers.  alpha_em = pi/512 is the reading
at the conformal-embedding scale (Layer 7, Singh ratio + WZW), used
inside the instanton vertex.  alpha(0) = pi^2/(256(2pi-1)) is the
reading at the IR pole-mass scale (Layer 7b), reached by passing
the (7,26) bridge self-interference: a marginal (h_bridge = 1)
primary with unit coupling (D^2_local = 1) in a topological coset
(c_coset = 0).  The same (7,26) bridge -- SAME sector, SAME
conformal weight h_bridge = 1 -- also gives Newton's constant via
182 scalar-like heat-kernel channels with xi_bridge = alpha_G2 *
(1/2) * h_bridge = 1/(48 pi).  EM and induced gravity are two
readings through one bridge, not two sectors.

================================================================================
Layer-to-module map
================================================================================

Layer | Module       | Role
------+--------------+-------------------------------------------------------
  0   | algebra.py   | Conformal embedding: Casimirs, central charges, indices.
  1   | scale.py     | Electroweak scale from instanton suppression of M_Pl.
  2   | leptons.py   | Brannen Z_3 closure -> charged-lepton masses.
  3   | quarks.py    | F_4 triality + WZW emergence -> six quark masses.
  4   | wzw.py       | SU(3)_3 modular tensor category: S-matrix, Verlinde.
  4b  | octonions.py | G_2/SU(3) Clebsch-Gordan -> B/A = sqrt(2).
  5   | ckm.py       | CKM matrix from D^(6) boundary structure.
  6   | pmns.py      | PMNS from the conjugation modular invariant.
  6b  | neutrinos.py | F_4 singlets -> rank-2 seesaw, m_1 = 0.
  7   | alpha_s.py   | Strong coupling from embedding-index chain.
  7b  | alpha_bridge | alpha(0) from one-loop self-interference of (7,26).
  8   | higgs.py     | Higgs mass from F_4(1) fusion + bridge threshold.
  10  | gravity.py   | Sakharov induced gravity from (7, 26) bridge.
  11  | baryogenesis | eta_B from G_2 instanton + Fano orientation.
  --  | scorecard.py | Global scorecard + emergence tree.

The structural algebraic identities that previously lived in
proofs.py (38/9, 217/18, 1165/12, 289/432, 32/27, 649/972,
N_vertex = 30, h_bridge = 1, D^2_local = 1, c_coset = 0) are now
asserted directly in the test suite (tests/test_scorecard.py)
and by inline asserts inside each module's derive() function.
The test suite is the certificate.
"""
