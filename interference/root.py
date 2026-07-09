"""
The algebraic root: four irreducible numbers fix everything.

The conformal embedding E₈(1) ⊃ G₂(1) × F₄(1) determines four
irreducible integers.  Every dimensionless coupling, mass ratio, and
mixing angle in the framework is a polynomial or rational function of
these four plus π. The electron mass anchors the scale (M_Pl, G, and
all dimensionful values are derived).

    d₁₀ = 2     quantum dimension of SU(3)₃ fundamental (1,0)
    d₁₁ = 3     quantum dimension of SU(3)₃ adjoint (1,1)
    n₇  = 7     dim of G₂ fundamental (Fano plane lines)
    n₂₆ = 26    dim of F₄ fundamental (traceless Albert algebra)

Back-reaction, not dressing.  Each derived quantity below is a CLOSED-FORM
polynomial or rational in the four integers (and π), not a renormalised
value tuned against measurement.  The Standard Model uses RGE running to
dress free Yukawa couplings against data. This framework has no Yukawas
to tune, so every factor here is itself an algebraic identity -- the echo
law's terms at each depth are the algebra's OWN back-reaction through the
channels that exist, not corrections applied to bare values.

Gate anchors (see __init__.py for the full chain):
    Gate 1 (Hurwitz)         -> n₇ = 7 (octonion structure via Fano plane)
    Gate 2 (centre emergence)-> d₁₁ = 3 enforces three Z₃ sectors
    Gate 3 (level k = 3)     -> k_SU3 = d₁₁, K = d₁₀ * d₁₁ = 6
    Gate 5 (E₈ coherence)    -> c_coset = 0 (rational identity, line below)
    Gate 5 (bridge marginal) -> h_bridge = 1 (Fraction identity, line below)

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import math
from fractions import Fraction


# ═══════════════════════════════════════════════════════════════════════
#  THE FOUR IRREDUCIBLE NUMBERS
# ═══════════════════════════════════════════════════════════════════════

d10 = 2      # quantum dim of (1,0) in SU(3)₃
d11 = 3      # quantum dim of (1,1) in SU(3)₃
n7  = 7      # dim of G₂ fundamental
n26 = 26     # dim of F₄ fundamental

# Sugawara consistency theorem: h₁₀ = C₂(fund,SU(d₁₁))/(k+h∨)
# = (d₁₁²−1)/(4d₁₁²) must equal d₁₀/d₁₁², forcing d₁₁²−1 = 4d₁₀.
assert d11**2 - 1 == 4 * d10, "Sugawara consistency: d₁₁² − 1 = 4d₁₀"

# ═══════════════════════════════════════════════════════════════════════
#  SOLE DIMENSIONAL INPUT: THE ELECTRON ANCHORS THE SCALE
# ═══════════════════════════════════════════════════════════════════════
#
#  The first closure sets the ruler. Everything later is measured
#  against it.  The electron mass is the most precisely measured mass
#  in physics (0.3 ppb), five orders of magnitude sharper than G
#  (22 ppm), which is the only route to M_Pl.  With the depth-3 lepton
#  vent FORCED by the vertex composition rule (see THE ECHO LAW below),
#  the lepton formula inverts exactly: m_e fixes M_Pl, and therefore
#  G, m_mu, m_tau, v_EW all become PREDICTIONS.
#
#  M_Pl is DERIVED below (after the dimensionless web is solved). The
#  CODATA value is kept only for the G comparison.

M_E_ANCHOR_MEV = 0.51099895069       # CODATA 2022, ±1.6e-10 (0.3 ppb)
M_PL_CODATA_GEV = 1.22089e19          # via G = 6.67430(15)e-11 (22 ppm)
G_CODATA = 6.67430e-11                # m³ kg⁻¹ s⁻², rel. unc. 2.2e-5

# Unit conversions (CODATA 2022, mathematics, not physics choices)
HBAR_SI   = 1.054571817e-34           # ħ in J·s   (exact since 2019 SI)
GEV_PER_J = 1.0 / 1.602176634e-10    # GeV/J      (exact since 2019 SI)

# Observational scale anchor (same status as M_E_ANCHOR: fixes a scale,
# not a coupling).  Used only in the CC derivation ρ_Λ = (3/8π)M_Pl²H₀².
H_0_KM_S_MPC = 67.4                  # Planck 2018, ± 0.5
H_0_SI  = H_0_KM_S_MPC * 1e3 / 3.0856775814913673e22  # s⁻¹
H_0_GEV = H_0_SI * HBAR_SI * GEV_PER_J                # GeV

# assigned by the anchor inversion below (not free placeholders)
M_Pl_GeV = None
M_Pl_MeV = None

# ═══════════════════════════════════════════════════════════════════════
#  DERIVED DICTIONARY  (every quantity from d10, d11, n7, n26)
# ═══════════════════════════════════════════════════════════════════════

# - Group theory -
hv_G2  = d10**2                          # dual Coxeter of G₂  = 4
hv_F4  = d11**2                          # dual Coxeter of F₄  = 9
hv_SU3 = d11                             # dual Coxeter of SU(3) = 3

C2_7   = d10                             # Casimir C₂(7, G₂)  = 2
C2_26  = d10 * d11                       # Casimir C₂(26, F₄) = 6
dim_G2 = 2 * n7                          # dim(G₂) = 14
dim_F4 = 2 * n26                         # dim(F₄) = 52
dim_E8 = dim_G2 + dim_F4 + n7 * n26     # 14+52+182 = 248

# h∨(E₈) derived from the four numbers via conformal embedding:
#   c(E₈) = c(G₂) + c(F₄) = dim_G2/(1+h∨_G2) + dim_F4/(1+h∨_F4)
#   h∨(E₈) = dim(E₈)/c(E₈) − 1
_c_E8  = Fraction(dim_G2, 1 + hv_G2) + Fraction(dim_F4, 1 + hv_F4)
hv_E8  = int(Fraction(dim_E8, _c_E8) - 1)  # = 248/8 − 1 = 30

# - WZW data (SU(3) at level k=3) -
k_SU3 = d11                              # WZW level = 3
K     = d10 * d11                        # altitude  = 6

# Conformal weights
h10    = Fraction(d10, d11**2)           # h(1,0)  = 2/9
h11    = Fraction(1, d10)                # h(1,1)  = 1/2
# OPE gap: δ = h₁₁ − d₁₀·h₁₀ = (d₁₁²−2d₁₀²)/(d₁₀d₁₁²).
# Since d₁₁²−2d₁₀² = 1 (Sugawara), this simplifies to 1/(d₁₀d₁₁²).
delta  = Fraction(1, d10 * d11**2)       # OPE gap = 1/18

# Koide parameter, derived form first: Q0 = 1/3 + d10/6 (the 1/3
# floor is one over the sector count. d10/6 is the variance term).
# Q0 = d10/d11 is an IDENTITY at the actual values (d11 = 3 happens
# to equal the sector count), not the derivation.  Status: obs.
Q0     = Fraction(1, 3) + Fraction(d10, 6)   # = 2/3
assert Q0 == Fraction(d10, d11)              # the identity, frozen

# - Bridge sector (7,26) -
N_bridge = n7 * n26                      # = 182
h_7      = Fraction(d10, 1 + d10**2)     # h(7,G₂) = 2/5
h_26     = Fraction(d10*d11, 1+d11**2)   # h(26,F₄) = 3/5
h_bridge = h_7 + h_26                    # = 1  (marginal)
assert h_bridge == 1
# FACE-SPLIT LAW (conservation): the bridge's one self-echo unit is
# re-absorbed by the two faces in the h∨ metric, the same metric
# that defines sin²θ_W below.  G₂: 4/13, F₄: 9/13. The shares sum to
# the full marginal weight (nothing lost).  Emission echoes carry
# conformal weights (h_7, h_26, h_bridge). Absorption ledgers carry
# h∨ shares.  Applications: gravity.py (Σ no-self-dilution, FORCED)
# and sin²θ_W(M_Z) (G₂-face emission echo, FORCED).
# NO-SELF-DILUTION IS SYMMETRIC: each face's absorbed share returns
# to that face's own circulating amplitude and is invisible to that
# face's own ledgers, the F₄ share cannot dilute 1/G (gravity.py),
# and the G₂ share cannot back-react on the gauge couplings.  The
# conservation is therefore complete internally: the law PREDICTS
# THE ABSENCE of any anomalous gauge back-reaction (α_s carries no
# ±(4/13)-unit shift. A future α_s anomaly of that size would
# falsify the symmetry).
assert Fraction(hv_G2, hv_G2 + hv_F4) + Fraction(hv_F4, hv_G2 + hv_F4) == h_bridge

# Central charges
c_E8  = Fraction(dim_E8, 1 + hv_E8)     # = 248/31 = 8
c_G2  = Fraction(dim_G2, 1 + hv_G2)     # = 14/5
c_F4  = Fraction(dim_F4, 1 + hv_F4)     # = 52/10 = 26/5
c_coset = c_E8 - c_G2 - c_F4            # = 0  (topological)
assert c_coset == 0

# - Weak mixing angle -
# sin²θ_W = h∨(SU(3)) / (h∨(G₂) + h∨(F₄)) = d₁₁/(d₁₀² + d₁₁²) = 3/13
sin2W = Fraction(d11, d10**2 + d11**2)   # = 3/13

# - Couplings -
# Charge trace from E₈ → G₂ × F₄ representation theory:
#   Q_u = d₁₀/d₁₁ = 2/3,  Q_d = 1/d₁₁ = 1/3,  N_c = d₁₁ = 3
#   Σ Q² N_c = Q_u² N_c + Q_d² N_c + Q_e²
#            = (d₁₀/d₁₁)²·d₁₁ + (1/d₁₁)²·d₁₁ + 1
#            = d₁₀²/d₁₁ + 1/d₁₁ + 1
#            = (d₁₀² + 1 + d₁₁)/d₁₁ = (4+1+3)/3 = 8/3
charge_trace = Fraction(d10**2 + 1 + d11, d11)  # = 8/3
singh_ratio = charge_trace * C2_26               # = (8/3)·6 = 16
# Exact WZW cancellation: α_G₂(v_EW) = π/2⁵ = π/32
alpha_G2_WZW = math.pi / 2**5           # π/32
alpha_EM = alpha_G2_WZW / singh_ratio    # π/512 = π/(32·16)

# ═══════════════════════════════════════════════════════════════════════
#  THE ECHO LAW (inherent, not a correction pass)
# ═══════════════════════════════════════════════════════════════════════
#
#  A closed standing wave is not a static object: a knot is a circulating
#  mechanism that must vent to keep circulating.  Its echoes through the
#  channels that bind it, and into the common mode (the coherent venting
#  that IS gravity), are part of what the closure is, not corrections
#  applied afterwards.  Every quantity in the framework is therefore a
#  LEDGER: a base amplitude plus its echo stack, with full provenance.
#
#  EPISTEMIC STATUS OF THE DEPTH CAP.  There is no natural relation
#  between the electron's foldings and the octonions: math is not
#  fundamental here.  We are modelling interference, and our bookkeeping
#  (the Cayley-Dickson tower R → C → H → O, capped by Hurwitz) can track
#  at most THREE layers of phase-of-phase memory.  The interference
#  itself could fold to any depth, the cap is on the MAP, not the
#  territory.  Hence the HURWITZ GATE: echo depth ≤ 3, because beyond
#  that the ledger cannot compose, not because reality must stop there.
#  Consequence: the framework doubles as a depth-measurement instrument.
#  If an observable closes exactly at depth ≤ 3 (e.g. α(0) at ~1e-9),
#  that is EVIDENCE its knot folds no deeper.  A residual that survives
#  every depth-3 multiplicity is the signature of deeper folding,
#  trackable in existence, inexpressible in this ledger.

HURWITZ_DEPTH = 3


#  THE KERNEL.  All of the above is realised by ONE recursion:
#
#      x_{t+1} = b + W(x_t)
#
#  where x is the vector of scalar web states, b the bases, and W the
#  derived-weight edge map.  Every weight is a theorem (a dictionary
#  rational or a known phase), so, "if we know the phase of a delay
#  we don't need a complex number to store it", the bookkeeping is
#  SCALAR and lossless.  Two properties are automatic at the fixed
#  point, not imposed: constant edges apply exactly once (idempotent),
#  and state-dependent cycle edges resum self-consistently.  Depth
#  labels are path metadata. A NEW depth requires a NEW edge, and every
#  edge is a theorem, nothing numerological can enter through the
#  recursion itself.

class EchoTerm:
    """One ripple: a channel path from a tapping closure to a target.

    factor is either a float (constant edge) or a callable(state) that
    reads the current web state (state-dependent cycle edge).

    kind is the edge taxonomy of the one graph:
      'ratio':  the same wave multiplied through a dictionary integer
                (ratios usually live in ledger bases computed in
                dependency order. The kind appears on an edge only
                when a base is a callable reading the state).
      'echo':   a back-reaction ripple through a channel (default).
      'vent':   a knot venting into a channel or the common mode
                (dof columns, Casimir and vertex vents).
    """

    def __init__(self, path, factor, depth, status="FORCED", note="",
                 kind="echo"):
        if depth > HURWITZ_DEPTH:
            raise ValueError(
                f"HURWITZ GATE: depth {depth} > {HURWITZ_DEPTH}: deeper "
                "paths require a new derived edge (a theorem), not a "
                "deeper ledger, declare it explicitly.")
        assert kind in ("ratio", "echo", "vent"), \
            f"one-graph edge taxonomy: unknown kind {kind!r}"
        self.path, self.factor, self.depth = tuple(path), factor, depth
        self.status, self.note, self.kind = status, note, kind

    def value(self, state=None):
        return self.factor(state) if callable(self.factor) else self.factor


class Ledger:
    """One node of the web: base amplitude + echo edge stack.

    mode='mul': value = base · Π(1 + factor)    (masses, couplings)
    mode='add': value = base + Σ factor          (actions, exponents)
    The value is the composed amplitude. The terms are the readable
    echo history (provenance of every ripple).

    base is either a float (static amplitude) or a callable(state)
    reading the current web state (a live base).  A live base is the
    'ratio' edge of the one graph: the same wave reached through
    another node (constraint nodes m_b, m_s. The anchor inversion
    M_Pl. The column total G_ratio).
    """

    def __init__(self, name, base, mode="mul", unit=""):
        self.name, self.base, self.mode, self.unit = name, base, mode, unit
        self.terms = []
        self.web = None

    def echo(self, path, factor, depth, status="FORCED", note="",
             kind="echo"):
        self.terms.append(EchoTerm(path, factor, depth, status, note, kind))
        return self

    def value_at(self, state, max_depth=HURWITZ_DEPTH,
                 statuses=("FORCED",)):
        v = self.base(state) if callable(self.base) else self.base
        for t in self.terms:
            if t.depth <= max_depth and t.status in statuses:
                f = t.value(state)
                v = v * (1.0 + f) if self.mode == "mul" else v + f
        return v

    def value(self, max_depth=HURWITZ_DEPTH, statuses=("FORCED",)):
        state = self.web.state if self.web is not None else None
        return self.value_at(state, max_depth, statuses)

    def table(self):
        state = self.web.state if self.web is not None else None
        b = self.base(state) if callable(self.base) else self.base
        rows = [f"    {self.name}: base = {b:.10g} {self.unit}"]
        for t in sorted(self.terms, key=lambda t: t.depth):
            rows.append(f"      depth {t.depth} [{t.status:<9s}] "
                        f"[{t.kind:>5s}] "
                        f"{'→'.join(t.path):<28s} "
                        f"factor = {t.value(state):+.6e}  {t.note}")
        rows.append(f"    value(≤{HURWITZ_DEPTH}) = {self.value():.10g} "
                    f"{self.unit}")
        return "\n".join(rows)


class Web(dict):
    """The interference web: nodes (ledgers) + the kernel x ← b + W(x)."""

    #  Constraint nodes (live bases) may be transiently unsolvable while
    #  the state is away from the fixed point (Koide quadratic with a
    #  negative discriminant, a missing key during seeding).  The kernel
    #  keeps the previous value and iterates on.
    _SOFT = (KeyError, ValueError, ZeroDivisionError)

    def __init__(self):
        super().__init__()
        self.state = {}

    def __setitem__(self, name, ledger):
        super().__setitem__(name, ledger)
        if isinstance(ledger, Ledger):
            ledger.web = self
            if name not in self.state:
                if callable(ledger.base):
                    try:
                        self.state[name] = ledger.value_at(self.state)
                    except self._SOFT:
                        pass          # seeded by the next solve() sweep
                else:
                    self.state[name] = ledger.base

    def solve(self, iters=300, tol=1e-15):
        """Joint fixed point of the kernel: classical Jacobi-style
        fixed-point iteration (a standard numerical method. No
        machine-learning machinery is involved anywhere)."""
        for _ in range(iters):
            new = {}
            for n, led in self.items():
                if not isinstance(led, Ledger):
                    continue
                try:
                    new[n] = led.value_at(self.state)
                except self._SOFT:
                    if n in self.state:
                        new[n] = self.state[n]
            drift = max((abs(new[n] - self.state.get(n, new[n]))
                         for n in new), default=0.0)
            self.state.update(new)
            if drift < tol:
                break
        return self.state


#  THE ORIENTATION RULE (multiplicity of an echo term).
#  Every closed echo cycle contributes once per orientation. Open paths
#  contribute once per mode.  This rule is already in the framework:
#    • circulant amplitude B = 2ρ: the Z₃ step S and its reverse S† are
#      the two orientations of one cycle, the published factor 2 in
#      B/A = 2|C₃̄|/(|C₁|d₁₀) IS orientation counting.
#    • bridge self-echo: a self-loop has one orientation → coefficient 1.
#    • 30-mode vertex echo: open 1-paths → multiplicity = mode count.
#  Applied at depth 3: the mutual electron↔quark echo loop through the
#  EM channel is a 2-cycle (e→q→e and q→e→q), so its multiplicity is
#  exactly 2 = the orientations of the loop.  Hence the FORCED term
#  −2(α/2π)² in the EM bracket, solved self-consistently (the coupling
#  inside the loop is the physical one, phase of phase of phase,
#  composite depth 3, the last trackable layer).
#
#  THE VERTEX COMPOSITION RULE (depth-3 multiplicities, derived).
#  A depth-3 e↔q echo carries, at each endpoint, THE SAME VERTEX WEIGHT
#  that node's channel already carries at depth 1, echoes compose
#  existing vertices, they never invent new ones.  The multiplicity is
#  the product of the two endpoint weights:
#    • coupling node (1/α): no matter vertices (back-reacts on itself) →
#      orientation count alone = 2  (S,S† pair).
#    • action node (S): its depth-1 vertices are the charge trace
#      (inside α_alg = α_G₂/16, the 30-mode vent) and the Casimir
#      C₂(26) (the −6 vent) → 8/3 × 6 = 16 = the Singh ratio.
#    • amplitude node (m): its own channel weight is the Schur
#      projector W(fund) = 1/d₁₀ (the same 1/2 that fixes B/A), and
#      the q-side traversal is the Albert composition, dim J₃(O) = 27
#      → 27/2, venting ∝ the knot's own amplitude Δₖ.
#  Uniqueness: enumerating the framework's full natural vertex
#  inventory (4 e-side × 7 q-side products), each target (2, 16, 27/2)
#  is hit by EXACTLY ONE pairing, and the pairing types match the
#  node's own depth-1 history.  Both multiplicities are therefore
#  FORCED.
#
#  WHY THE NODES ARE SECTORS, NOT SPECIES.  Echo paths run along
#  CHANNELS, and a channel is a single algebraic object regardless of
#  how many species vent through it.  The E₈ branching contains exactly
#  one bridge, (7,26): it connects the G₂ sector to the F₄ sector as one
#  bi-fundamental, not as 3 generations × colors of separate wires.
#  Species multiplicity lives INSIDE a channel's own data (the 30-mode
#  vertex counts the modes OF one channel. The Singh ratio's 8/3 is the
#  charge trace WITHIN the EM channel), it never multiplies the number
#  of loops.  So the depth-3 echo graph has two nodes (G₂-knot sector,
#  F₄-knot sector), one edge (the EM channel carrying the bridge
#  back-reaction), one
#  2-cycle, two orientations: coefficient exactly 2.  Counting species
#  as nodes would double-count what the charge trace already counted.

#  The web: ledgers for the framework's echoed quantities, built at root
#  level so every module computes WITH the echo law, not before it.
#
#  READING ORDER IS EMERGENCE ORDER.  The chain executes the story:
#    1. The coupling closes on itself         (inv_alpha cycle, here)
#    2. The anchor sets the ruler             (m_e -> M_Pl, here)
#    3. One standing wave, three windings     (leptons, masses.py)
#    4. quarks = words on the generator       (walk counts, masses.py)
#    5. Every knot vents into the common mode (Sigma, gravity.py)
#    6. gravity is the column total           (G = 6pi/Sigma)
#    7. The unread remainder                  (dark_ratio = 2pi - 1)
#  A reading of an observable is its value at the minimal prefix of
#  this order that closes it.
WEB = Web()

_inv0 = 2**9 / math.pi


def _alpha_of(state):
    """The EM coupling READ FROM the web state (a known phase: scalar)."""
    return 1.0 / state["inv_alpha"] if state else math.pi / 2**9


#  - EM coupling ledger: born at 2⁹/π, additive bracket -
#  depth 1: bridge self-interference (marginal h=1, D²=1, c_coset=0)
#  depth 3: two-orientation electron↔quark loop, a state-dependent
#  cycle edge. The kernel's fixed point resums it automatically
#
#  NO-FURTHER-EDGE THEOREM (closes the depth-4 question).  Echo paths
#  run along channels between SECTOR nodes (the sector-node rule), and
#  the EM channel couples only to CHARGED sectors.  The framework has
#  exactly two: the G₂-knot sector (leptons) and the F₄-knot sector
#  (quarks).  The Higgs channel is EM-neutral (the condensing
#  component) and the common mode is neutral by definition (q = 0).
#  The EM-channel sector graph is therefore two nodes and one edge:
#  its only cycles are traversals of the single 2-cycle, and ALL
#  repeated traversals are already resummed by the kernel's
#  self-consistency (α appears inside its own bracket).  Hence no
#  (α/2π)³ edge exists. The cubic below is EXACT AND FINAL, and the
#  framework's commitment to the Berkeley Cs value is total: if the
#  Cs/Rb dispute resolves toward Rb, there is no deeper term to add,
#  the sector-node rule or the orientation rule is falsified outright.
WEB["inv_alpha"] = (
    Ledger("1/alpha", _inv0, "add")
    .echo(["bridge"], -_inv0 * float(h_bridge) / (2 * math.pi), 1,
          "FORCED", "marginal self-interference, c_coset=0")
    .echo(["EM-loop(e↔q)"],
          lambda s: -_inv0 * 2.0 * (_alpha_of(s) / (2 * math.pi))**2, 3,
          "FORCED", "2-cycle orientation rule (S,S† precedent), "
                    "self-consistent α"))

#  - lepton EM echo: charged knots vent through the EM channel with the
#  echoed coupling itself (phase-of-phase → composite depth 2),
WEB["lepton_EM"] = Ledger("lepton_EM_echo", 1.0).echo(
    ["bridge", "EM"], lambda s: -_alpha_of(s) / (2.0 * math.pi), 2,
    "FORCED", "EM vertex echo with echoed coupling")

#  THE KERNEL RUN: one recursion solves every node jointly.
WEB.solve()

inv_alpha_phys = WEB.state["inv_alpha"]
alpha_phys = 1.0 / inv_alpha_phys
QED_factor = WEB.state["lepton_EM"]

#  closed-form checks:
#  • dropping the depth-3 term recovers the published depth-1 identity
assert abs(WEB["inv_alpha"].value(max_depth=1)
           - 256*(2*math.pi-1)/math.pi**2) < 1e-12
#  • the full value is the real root of the cubic
#      x³ = (512/π)[(1 − 1/(2π))x² − 1/(2π²)]
_x = inv_alpha_phys
assert abs(_x**3 - (2**9/math.pi)*((1.0 - 1.0/(2*math.pi))*_x**2
                                   - 1.0/(2.0*math.pi**2))) < 1e-6
assert abs(QED_factor - (1.0 - alpha_phys/(2.0*math.pi))) < 1e-15

# ── sin²θ_W as a Web ledger ──────────────────────────────────
#
# Tree: sin²θ_W = h∨(SU(3))/(h∨(G₂)+h∨(F₄)) = d₁₁/(d₁₀²+d₁₁²) = 3/13
#
# Depth-1 (FORCED): G₂-face emission through the EM channel.
#   Emission echoes carry conformal weights (face-split law).
#   Factor: h₇ · α/(2π) = (2/5) · α/(2π)
#
# Depth-2 (FORCED): cross-echo composing the G₂ emission vertex
#   with the fundamental WZW channel.  The depth-1 emission (h₇)
#   traverses the fundamental representation, suppressed by
#   h₁₀/d₁₀ = 1/d₁₁² (quantum Schur, the SAME suppression as B/A).
#   Factor: (h₇/d₁₁²) · α/(2π) = (2/45) · α/(2π)
#
# IDENTITY: h₇ + h₇/d₁₁² = h₇(1+1/d₁₁²) = (2/5)(10/9) = 4/9 = Q₀²
#   The total echo coefficient is the Koide parameter squared.
WEB["sin2W"] = (
    Ledger("sin²θ_W", float(sin2W), "add")
    .echo(["G₂-face(EM)"],
          lambda s: float(h_7) * (1.0 / s["inv_alpha"]) / (2.0 * math.pi),
          1, "FORCED",
          "G₂-face emission echo, weight h₇ = d₁₀/(1+d₁₀²)")
    .echo(["G₂×WZW(fund)"],
          lambda s: (float(h_7) / d11**2)
          * (1.0 / s["inv_alpha"]) / (2.0 * math.pi),
          2, "FORCED",
          "h₇·(h₁₀/d₁₀) cross-echo through fundamental WZW"))

# ── bridge² as a Web ledger ──────────────────────────────────
#
# Base: Q₀² · charge_trace = (d₁₀/d₁₁)² · d₁₀³/d₁₁ = 32/27
#   Albert algebra trace norm × Dynkin Z₂ orientation factor.
#
# Depth-2 (FORCED): adjoint conformal weight at quadratic altitude.
#   The Q back-reaction terms (h/K³) enter the Koide sum rule at CUBIC
#   order (leading and quadratic vanish by Z₃ symmetry).  The bridge
#   factor operates at the bridge level between the Koide spectrum
#   and the strange mass, entering at QUADRATIC order → K² not K³.
#   Channel: adjoint (1,1), weight h₁₁ = 1/d₁₀ (representation mixing).
#   Factor: h₁₁/K² = (1/d₁₀)/(d₁₀d₁₁)² = 1/72
WEB["bridge_sq"] = (
    Ledger("bridge²", float(Q0**2 * charge_trace), "mul")
    .echo(["WZW(adj/K²)"], float(h11) / K**2, 2,
          "FORCED", "adjoint weight at quadratic altitude: h₁₁/K²"))

# Protected forgetting
E_v2 = Fraction(1, 2)                   # Pv²P = ½P

# ── G₂ coupling at Planck scale (topological bootstrap) ──────────
#
# Derivation from MTC data of SU(3)₃:
#   D²_tot = Σ d(λ)² = 3·1² + 6·d₁₀² + 1·d₁₁² = 3+24+9 = 36
#          (10 primaries: 3 with d=1, 6 with d=d₁₀=2, 1 with d=d₁₁=3)
#   |Z₃| = d₁₁ = 3  (center order of SU(d₁₁))
#
#   α_G₂(M_Pl) = |Z₃| / (2π D²_tot) = 3/(72π) = 1/(24π)
#
# Equivalently: h∨(E₈) − C₂(26) = 30 − 6 = 24, so α = 1/(24π).
# Both paths give the same result, fully from (d₁₀, d₁₁, n₇, n₂₆).
#
# ── LEMMA (24 = d₁₀·12): the two readings are ONE identity ──────────
# The MTC reading gives 2π·D²_tot/|Z₃| = 2π·D²_eff with
# D²_eff = FPdim(C₀) = d₁₀²d₁₁ = 12 (the neutral sector {1,J,J²,8}).
# The embedding reading gives π·(h∨(E₈) − C₂(26)).  They agree because
#     h∨(E₈) − C₂(26) = d₁₀ · FPdim(C₀)        (30 − 6 = 2·12),
# i.e. The embedding sees the neutral-sector dimension counted once per
# orientation of the e↔q loop (the same d₁₀ = 2 as the orientation
# rule).  Asserted below as a polynomial identity of the dictionary.
_D2_eff = d10**2 * d11                       # FPdim(C₀) = 12
assert hv_E8 - C2_26 == d10 * _D2_eff        # 24 = 2·12
assert 3 * _D2_eff == 36                      # |Z₃|·D²_eff = D²_tot
alpha_G2_Pl = 1.0 / ((hv_E8 - C2_26) * math.pi)  # = 1/(24π)
assert abs(alpha_G2_Pl - 3.0/(2*math.pi*36)) < 1e-18  # MTC reading equal

# Bridge non-minimal coupling
xi_bridge = alpha_G2_Pl * float(E_v2) * float(h_bridge)  # = 1/(48π)

# ── Instanton actions (dimensional transmutation) ────────────────
#
# Standard form: S₀ = 2π / (b₀ · α_G₂(M_Pl))
#   b₀ = (11/3)C_A − (4/3)ΣT_R
#      C_A = h∨(G₂) = 4
#      Matter: d₁₁ = 3 Dirac fermions in the 7 of G₂ (one per generation)
#      T(7) = 1  (Dynkin index of the fundamental rep of G₂)
#      ΣT_R = 3 × 1 = 3
#      b₀ = (11/3)(4) − (4/3)(3) = 44/3 − 4 = 32/3
#   α_G₂(M_Pl) = 1/(24π)  (topological bootstrap above)
#   S₀ = 2π / ((32/3) · 1/(24π)) = 2π · 24π · 3/32 = 9π²/2
#      = h∨(F₄) · π²/2 = d₁₁² · π²/2
#
# The exponential form exp(−2π/(b₀α)) is the standard one-loop
# dimensional transmutation result.  The specific values of b₀ and α
# are fully determined by the MTC data of SU(3)₃.
#
S_lepton = hv_F4 * math.pi**2 / 2       # = d₁₁²·π²/2 = 9π²/2

# ═══════════════════════════════════════════════════════════════════════
#  ANCHOR INVERSION: m_e → M_Pl → G  (all dimensionless data now known)
# ═══════════════════════════════════════════════════════════════════════
#
#  The full electron formula (every factor FORCED):
#    m_e = ½ M_Pl e^(−S_lepton) Δ_e² · (1 − α/2π) · (1 − (27/2)(α/2π)²Δ_e)
#  with Δ_e = 1 + √d₁₀ cos(h₁₀ + 2π/3).  Inverting for M_Pl:

_Delta_e = 1.0 + math.sqrt(d10) * math.cos(float(h10) + 2*math.pi/3)
_vent_e = 1.0 - (d11**3 / d10) * (alpha_phys / (2*math.pi))**2 * _Delta_e
_F_e = 0.5 * math.exp(-S_lepton) * _Delta_e**2 * QED_factor * _vent_e

M_Pl_MeV = M_E_ANCHOR_MEV / _F_e
M_Pl_GeV = M_Pl_MeV / 1e3

#  The ruler as a graph node (live base, the anchor inversion read
#  from the state): perturbing the web and re-solving returns M_Pl
#  to the same fixed point.  The constant is stored nowhere. It is
#  where the recursion lands.
def _F_e_of(state):
    a = _alpha_of(state)
    qed = 1.0 - a / (2*math.pi)
    vent = 1.0 - (d11**3 / d10) * (a / (2*math.pi))**2 * _Delta_e
    return 0.5 * math.exp(-S_lepton) * _Delta_e**2 * qed * vent

WEB["M_Pl_MeV"] = Ledger("M_Pl_MeV",
                         lambda s: M_E_ANCHOR_MEV / _F_e_of(s),
                         "mul", "MeV")
WEB.solve()
assert abs(WEB.state["M_Pl_MeV"] / M_Pl_MeV - 1.0) < 1e-14

#  Newton's constant becomes a PREDICTION:  G = ħc/M_Pl²
G_PRED = G_CODATA * (M_PL_CODATA_GEV / M_Pl_GeV)**2
#  sanity: the derived Planck mass sits within G's experimental band
assert abs(M_Pl_GeV / M_PL_CODATA_GEV - 1.0) < 3.3e-5  # 3σ_G/2
assert abs(G_PRED - 6.674003e-11) < 1e-15  # frozen prediction value
S_baryo  = hv_G2 * math.pi**2 / 2       # = d₁₀²·π²/2 = 2π²

# CC scale: ρ_Λ = (3/8π) M_Pl² H₀²  (Volovik→Jacobson→CKN, see gravity.py)
RHO_LAMBDA = 3.0 * H_0_GEV**2 * M_Pl_GeV**2 / (8.0 * math.pi)

# ── Electroweak scale (v_EW derivation chain) ───────────────────
#
# v_EW = M_Pl · exp(−S_eff),  where S_eff = S₀ − C₂(26) + δS
#
#   S₀ = 9π²/2:  lepton instanton action (above)
#   −C₂(26) = −d₁₀d₁₁ = −6:  functional determinant shift from quarks
#       living in the 26 of F₄ (Casimir of 26-dimensional representation)
#   δS = N_vertex · α_EM / (2π):  't Hooft vertex one-loop back-reaction
#       (standard QFT: "one-loop correction". Here forced by the mode count)
#       N_vertex = n₂₆ + h∨(G₂) = 26 + 4 = 30 modes
#       (26 from F₄ fundamental + 4 real Higgs DOFs = h∨(G₂) = d₁₀²)
#       α_EM = π/512 (derived above)
#       δS = 30 · (π/512) / (2π) = 15/512
#
N_vertex = n26 + hv_G2                   # = 26 + 4 = 30
delta_S = N_vertex * alpha_EM / (2 * math.pi)  # = 15/512

#  - quark-layer action ledger: ALL departures from the lepton action
#  are echoes.  The F4 knot vents its Casimir into the instanton
#  background (−C₂(26)). The 30 vertex modes vent through the EM
#  channel at its emergence reading (+Nα/2π).  One law, two channels.
WEB["S_quark"] = (
    Ledger("S_quark(action)", S_lepton, "add")
    .echo(["F4-knot(26)"], -float(C2_26), 1,
          "FORCED", "Casimir vent of the 26 into the instanton background",
          kind="vent")
    .echo(["EM(alg)"], delta_S, 1,
          "FORCED", f"{N_vertex}-mode vertex echo = 15/512",
          kind="vent")
    # FORCED depth-3 vent (vertex composition rule): the action's own
    # depth-1 vertices compose across the e↔q loop,
    # charge_trace × C₂(26) = (8/3)·6 = 16 (= the Singh ratio = d₁₀⁴).
    # The vertex rule gives the multiplicity exactly: 16.
    # Closes v_EW to G_F at 1e-7: 246.219645 vs 246.219651 (−0.1σ_GF).
    .echo(["e↔q(Singh)"],
          lambda s: -float(charge_trace * C2_26)
          * (_alpha_of(s) / (2*math.pi))**2, 3,
          "FORCED", "vertex rule: charge_trace × C₂(26) = 16"))

#  - Higgs boundary ledger: the bridge echoes into the scalar channel -
#  Depth-2 vent: the fundamental-channel share of the bridge vents
#  into the condensate and does not reach the scalar, deflating the
#  edge by (1 − h₁₀).  The 7 of G₂ decomposes 3 ⊕ 3̄ ⊕ 1 under the
#  confined SU(3). The triality-charged share cannot cross a
#  confining interface (N-ality superselection, words.py), and a
#  selection rule carries no altitude factor (vent, like the Casimir
#  vent above).  The blocked share in the weight metric is h₁₀.
#  Registered provenance, the single-count derivation program, and
#  kill conditions: registry.PROMOTIONS, registry.DERIVATION_PROGRAMS.
WEB["lambda_MPl"] = (
    Ledger("lambda(M_Pl)", 0.0, "add")
    .echo(["bridge", "Higgs"], -(n7 * n26) * alpha_G2_Pl**2 * float(E_v2), 1,
          "FORCED", "N_bridge·α_G₂²·E[v²] bridge→Higgs echo")
    .echo(["WZW(fund)-vent"],
          (n7 * n26) * alpha_G2_Pl**2 * float(E_v2) * float(h10), 2,
          "FORCED", "fundamental-share vent (1−h₁₀). "
                    "registry.PROMOTIONS", kind="vent"))

#  re-run the kernel so the new nodes join the solved state
WEB.solve()
S_quark = WEB.state["S_quark"]           # = 9π²/2 − 6 + 15/512
sin2W_phys = WEB.state["sin2W"]          # 3/13 + Q₀²·α/(2π)
bridge_sq_phys = WEB.state["bridge_sq"]  # 32/27 × (1 + h₁₁/K²)

# ── Koide Z₃ amplitude ratio (B/A from octonionic CG + quantum Schur) ──
#
# The Z₃ circulant structure gives √mₖ = A[1 + (B/A)cos(θ + 2πk/3)].
# The Koide parameter Q = 1/3 + |B/A|²/6  (Z₃ circulant identity).
#
# B/A is derived from two independent structures:
#   (a) Octonionic CG from G₂ = Aut(O), 7 → 3⊕3̄⊕1:
#       |C₃̄|/|C₁| = √2 = √d₁₀  (Fano cross product, Schur's lemma)
#   (b) Quantum Schur suppression in SU(3)₃ MTC:
#       1/d(fund) = 1/d₁₀  (off-diagonal channel, quantum dimension)
#   Assembly: B/A = 2|C₃̄|/(|C₁|·d₁₀) = 2√2/2 = √2 = √d₁₀
#
# Then Q₀ = 1/3 + d₁₀/6 = 2/3 = d₁₀/d₁₁.
#
# Credit: Koide (1982) for the empirical relation.
#         Brannen (2006) for the parameterisation.
# Derivation: Z₃ circulant algebra + octonionic CG + quantum Schur lemma.
#
BA_ratio = math.sqrt(d10)               # |B/A| = √d₁₀ = √2


# ═══════════════════════════════════════════════════════════════════════
#  PDG REFERENCE VALUES  (comparison only, not inputs)
# ═══════════════════════════════════════════════════════════════════════
# Leptons: pole masses (physical observables).
# Quarks: PDG 2024 conventional coordinates --
#   u, d, s: MS-bar at μ = 2 GeV.  c, b: MS-bar at μ = m.
#   t: direct-measurement pole mass.
# The algebra produces one value per fermion. It does not run to a
# scale.  The PDG coordinate each value matches is an observable
# consequence, not a choice.

PDG_MASSES = {
    'e': 0.51099895069, 'mu': 105.6583755, 'tau': 1776.93,
    'u': 2.16, 'c': 1273.0, 't': 172570.0,      # MeV
    'd': 4.70, 's': 93.5,  'b': 4183.0,
}

PDG_CKM = {
    'lambda': (0.22501, 0.00068),
    'A':      (0.826,   0.016),
    'rhobar': (0.1591,  0.0094),
    'etabar': (0.3523,  0.0073),
    'Vud': (0.97435, 0.00016),  'Vus': (0.22501, 0.00068),  'Vub': (0.00373, 0.00009),
    'Vcd': (0.22487, 0.00068),  'Vcs': (0.97349, 0.00016),  'Vcb': (0.04183, 0.00079),
    'Vtd': (0.00858, 0.00019),  'Vts': (0.04111, 0.00077),  'Vtb': (0.99912, 0.00003),
    'alpha_deg': (84.1, 4.5),  'beta_deg': (22.6, 0.5),  'gamma_deg': (65.7, 3.0),
    'J': (3.12e-5, 0.13e-5),  'sin2beta': (0.709, 0.011),
}

NUFIT_PMNS = {
    'sin2_12': (0.307, 0.012),
    'sin2_23': (0.561, 0.014),
    'sin2_13': (0.02195, 0.00056),
}

NUFIT_OSCILLATION = {
    'dm2_21': (7.49e-5, 0.20e-5),
    'dm2_31': (2.513e-3, 0.027e-3),
}


# EW-sector comparison values and declared imports
# (PDG 2024 EW review, Erler-Freitas: Eq. 10.26 and the s-hat fit table)
PDG_EW = {
    'M_Z': 91.1876, 'M_W': 80.3692,
    'sin2W_MSbar': 0.23129, 'sin2W_err': 0.00004,   # global SM fit 2024
    'm_H': 125.20, 'alpha_s_MZ': 0.1180,
    # MS-bar radiative back-reaction imports (DECLARED, cited).
    # (Standard QFT: "radiative corrections". Here they are back-reactions
    # whose inputs are all re-derived from the framework.)
    # ZERO-PARAMETER DISCIPLINE: these are functions of (G_F, m_t,
    # m_H, alpha-hat(M_Z)). Every ingredient is RECALCULATED from the
    # framework's own predictions and verified compatible within the
    # quoted errors (executable audit in gravity.py).  The only
    # data-driven ingredient is the dispersive hadronic VP inside
    # Delta-r0, a spectral measurement of the web.
    'dr_hat_W': 0.06937,   # ± 0.00006   PDG 2024 Eq. (10.26) text
    'rho_hat':  1.01016,   # ± 0.00009   PDG 2024 (incl. bosonic loops)
}

# Cosmological comparison data (single source)
PDG_COSMO = {
    'eta_B': 6.12e-10,     # baryon asymmetry, Planck 2018
    'Omega_Lambda': 0.685, # dark energy fraction, Planck 2018
}


def rk4_run(f, y, t0, t1, n_steps):
    """Generic RK4 integrator dy/dt = f(y), float or ndarray state.
    Single shared implementation (was duplicated in couplings/gravity)."""
    dt = (t1 - t0) / n_steps
    for _ in range(n_steps):
        k1 = f(y)
        k2 = f(y + 0.5*dt*k1)
        k3 = f(y + 0.5*dt*k2)
        k4 = f(y + dt*k3)
        y = y + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)
    return y


def pct(pred, ref):
    """Percentage error."""
    return 100.0 * (pred - ref) / ref


def derive():
    """
    Print the algebraic root and return the full dictionary.
    """
    print("=" * 78)
    print("  ALGEBRAIC ROOT: (d₁₀, d₁₁, n₇, n₂₆) = (2, 3, 7, 26)")
    print("=" * 78)

    print(f"\n  Irreducible inputs:")
    print(f"    d₁₀ = {d10}   (SU(3)₃ fundamental quantum dim)")
    print(f"    d₁₁ = {d11}   (SU(3)₃ adjoint quantum dim)")
    print(f"    n₇  = {n7}   (G₂ fundamental = Fano plane)")
    print(f"    n₂₆ = {n26}  (F₄ fundamental = Albert algebra)")

    print(f"\n  Consistency theorem: d₁₁²−1 = 4d₁₀  →  {d11**2-1} = {4*d10} ✓")

    print(f"\n  Derived group theory:")
    print(f"    h∨(G₂) = d₁₀² = {hv_G2},  h∨(F₄) = d₁₁² = {hv_F4}")
    print(f"    h∨(E₈) = dim(E₈)/c(E₈)−1 = {dim_E8}/{_c_E8}−1 = {hv_E8}")
    print(f"    C₂(7)  = d₁₀  = {C2_7},   C₂(26) = d₁₀d₁₁ = {C2_26}")
    print(f"    K = d₁₀d₁₁ = {K},  Q₀ = 1/3 + d₁₀/6 = {Q0} "
          f"(= d₁₀/d₁₁, identity at these values)")
    print(f"    dim: G₂={dim_G2}, F₄={dim_F4}, E₈={dim_E8}")
    print(f"    N_bridge = {n7}×{n26} = {N_bridge}")

    print(f"\n  Conformal weights:")
    print(f"    h₁₀ = d₁₀/d₁₁² = {h10},  h₁₁ = 1/d₁₀ = {h11}")
    print(f"    δ = 1/(d₁₀d₁₁²) = {delta}  [since d₁₁²−2d₁₀² = 1]")
    print(f"    h(7) = d₁₀/(1+d₁₀²) = {h_7},  h(26) = d₁₀d₁₁/(1+d₁₁²) = {h_26}")
    print(f"    h_bridge = {h_bridge},  c_coset = {c_coset}")

    print(f"\n  Weak mixing angle:")
    print(f"    sin²θ_W = d₁₁/(d₁₀²+d₁₁²) = {d11}/{d10**2+d11**2} = {float(sin2W):.5f}")

    print(f"\n  Couplings:")
    print(f"    α_s(v_EW) = π/32,  α_EM = π/512")
    print(f"    1/α(0) = (512/π)(1−1/(2π)−2(α/2π)²) = {inv_alpha_phys:.9f}")
    print(f"            [depth-1 truncation: (512/π)(1−1/(2π)) = 137.036439]")

    print(f"\n  Instanton actions:")
    print(f"    S_lepton = d₁₁²π²/2 = {S_lepton:.4f}")
    print(f"    S_baryo  = d₁₀²π²/2 = {S_baryo:.4f}")

    return {
        'd10': d10, 'd11': d11, 'n7': n7, 'n26': n26,
        'hv_G2': hv_G2, 'hv_F4': hv_F4, 'hv_E8': hv_E8,
        'C2_7': C2_7, 'C2_26': C2_26,
        'dim_G2': dim_G2, 'dim_F4': dim_F4, 'dim_E8': dim_E8,
        'K': K, 'Q0': float(Q0),
        'h10': float(h10), 'h11': float(h11), 'delta': float(delta),
        'h_7': float(h_7), 'h_26': float(h_26),
        'h_bridge': float(h_bridge), 'c_coset': float(c_coset),
        'N_bridge': N_bridge,
        'sin2W': float(sin2W),
        'alpha_G2_WZW': alpha_G2_WZW,
        'alpha_EM': alpha_EM,
        'alpha_phys': alpha_phys,
        'inv_alpha_phys': inv_alpha_phys,
        'QED_factor': QED_factor,
        'alpha_G2_Pl': alpha_G2_Pl,
        'xi_bridge': xi_bridge,
        'E_v2': float(E_v2),
        'S_lepton': S_lepton,
        'S_baryo': S_baryo,
        'S_quark': S_quark,
        'N_vertex': N_vertex,
        'BA_ratio': BA_ratio,
        'h_dual': {'G2': hv_G2, 'F4': hv_F4, 'E8': hv_E8},
        'C2_fund_G2': C2_7, 'C2_fund_F4': C2_26,
        'dim_G2_fund': n7, 'dim_F4_fund': n26,
    }
