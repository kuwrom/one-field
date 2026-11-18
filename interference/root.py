"""
The algebraic root: depth 3 searchlight vocabulary.

The conformal embedding E₈(1) ⊃ G₂(1) × F₄(1) determines four
integers at the depth 3 tracking layer.  Every dimensionless coupling,
mass ratio, and mixing angle in the framework is a closed form
expression in these four plus π.  The electron mass sets the scale.
M_Pl, G, and all other dimensionful values are derived.

    d₁₀ = 2     quantum dimension of SU(3)₃ fundamental (1,0)
    d₁₁ = 3     quantum dimension of SU(3)₃ adjoint (1,1)
    n₇  = 7     dim of G₂ fundamental (Fano plane lines)
    n₂₆ = 26    dim of F₄ fundamental (traceless Albert algebra)

Each quantity has one reading, and its backreaction is part of that
reading.  Each derived quantity below is a CLOSED FORM expression in the
four integers (and π), not a value tuned against measurement.
``Ledger.base`` is the zero edge term of the recursion, not a bare
physical state.  Echo terms are the same quantity's backreaction
through the web's own channels.  Standard RGE coordinates appear only
when a derived reading is transported for comparison.

Gate anchors (see __init__.py for the full chain):
    Gate 1 (normed ledger criterion + Hurwitz classification)
                             -> n₇ = 7 (octonion structure via Fano plane)
    Gate 2 (centre emergence)-> d₁₁ = 3 enforces three Z₃ sectors
    Gate 3 (level k = 3)     -> k_SU3 = d₁₁, K = d₁₀ * d₁₁ = 6
    Gate 5 (E₈ coherence)    -> c_coset = 0 (rational identity, line below)
    Gate 5 (bridge weight)   -> h_bridge = 1 (Fraction identity, line below)

Reference:
    Kahsay, Kibrom Kidane (2026). The Octavian Coherence Gate: The Four
    Irreducible Integers of the E₈(1) ⊃ G₂(1) × F₄(1) Conformal Embedding.
    Zenodo. https://doi.org/10.5281/zenodo.20493955
"""

import cmath
import math
from fractions import Fraction
from . import wiring as _wiring
from . import boundary as _boundary


# ═══════════════════════════════════════════════════════════════════════
#  DEPTH-3 SEARCHLIGHT VOCABULARY (read from wiring.py)
# ═══════════════════════════════════════════════════════════════════════

WIRING_DATA = _wiring.canonical()
d10 = round(WIRING_DATA["d10"])   # quantum dim of (1,0) in SU(3)₃
d11 = round(WIRING_DATA["d11"])   # quantum dim of (1,1) in SU(3)₃
n7  = WIRING_DATA["n7"]            # dim of G₂ fundamental
n26 = WIRING_DATA["n26"]           # dim of F₄ fundamental

# Sugawara CONSISTENCY CHECK: h₁₀ =
# C₂(fund,SU(d₁₁))/(k+h∨) = (d₁₁²−1)/(4d₁₁²) must equal d₁₀/d₁₁²,
# i.e. d₁₁²−1 = 4d₁₀.  The closing values satisfy it, but it does not
# single them out: (6,5), (12,7), ... satisfy it too.  What selects
# (2,3) is the embedding gate together with the level-3 gate.  The
# embedding gate (embedding_uniqueness.py: E₈(1) ⊃ G₂(1)×F₄(1)) fixes
# n₇ = 7, n₂₆ = 26.  The level-3 gate (h(J) = 1 forces k = 3) gives the
# MTC quantum dimensions d₁₀ = 2, d₁₁ = 3, computed in
# wiring.su3_qdims.  The assert below checks this closure's consistency.
assert d11**2 - 1 == 4 * d10, "Sugawara consistency: d₁₁² − 1 = 4d₁₀"

# ═══════════════════════════════════════════════════════════════════════
#  SOLE DIMENSIONAL INPUT: THE ELECTRON SETS THE SCALE
# ═══════════════════════════════════════════════════════════════════════
#
#  The first closure sets the unit.  Everything later is measured in
#  that unit.  The electron mass is the most precisely measured mass in
#  physics (0.3 ppb), five orders of magnitude sharper than G (22 ppm),
#  and G is the only direct route to M_Pl.  With the depth 4 lepton
#  response (see THE ECHO LAW below), the lepton formula inverts
#  exactly: m_e fixes M_Pl, and m_mu, m_tau, v_EW become PREDICTIONS.
#  For G the statement is narrower.  The PREDICTED quantity is M_Pl
#  (+22 ppm above CODATA's Planck mass).  The printed
#  G_PRED = common_mode·ħc/M_Pl² reads M_Pl through the common mode
#  (THE COMMON MODE, below).  It is not a second, independent gravity
#  result.
#
#  THE PRECISION BLOCK AND WHAT ENDS IT.  G = 6.673936e-11 (−2.43σ of
#  CODATA, +0.17σ with a spread added to each input:
#  tests/probes/common_mode_count.py C6).  m_mu = 105.6583758 MeV
#  (ppb, +0.12σ).  m_tau = 1776.9092814 MeV (−0.23σ).
#  v_EW = 246.219645 GeV against G_F (+0.08σ at 1e-7).
#  M_Pl = 1.2209171e19 GeV (output).  An improved G that confirms the
#  CODATA central value at ≤ 1e-5 kills the chain.  So does any of
#  these precision values moving outside its stated pull (CODATA/PDG
#  updates and Belle II for m_tau decide).
#
#  M_Pl is DERIVED below (after the dimensionless web is solved). The
#  CODATA value is kept only for the G comparison.

M_E_ANCHOR_MEV = 0.51099895069       # CODATA 2022, ±1.6e-10 (0.3 ppb)
M_PL_CODATA_GEV = 1.22089e19          # via G = 6.67430(15)e-11 (22 ppm)
G_CODATA = 6.67430e-11                # m³ kg⁻¹ s⁻², rel. unc. 2.2e-5
G_CODATA_SIGMA = 0.00015e-11          # CODATA 2022 standard uncertainty
G_REL_SIGMA = G_CODATA_SIGMA / G_CODATA

# Unit conversions (CODATA 2022, mathematics, not physics choices)
HBAR_SI   = 1.054571817e-34           # ħ in J·s   (exact since 2019 SI)
GEV_PER_J = 1.0 / 1.602176634e-10    # GeV/J      (exact since 2019 SI)
C_SI      = 299792458.0               # c in m/s   (exact)

# The comparison constant is CODATA's Planck mass as printed, held to
# its last digit against sqrt(hbar c^5/G).
assert abs(M_PL_CODATA_GEV / ((HBAR_SI * C_SI**5 / G_CODATA)**0.5 * GEV_PER_J)
           - 1.0) < 4.1e-7

# Comparison scale for the Friedmann critical density
# ρ_crit = (3/8π)M_Pl²H₀², bookkeeping that no prediction reads.
H_0_KM_S_MPC = 67.4                  # Planck 2018, ± 0.5
H_0_SI  = H_0_KM_S_MPC * 1e3 / 3.0856775814913673e22  # s⁻¹
H_0_GEV = H_0_SI * HBAR_SI * GEV_PER_J                # GeV

# assigned by the anchor inversion below
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
# The normalized E8 Cartan current trace is computed directly from all
# 240 roots by wiring.py.  One representative per conjugate +/- pair gives
# sum_(alpha>0)(alpha.q)^2 = 30 for every unit current direction q.
vertex_charge_trace = round(WIRING_DATA["vertex_charge_trace"])
assert vertex_charge_trace == hv_E8

# - WZW data (SU(3) at level k=3) -
k_SU3 = d11                              # WZW level = 3
K     = d10 * d11                        # altitude  = 6

# Conformal weights
h10    = Fraction(d10, d11**2)           # h(1,0)  = 2/9
h11    = Fraction(1, d10)                # h(1,1)  = 1/2
# OPE gap: δ = h₁₁ − d₁₀·h₁₀ = (d₁₁²−d₁₀³)/(d₁₀d₁₁²).
# Since d₁₁²−d₁₀³ = 1 at (2, 3), this simplifies to 1/(d₁₀d₁₁²).
assert d11**2 - d10**3 == 1
delta  = Fraction(1, d10 * d11**2)       # OPE gap = 1/18

# Koide parameter, downstream of the octonionic/MTC channel map:
# |B/A|^2 = d10, and the circulant identity gives
# Q0 = 1/3 + |B/A|^2/6 = 1/3 + d10/6 (the 1/3 floor is one over
# the sector count, d10/6 is the variance term). At d10 = 2 and d11 = 3,
# Q0 also equals d10/d11.
Q0     = Fraction(1, 3) + Fraction(d10, 6)   # = 2/3
assert Q0 == Fraction(d10, d11)              # the identity, frozen

# - Bridge sector (7,26) -
# TWO LEVELS COEXIST IN ONE WEB.  This follows from the embedding
# E₈(1) ⊃ G₂(1) × F₄(1) with SU(3) at level 3 inside G₂.  h_7 and
# h_26 below are LEVEL-1 weights (denominators 1 + h∨).  h10 and h11
# above are LEVEL-3 SU(3) weights (k + h∨ = 6).  The α(0) chain and
# h_bridge use the level-1 weights.  The Koide/mass chain uses the
# level-3 weights.  Any ledger that combines the two (e.g. the sin²θ_W
# echo coefficient h_7·(1 + 1/d₁₁²) below) mixes levels, which the
# embedding permits, and says so where it does.
N_bridge = n7 * n26                      # = 182
h_7      = Fraction(d10, 1 + d10**2)     # h(7,G₂) = 2/5   (level 1)
h_26     = Fraction(d10*d11, 1+d11**2)   # h(26,F₄) = 3/5  (level 1)
h_bridge = h_7 + h_26                    # = 1  (dimension one)
assert h_bridge == 1
# APPORTIONMENT LAW (conservation).  The two faces reabsorb the
# bridge's one self-echo unit in the h∨ metric, the same metric that
# defines sin²θ_W below: G₂ takes 4/13 and F₄ takes 9/13.  The shares
# sum to the full unit weight, so nothing is lost.  Emission echoes
# carry conformal weights (h_7, h_26, h_bridge).  Absorption ledgers
# carry h∨ shares.  Applications: gravity.py (Σ nondilution, FORCED)
# and sin²θ_W(M_Z) (G₂ face emission echo, FORCED).
# NONDILUTION IS SYMMETRIC.  Each face's absorbed share returns to that
# face's own circulating amplitude, and that face's own ledgers do not
# see it.  So the F₄ share cannot dilute 1/G (gravity.py), and the G₂
# share cannot feed back into the gauge couplings.  The conservation is
# therefore complete internally.  The law PREDICTS THE ABSENCE of any
# anomalous gauge backreaction: α_s carries no ±(4/13)-unit shift, and
# a future α_s anomaly of that size would falsify the symmetry.
assert Fraction(hv_G2, hv_G2 + hv_F4) + Fraction(hv_F4, hv_G2 + hv_F4) == h_bridge

# Central charges
c_E8  = Fraction(dim_E8, 1 + hv_E8)     # = 248/31 = 8
c_G2  = Fraction(dim_G2, 1 + hv_G2)     # = 14/5
c_F4  = Fraction(dim_F4, 1 + hv_F4)     # = 52/10 = 26/5
c_coset = c_E8 - c_G2 - c_F4            # = 0  (no coset Virasoro sector)
assert c_coset == 0

# - Weak mixing angle -
# sin²θ_W = h∨(SU(3)) / (h∨(G₂) + h∨(F₄)) = d₁₁/(d₁₀² + d₁₁²) = 3/13
sin2W = Fraction(d11, d10**2 + d11**2)   # = 3/13
# IDENTIFICATION: dual Coxeter numbers (a Casimir quantity) are read
# as a gauge coupling ratio (a Lagrangian quantity).  The h∨ metric is
# stated, not computed.  It is CROSS BRACED, not fixed per observable.
# Two observables use the same h∨ metric: sin²θ_W here (−0.11σ) and the
# prediction that α_s carries no 4/13 shift (couplings.py).  Changing
# the metric to help one would change the other.

# - Couplings -
# Charge trace, an IDENTIFICATION of the four integers with the SM
# charges and the colour count:
#   Q_u = d₁₀/d₁₁ = 2/3,  Q_d = 1/d₁₁ = 1/3,  N_c = d₁₁ = 3
#   Σ Q² N_c = Q_u² N_c + Q_d² N_c + Q_e²
#            = (d₁₀/d₁₁)²·d₁₁ + (1/d₁₁)²·d₁₁ + 1
#            = d₁₀²/d₁₁ + 1/d₁₁ + 1
#            = (d₁₀² + 1 + d₁₁)/d₁₁ = (4+1+3)/3 = 8/3
charge_trace = Fraction(d10**2 + 1 + d11, d11)  # = 8/3
# The dictionary as forms in the quantum dimensions d(1,0) = d10 and
# d(1,1) = d11, so a reader can see which primary each charge reads.
CHARGE_FORMS = {"Q_u": lambda d_f, d_adj: d_f / d_adj,
                "Q_d": lambda d_f, d_adj: 1 / d_adj,
                "N_c": lambda d_f, d_adj: d_adj}
Q_u, Q_d, N_c = (CHARGE_FORMS[k](Fraction(d10), Fraction(d11)) for k in ("Q_u", "Q_d", "N_c"))
assert charge_trace == Q_u**2 * N_c + Q_d**2 * N_c + 1
singh_ratio = charge_trace * C2_26               # = (8/3)·6 = 16
# Exact WZW cancellation: α_G₂(μ*) = π/2⁵ = π/32 at the
# gauge matching scale μ* = M_Pl·e^{−(S₀−C₂(26))} ≈ 253.5 GeV, NOT
# at v_EW = 246.2 GeV (couplings.py derives and holds the scale).
# π/32 is the one loop endpoint of the G₂ lever arm, asserted below.
alpha_G2_WZW = math.pi / 2**5           # π/32, endpoint of the G₂ lever arm
alpha_EM = alpha_G2_WZW / singh_ratio    # π/512 = π/(32·16)
# IDENTIFICATION, narrow scope.  Dividing an embedded U(1) coupling by
# a charge trace over the multiplet is the STANDARD embedding
# normalisation (the same operation as the GUT hypercharge factor).
# What is identified here is WHICH trace applies,
# charge_trace × C₂(26) = (8/3)·6 = 16, not that a trace division
# happens.  All downstream EM precision results (α(0), v_EW→G_F,
# sin²θ_W, M_W, M_Z) inherit this trace.

# ═══════════════════════════════════════════════════════════════════════
#  THE ECHO LAW
# ═══════════════════════════════════════════════════════════════════════
#
#  A closed standing wave is not static.  A knot circulates, and
#  whatever does not cancel locally vents.  Its echoes run through the
#  channels that bind it and into the common mode (the coherent venting
#  that IS gravity).  These echoes are part of the closure itself.
#  Every quantity in the framework is
#  therefore a LEDGER: a base amplitude plus its stack of echoes, each
#  with its source recorded.
#
#  The octonionic searchlight resolves the closed a→b→c→a loop at
#  depth 3. The scalar recursive web can then refine it at depth 4 and
#  beyond through registered derived edges.

SEARCHLIGHT_DEPTH = 3


#  THE KERNEL.  All of the above is realised by ONE recursion:
#
#      x_{t+1} = b + W(x_t)
#
#  where x is the vector of scalar web states, b the bases, and W the
#  edge map, each edge carrying its status.  Derived weights are
#  dictionary rationals or known phases.  The weights that are
#  identified rather than derived are the mixing coordinates.  Because
#  "if we know the phase of a delay we don't need a complex number to
#  store it", the bookkeeping is SCALAR and loses nothing.  Two
#  properties hold automatically at the fixed point: constant edges
#  apply exactly once (idempotent), and state dependent cycle edges
#  resum self-consistently.  Depth labels record each path and its
#  ancestry.  A NEW depth requires a NEW edge, and every post
#  searchlight edge must name the depth-(n-1) ledger it reads.  Every
#  edge is a theorem, so nothing numerological can enter through the
#  recursion itself.

class EchoTerm:
    """One echo: a channel path from the closure it starts at to a target.

    factor is either a float (constant edge) or a callable(state) that
    reads the current web state (state dependent cycle edge).

    kind is the edge taxonomy of the one graph:
      'ratio':  the same wave multiplied through a dictionary integer
                (ratios usually live in ledger bases computed in
                dependency order. The kind appears on an edge only
                when a base is a callable reading the state).
      'echo':   a backreaction through a channel (default).
      'vent':   a knot venting into a channel or the common mode
                (dof columns, Casimir and vertex vents).

    WHY THE PROTON IS STABLE.  The conserved quantity is the
    fermion line terminus.  Every edge in this taxonomy composes
    amplitudes AT a species node, and no edge kind changes a knot's
    species.  So baryon and lepton number are conserved node labels of
    the graph.  A B-violating channel would need an identity changing
    edge kind.  None exists, and none can be added without a new
    theorem.  The assertion below enforces this.  Lowest
    admissible B-violating edge: NONE, at any web depth.  Against
    tau_p > 1.6e34 yr (Super-K) the framework predicts no decay at
    any sensitivity, and any confirmed proton decay kills it.
    """

    def __init__(self, path, factor, depth, status="FORCED", note="",
                 kind="echo", parent=None):
        if not isinstance(depth, int) or depth < 0:
            raise ValueError("web depth must be a nonnegative integer")
        assert kind in ("ratio", "echo", "vent"), \
            f"one graph edge taxonomy: unknown kind {kind!r}"
        self.path, self.factor, self.depth = tuple(path), factor, depth
        self.status, self.note, self.kind = status, note, kind
        self.parent = parent

    def value(self, state=None):
        return self.factor(state) if callable(self.factor) else self.factor


class Ledger:
    """One node of the web: base amplitude + echo edge stack.

    mode='mul': value = base · Π(1 + factor)    (masses, couplings)
    mode='add': value = base + Σ factor          (actions, exponents)
    The value is the composed amplitude.  The terms are the readable
    echo history (the source of every echo).

    base is either a float (static amplitude) or a callable(state)
    reading the current web state (a live base).  A live base is the
    'ratio' edge of the one graph: the same wave reached through
    another node (constraint nodes m_b, m_s. The anchor inversion
    M_Pl. The column total G_ratio).

    web_depth records the recursive interference depth of this reading.
    """

    def __init__(self, name, base, mode="mul", unit="", web_depth=0):
        if not isinstance(web_depth, int) or web_depth < 0:
            raise ValueError("web depth must be a nonnegative integer")
        self.name, self.base, self.mode, self.unit = name, base, mode, unit
        self.web_depth = web_depth
        self.terms = []
        self.web = None

    def echo(self, path, factor, depth, status="FORCED", note="",
             kind="echo", parent=None):
        self.terms.append(EchoTerm(path, factor, depth, status, note, kind,
                                   parent))
        return self

    def value_at(self, state, max_depth=None,
                 statuses=("FORCED",)):
        v = self.base(state) if callable(self.base) else self.base
        for t in self.terms:
            if (max_depth is None or t.depth <= max_depth) \
                    and t.status in statuses:
                f = t.value(state)
                v = v * (1.0 + f) if self.mode == "mul" else v + f
        return v

    def value(self, max_depth=None, statuses=("FORCED",)):
        state = self.web.state if self.web is not None else None
        return self.value_at(state, max_depth, statuses)

    def table(self):
        state = self.web.state if self.web is not None else None
        b = self.base(state) if callable(self.base) else self.base
        rows = [f"    {self.name}: base = {b:.10g} {self.unit} "
                f"[web depth {self.web_depth}]"]
        for t in sorted(self.terms, key=lambda t: t.depth):
            ancestry = f" parent={t.parent}" if t.parent else ""
            rows.append(f"      depth {t.depth} [{t.status:<9s}] "
                        f"[{t.kind:>5s}] "
                        f"{'→'.join(t.path):<28s} "
                        f"factor = {t.value(state):+.6e}{ancestry}  "
                        f"{t.note}")
        max_depth = max((t.depth for t in self.terms), default=0)
        rows.append(f"    value(≤{max_depth}) = {self.value():.10g} "
                    f"{self.unit}")
        return "\n".join(rows)


class Web(dict):
    """The interference web: nodes (ledgers) + the kernel x ← b + W(x)."""

    #  Constraint nodes (live bases) may be transiently unsolvable while
    #  the state is away from the fixed point (Koide quadratic with a
    #  negative discriminant, a missing key during seeding).  The kernel
    #  keeps the previous value and iterates on.
    _SOFT = (KeyError, ValueError, ZeroDivisionError)

    class _ReadTracker(dict):
        """State snapshot that records which ledger keys an edge reads."""

        def __init__(self, state):
            super().__init__(state)
            self.reads = set()

        def __getitem__(self, key):
            self.reads.add(key)
            return super().__getitem__(key)

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

    def validate_depth_ancestry(self):
        """Require every post searchlight term to read an explicit parent.

        The octonionic searchlight resolves through depth 3. A depth-n
        recursive web term is therefore meaningful only when it names a
        ledger at depth n-1. This turns depth from a decorative tag into
        a checked dependency chain.
        """
        for node, ledger in self.items():
            if not isinstance(ledger, Ledger):
                continue
            for term in ledger.terms:
                if term.depth <= SEARCHLIGHT_DEPTH:
                    continue
                if not term.parent:
                    raise ValueError(
                        f"{node}: depth {term.depth} term "
                        f"{'->'.join(term.path)} has no parent")
                parent = self.get(term.parent)
                if not isinstance(parent, Ledger):
                    raise ValueError(
                        f"{node}: missing ledger parent {term.parent!r}")
                if parent.web_depth != term.depth - 1:
                    raise ValueError(
                        f"{node}: depth {term.depth} term requires a "
                        f"depth {term.depth - 1} parent; {term.parent!r} "
                        f"has depth {parent.web_depth}")
                if not callable(term.factor):
                    raise ValueError(
                        f"{node}: post searchlight term must read its "
                        "parent through a state dependent edge")
                tracker = self._ReadTracker(self.state)
                try:
                    term.value(tracker)
                except self._SOFT:
                    # A transient numerical failure remains solve()'s job;
                    # key access up to that point is still observable here.
                    pass
                if term.parent not in tracker.reads:
                    raise ValueError(
                        f"{node}: depth {term.depth} term declares parent "
                        f"{term.parent!r} but does not read it")

    def solve(self, iters=300, tol=1e-15):
        """Joint fixed point of the kernel: classical Jacobi style
        fixed point iteration (a standard numerical method. No
        machine learning machinery is involved anywhere).

        FAILURES ARE LOUD.  Soft exceptions are tolerated only as
        TRANSIENTS while the state is away from the fixed point (see
        _SOFT above).  At exit, three conditions are enforced rather
        than silently absorbed: (i) the iteration converged (drift <
        tol), (ii) every ledger node holds a value, (iii) every
        ledger node evaluated cleanly on the final sweep.  Without
        (iii), a node still raising at the fixed point would carry a
        stale value unnoticed."""
        self.validate_depth_ancestry()
        drift = float("inf")
        failing = {}
        for _ in range(iters):
            new = {}
            failing = {}
            for n, led in self.items():
                if not isinstance(led, Ledger):
                    continue
                try:
                    new[n] = led.value_at(self.state)
                except self._SOFT as e:
                    failing[n] = repr(e)
                    if n in self.state:
                        new[n] = self.state[n]
            drift = max((abs(new[n] - self.state.get(n, new[n]))
                         for n in new), default=0.0)
            self.state.update(new)
            if drift < tol:
                break
        if drift >= tol:
            raise RuntimeError(
                f"web did not converge: drift {drift:.3e} >= tol {tol:.0e} "
                f"after {iters} iterations")
        if failing:
            raise RuntimeError(
                "constraint nodes still failing AT the fixed point "
                f"(stale values would be invisible): {failing}")
        unsolved = [n for n, led in self.items()
                    if isinstance(led, Ledger) and n not in self.state]
        if unsolved:
            raise RuntimeError(f"ledger nodes never solved: {unsolved}")
        return self.state


#  THE ORIENTATION RULE (multiplicity of an echo term).
#  Every closed echo cycle contributes once per orientation.  Open paths
#  contribute the trace of their channel response on the active space.
#  An identity response gives the mode count.  A current insertion
#  gives its quadratic charge trace.  The framework already uses this
#  rule:
#    • circulant amplitude B = 2ρ: the Z₃ step S and its reverse S† are
#      the two orientations of one cycle.  The published factor 2 in
#      B/A = 2|C₃̄|/(|C₁|d₁₀) IS orientation counting.
#    • bridge self-echo: a self-loop has one orientation → coefficient 1.
#    • EM action vertex: wiring.py evaluates the E8 current operator on
#      the 120 positive root modes and obtains Tr_+(Q²)=30=h∨(E8).
#  The searchlight finds the mutual electron↔quark carrier at depth 3.
#  Reading that carrier back through formed observables is the first
#  recursive layer beyond the searchlight, web depth 4.  The carrier is
#  a 2-cycle (e→q→e and q→e→q), so its multiplicity is exactly 2, the
#  number of orientations of the loop.  Hence the FORCED term −2(α/2π)²
#  in the EM bracket.  The coupling inside the loop is the physical one,
#  and the bracket is solved self-consistently.
#
#  THE VERTEX COMPOSITION RULE (depth 4 response assignment).
#  A depth 4 reading of the e↔q carrier has, at each endpoint, THE SAME
#  VERTEX WEIGHT that node's channel already carries at depth 1.  Echoes
#  compose existing vertices.  They never invent new ones.  The
#  multiplicity is the product of the two endpoint weights:
#    • coupling node (1/α): no matter vertices (backreacts on itself) →
#      orientation count alone = 2  (S,S† pair).
#    • action node (S): its depth 1 vertices are the charge trace
#      (inside α_alg = α_G₂/16, the E8 current trace vent) and the Casimir
#      C₂(26) (the −6 vent) → 8/3 × 6 = 16 = the Singh ratio.
#    • amplitude node (m): its own channel weight is the Schur
#      projector W(fund) = 1/d₁₀ (the same 1/2 that fixes B/A), and
#      the q-side traversal is the Albert composition, dim J₃(O) = 27
#      → 27/2, venting ∝ the knot's own amplitude Δₖ.
#  THE CLOSURE (tests/probes/wiring_theorem_probe.py): at searchlight
#  depth 3, F4 equivariance alone permits T=a P1+P26. At recursive depth
#  4 the path is an Albert automorphism: it preserves the unique unit and
#  its reverse is its inverse. Hence a=1, the junction is g^-1 g=1, and
#  T=P1+P26. The endpoint product 27/2 is therefore derived. The 16
#  is the unique grade of wiring_theorem_probe.py PART D.
#
#  WHY THE NODES ARE SECTORS, NOT SPECIES.  Echo paths run along
#  CHANNELS, and a channel is a single algebraic object however many
#  species vent through it.  The E₈ branching contains exactly one
#  bridge, (7,26).  It connects the G₂ sector to the F₄ sector as one
#  bifundamental, not as 3 generations × colors of separate wires.
#  Species multiplicity lives INSIDE a channel's own data: the E8
#  current index 30 and the Singh ratio 8/3 are traces within the EM
#  channel.  They never multiply the number of loops.  So the carrier
#  graph found at searchlight depth 3 has two nodes (G₂-knot sector,
#  F₄-knot sector), one edge (the EM channel carrying the bridge
#  backreaction), one 2-cycle, and two orientations: coefficient
#  exactly 2.  Counting species as nodes would double count what the
#  charge trace already counts.

#  The web: ledgers for the framework's echoed quantities.  They are
#  built in root so that every module computes WITH the echo law, not
#  before it.
#
#  READING ORDER IS EMERGENCE ORDER.  The chain runs in this order:
#    1. The coupling closes on itself         (inv_alpha cycle, here)
#    2. The anchor sets the unit              (m_e -> M_Pl, here)
#    3. One standing wave, three windings     (leptons, masses.py)
#    4. quarks = words on the generator       (walk counts, masses.py)
#    5. Every knot vents into the common mode (Sigma, gravity.py)
#    6. gravity closes on the column total    (gates on Sigma, G is step 2)
#    7. The unread remainder                  (dark_ratio = 2pi - 1)
#  The reading of an observable is its value at the shortest initial
#  part of this order that closes it.
WEB = Web()

_inv0 = 1.0 / alpha_EM                   # 512/π, the algebraic zero edge base


def _alpha_of(state):
    """The EM coupling READ FROM the web state (a known phase: scalar)."""
    return 1.0 / state["inv_alpha"] if state else math.pi / 2**9


# The searchlight closes the unique charged sector 2-cycle at depth 3.
# Depth-4 ledgers below do not recreate that loop: they explicitly read
# this carrier with their channel specific response weights.
WEB["eq_loop"] = Ledger(
    "eq_loop(e<->q carrier)", 0.0, "add", web_depth=3).echo(
    ["e<->q 2-cycle"],
    lambda s: (_alpha_of(s) / (2.0 * math.pi))**2, 3,
    "FORCED", "unique charged sector carrier L = (alpha/2pi)^2")


#  - EM coupling ledger: base 2⁹/π, additive bracket -
#  depth 1: bridge self-interference (h=1, D²=1, c_coset=0)
#  depth 4: recursive reading of the depth 3 electron↔quark carrier,
#  a state dependent cycle edge. The kernel's fixed point resums it.
#
#  EM CHANNEL FINALITY (not a global web depth cap).  Echo paths run
#  along channels between SECTOR nodes (the sector node rule), and the
#  EM channel couples only to CHARGED sectors.  The framework has
#  exactly two: the G₂-knot sector (leptons) and the F₄-knot sector
#  (quarks).  The Higgs channel is EM neutral (the condensing
#  component).  The common mode is neutral (q = 0), and no node reads
#  it back (THE COMMON MODE, below).  The EM channel's sector graph is
#  therefore two nodes and one edge.  Its only cycles are traversals of
#  the single 2-cycle, and the kernel's self-consistency already
#  resums ALL repeated traversals (α appears inside its own bracket).
#  Hence no separate (α/2π)³ edge exists (tests/probes/em_finality.py).
#  The cubic below is EXACT AND FINAL for this channel.  Other sectors
#  may carry later web depths.  The framework is fully committed to the
#  Berkeley Cs value.  If the Cs/Rb dispute resolves toward Rb, there
#  is no deeper term to add, and the sector node rule or the
#  orientation rule is falsified outright.
WEB["inv_alpha"] = (
    Ledger("1/alpha", _inv0, "add", web_depth=4)
    .echo(["bridge"],
          -_inv0 * _wiring.normalized_current_response(h_bridge), 1,
          "FORCED", "normalized bridge self-response, c_coset=0")
    .echo(["read:eq_loop"],
          lambda s: -_inv0 * 2.0 * s["eq_loop"], 4,
          "FORCED", "2-cycle orientation rule (S,S† precedent), "
                    "self-consistent α", parent="eq_loop"))

#  - lepton EM echo: charged knots vent through the EM channel with the
#  echoed coupling itself (a phase of a phase, so composite depth 2).
#  1/α's depth one edge reads the same vertex
#  (tests/probes/frame_coefficients.py F2).
WEB["lepton_EM"] = Ledger(
    "lepton_EM_echo", 1.0, web_depth=2).echo(
    ["bridge", "EM"], lambda s: -_alpha_of(s) / (2.0 * math.pi), 2,
    "FORCED", "EM vertex echo with echoed coupling")

#  THE KERNEL RUN: one recursion solves every node jointly.
WEB.solve()

inv_alpha_phys = WEB.state["inv_alpha"]
alpha_phys = 1.0 / inv_alpha_phys
QED_factor = WEB.state["lepton_EM"]

#  Closed form checks of the construction above.  The cubic is the
#  ledger's own fixed point equation written out.  The depth 1 identity
#  is its truncation.
#  • dropping the depth 4 response recovers the published depth 1 identity
assert abs(WEB["inv_alpha"].value(max_depth=1)
           - 256*(2*math.pi-1)/math.pi**2) < 1e-12
#  • the full value is the real root of the cubic
#      x³ = (512/π)[(1 − 1/(2π))x² − 1/(2π²)]
_x = inv_alpha_phys
assert abs(_x**3 - (2**9/math.pi)*((1.0 - 1.0/(2*math.pi))*_x**2
                                   - 1.0/(2.0*math.pi**2))) < 1e-6
assert abs(QED_factor - (1.0 - alpha_phys/(2.0*math.pi))) < 1e-15

# # sin²θ_W as a Web ledger #################
#
# Base: sin²θ_W = h∨(SU(3))/(h∨(G₂)+h∨(F₄)) = d₁₁/(d₁₀²+d₁₁²) = 3/13
#
# Depth-1 (FORCED): G₂ face emission through the EM channel.
#   Emission echoes carry conformal weights (apportionment law).
#   Factor: h₇ · α/(2π) = (2/5) · α/(2π)
#
# Depth-2 (FORCED): cross echo composing the G₂ emission vertex
#   with the fundamental WZW channel.  The depth 1 emission (h₇)
#   traverses the fundamental representation, suppressed by
#   h₁₀/d₁₀ = 1/d₁₁² (quantum Schur, the SAME suppression as B/A).
#   Factor: (h₇/d₁₁²) · α/(2π) = (2/45) · α/(2π)
#
# IDENTITY: h₇ + h₇/d₁₁² = h₇(1+1/d₁₁²) = (2/5)(10/9) = 4/9 = Q₀²
#   The total echo coefficient is the Koide parameter squared.
WEB["sin2W"] = (
    Ledger("sin²θ_W", float(sin2W), "add", web_depth=2)
    .echo(["G₂-face(EM)"],
          lambda s: float(h_7) * (1.0 / s["inv_alpha"]) / (2.0 * math.pi),
          1, "FORCED",
          "G₂ face emission echo, weight h₇ = d₁₀/(1+d₁₀²)")
    .echo(["G₂×WZW(fund)"],
          lambda s: (float(h_7) / d11**2)
          * (1.0 / s["inv_alpha"]) / (2.0 * math.pi),
          2, "FORCED",
          "h₇·(h₁₀/d₁₀) cross echo through fundamental WZW"))

# # bridge² as a Web ledger #################
#
# Base: Q₀² · charge_trace = (d₁₀/d₁₁)² · d₁₀³/d₁₁ = 32/27
#   Albert algebra trace norm × Dynkin Z₂ orientation factor.
#
# Depth-2 (FORCED): adjoint conformal weight at quadratic altitude.
#   The Q backreaction terms (h/K³) enter the Koide sum rule at CUBIC
#   order (leading and quadratic vanish by Z₃ symmetry).  The bridge
#   factor operates at the bridge level between the Koide spectrum
#   and the strange mass, entering at QUADRATIC order → K² not K³.
#   Channel: adjoint (1,1), weight h₁₁ = 1/d₁₀ (representation mixing).
#   Factor: h₁₁/K² = (1/d₁₀)/(d₁₀d₁₁)² = 1/72
WEB["bridge_sq"] = (
    Ledger("bridge²", float(Q0**2 * charge_trace), "mul",
           web_depth=2)
    .echo(["WZW(adj/K²)"], float(h11) / K**2, 2,
          "FORCED", "adjoint weight at quadratic altitude: h₁₁/K²"))

# Protected forgetting: PER MODE second moment of the Z₃ harmonic.
#   v_k = B cos(θ + 2πk/3),  B = √d₁₀,  θ = h₁₀
#   E[v²] = |v̂₁|² = B²/4 = 1/2      one Z₃ channel (what the bridge sees)
#   Pv²P  = Σv_k²/3 = B²/2 = 1      total, both channels: Pv²P = 2·E[v²]
# The unique (7,26) connector carries the G2/Z3 response tensored with
# identity on the irreducible F4 26.  wiring.py evaluates its normalized
# Schur trace over all 182 copies, so no independent bridge amplitude is
# introduced and the ledger consumes E[v²] = 1/2.
# The two contractions are distinct and easily conflated: "Pv²P = ½P"
# names the TOTAL, not the per mode power.  See protected_forgetting.py
# and tests/probes/circulant_normalization.py.
E_v2 = Fraction(WIRING_DATA["bridge_mode_moment"]).limit_denominator()
assert E_v2 == Fraction(1, 2)

# # G₂ coupling at the Planck scale ###############
#
# The wiring scan selects the SU(3)_3 extension current J=(3,0) by
# h_J=k/3=1.  A dimension one current has its local normalisation fixed
# by current algebra. Unlike a W-generator it introduces no independent
# coupling.  The same scan computes the asymptotic charge probability
# of the pointed (simple current) extension channel,
#
#   P = Σ_{a∈Z₃} |S_{0a}|² = Σ_{a∈Z₃} d_a²/D²_tot = 3/36 = 1/12,
#
# and evaluates the normalized current response h_J/(2π).  By the
# engine's channel coupling definition, the physical channel coupling
# is closure probability times local response:
#
#   α_G₂(M_Pl) = P·h_J/(2π) = 1/(2π·d₁₀²d₁₁) = 1/(24π).
#
# This is the same normalized dimension one current response used by
# the bridge self-echo.  The algebraic coupling above uses P and h_J.
#
# # One derivation, three forms ###################
# |Z₃|/D²_tot = 1/FPdim(C₀) holds for ANY faithfully Z₃-graded fusion
# category (each graded component carries FPdim(C)/3, Gelaki-Nikshych),
# and h∨(E₈) − C₂(26) = 24 = d₁₀·FPdim(C₀) is a dictionary identity
# asserted here.  This is one derivation written three ways.
#
# # The RG form of the same measurement ###
# 1/α = 24π is exactly equivalent, through the derived b₀ = 32/3, to
#
#   S₀ = 2π/(b₀ α) = h∨(F₄)·π²/2 = 9π²/2,
#
# i.e. to M_Pl·e^{−S₀} being the one loop RG invariant scale of G₂.
# Equivalently  h∨(E₈) − C₂(26) = b₀·h∨(F₄)/4   (24 = (32/3)(9)/4).
# Since m_e = ½·M_Pl·e^{−S₀}·(…), the m_e→M_Pl→G chain measures this.
# The derived M_Pl sits +22 ppm above CODATA.  The (smaller) S₀ that
# would put M_b on CODATA's Planck mass implies 1/α = 75.398187,
# against 24π = 75.398224, a relative agreement of 4.9e-7.  That 4.9e-7
# is the +22 ppm Planck mass agreement divided by S₀ = 44.41.  It is the
# one G comparison rewritten in RG variables.  The only other downstream quantity that
# moves with α is the Higgs.  Perturbing alpha moves the transported
# Higgs running coordinate (through lambda(M_Pl) = -N_bridge*alpha²*E[v²]).
# The heat kernel sum Sigma is a count and does not move with α
# (boundary.py, THE AMOUNT'S COUPLING).
#
# The F₄/lepton action is derived below by one loop dimensional
# transmutation, S_lepton = 2π/(b₀α). The later G₂ baryogenesis
# reading is not obtained by equating a static texture energy with its
# target space volume.  It is the same action, transported to the G₂
# face by the h-dual response metric.  The exact web derivation
# follows below.
_D2_eff = d10**2 * d11                       # FPdim(C₀) = 12
assert hv_E8 - C2_26 == d10 * _D2_eff        # 24 = 2·12
assert abs(WIRING_DATA["D2_tot"] - 3 * _D2_eff) < 1e-12
assert abs(WIRING_DATA["pointed_weight"] - 1.0/_D2_eff) < 1e-12
alpha_G2_Pl = WIRING_DATA["alpha_G2_Pl"]
assert abs(alpha_G2_Pl - 1.0/(24.0*math.pi)) < 1e-15
assert Fraction(hv_E8 - C2_26) == Fraction(32, 3) * hv_F4 / 4  # RG form

# The nonminimal coupling of the bridge label.  The amount is
# minimally coupled and does not carry it (boundary.py).
# The bridge moment is generated by the W6 Schur trace in wiring.py.
xi_bridge = alpha_G2_Pl * float(E_v2) * float(h_bridge)  # = 1/(48π)

# # Instanton actions (dimensional transmutation) ########
#
# Standard form: S₀ = 2π / (b₀ · α_G₂(M_Pl))
#   b₀ = (11/3)C_A − (4/3)ΣT_R
#      C_A = h∨(G₂) = 4
#      Matter: d₁₁ = 3 Dirac fermions in the 7 of G₂ (one per generation)
#      T(7) = 1  (Dynkin index of the fundamental rep of G₂)
#      ΣT_R = 3 × 1 = 3
#      b₀ = (11/3)(4) − (4/3)(3) = 44/3 − 4 = 32/3
#   α_G₂(M_Pl) = 1/(24π)  (P·h_J/(2π), the normalized current response above)
#   S₀ = 2π / ((32/3) · 1/(24π)) = 2π · 24π · 3/32 = 9π²/2
#      = h∨(F₄) · π²/2 = d₁₁² · π²/2
#
# The exponential form exp(−2π/(b₀α)) is the standard one loop
# dimensional transmutation result.  The code computes that path
# directly.  The F₄ expression is then an asserted equality, not the
# source.
T_7 = Fraction(1, 1)
n_dirac_G2 = d11
# The running reads the light G2-charged rows of the one boundary list
# (boundary.py), the list the induced stiffness reads (gravity.py).
b0_G2 = _boundary.b0(_boundary.rows(d10, d11, n7, n26))
assert b0_G2 == Fraction(11, 3) * hv_G2 - Fraction(4, 3) * n_dirac_G2 * T_7
assert b0_G2 == Fraction(32, 3)
S_lepton = 2.0 * math.pi / (float(b0_G2) * alpha_G2_Pl)
assert abs(S_lepton - hv_F4 * math.pi**2 / 2.0) < 1e-12
# The G₂ lever arm ends at μ* = M_Pl e^{−(S_lepton − C₂(26))}.  One
# loop running from 1/α_G₂(M_Pl) over that arm gives exactly π/32, the
# value of alpha_G2_WZW above.
_inv_alpha_G2_mu_star = (1.0 / alpha_G2_Pl
                         - float(b0_G2) / (2.0 * math.pi)
                         * (S_lepton - float(C2_26)))
assert abs(_inv_alpha_G2_mu_star - 1.0 / alpha_G2_WZW) < 1e-9, \
    "π/32 is the one loop endpoint of the G₂ lever arm"
assert abs(_inv0 - 1.0 / alpha_EM) == 0.0

# The action has two face readings, not two unrelated normalizations.
# The F4-face coefficient is the one loop b0 above. Transport to the G2
# face uses the already selected h-dual bridge metric. Its response
# coefficient closes independently on the embedding invariant 24:
#
#   b_G2-face = b0(G2) h∨(F4)/h∨(G2) = 24
#             = h∨(E8) - C2(26) = d10 FPdim(C0).
#
# Therefore the recursive G2-face reading is
#   S_baryo = 2pi/(24 alpha_G2) = 2pi^2 = (4/9) S_lepton.
# This is the depth 4 web action that baryogenesis uses.
b_baryo_face = b0_G2 * Fraction(hv_F4, hv_G2)
assert b_baryo_face == hv_E8 - C2_26 == d10 * _D2_eff == 24
S_baryo_from_response = 2.0 * math.pi / (
    float(b_baryo_face) * alpha_G2_Pl)
assert abs(S_baryo_from_response - Fraction(hv_G2, hv_F4)
           * S_lepton) < 1e-12

WEB["S_lepton"] = Ledger(
    "S_lepton(action,F4-face)", S_lepton, "add", web_depth=3)
WEB["S_baryo"] = Ledger(
    "S_baryo(action,G2-face)",
    lambda s: float(Fraction(hv_G2, hv_F4)) * s["S_lepton"],
    "add", web_depth=4)
WEB.solve()
S_baryo = WEB.state["S_baryo"]
assert abs(S_baryo - S_baryo_from_response) < 1e-12
assert abs(S_baryo - 2 * math.pi**2) < 1e-12

# ═══════════════════════════════════════════════════════════════════════
#  ANCHOR INVERSION: m_e → M_Pl → G  (all dimensionless data in hand)
# ═══════════════════════════════════════════════════════════════════════
#
#  The full electron formula (every factor FORCED):
#    m_e = ½ M_Pl e^(−S_lepton) Δ_e² · (1 − α/2π) · (1 − (27/2)(α/2π)²Δ_e)
#  with Δ_e = 1 + √d₁₀ cos(h₁₀ + 2π/3).  Inverting for M_Pl:

# Trace lift normalisation, ONE symbol for the prefactor in every
# A = C_EFF·M_Pl·e^{−S₀} amplitude (the same quantity appears in
# higgs_ew.py's Lambda_conf_GeV).  C_EFF is an identification:
# C_EFF = 1/d₁₀ = W(fund), the SAME quantum Schur channel weight
# used by the B/A channel map (octonions.py projector lemma).  What is
# identified is that the amplitude prefactor IS the fundamental
# channel's weight, not the value of the weight itself.  The frame
# reads the closure's whole common mode, C_EFF·P Δ² P = 1
# (tests/probes/frame_coefficients.py F1).
C_EFF = 1.0 / d10                        # = 1/2
assert abs(C_EFF * (1 + 2 * E_v2) - 1) < 1e-15

_Delta_e = 1.0 + math.sqrt(d10) * math.cos(float(h10) + 2*math.pi/3)
_vent_e = 1.0 - (d11**3 / d10) * (alpha_phys / (2*math.pi))**2 * _Delta_e
_F_e = C_EFF * math.exp(-S_lepton) * _Delta_e**2 * QED_factor * _vent_e

M_Pl_MeV = M_E_ANCHOR_MEV / _F_e
M_Pl_GeV = M_Pl_MeV / 1e3
# IDENTIFICATION: the boundary.  The inversion computes M_b, the scale
# at which alpha_G2 = 1/(24pi) holds, and reads it as the full Planck
# mass sqrt(hbar c/G).  M_b cancels from v_EW/m_e and from alpha_s, so
# G is its only test.  Under the induced ledger's 1/G = Sigma/(2pi t0),
# the same identification gives t0 = Sigma/(2pi) M_b^-2.  The regulator
# is the boundary written in proper time form, not a second unknown.
# On the lattice of the engine's shares, the only fusion closed set
# that contains every closure is the whole, c = 1
# (tests/probes/boundary_lattice.py).  On the depth 4 web the common
# mode reads it as G M_b^2 = 3 Q_l - 1 (THE COMMON MODE, below).

#  The unit as a graph node (a live base: the anchor inversion read
#  from the state).  Perturbing the web and solving again returns M_Pl
#  to the same fixed point.  The constant is stored nowhere.  It is the
#  fixed point of the recursion.
#  The lepton circulant's gaps Δₖ = 1 + √d₁₀ cos(h₁₀ + 2πk/3), and the
#  depth 4 term each rate carries, −(d₁₁³/d₁₀)·L·Δₖ, have one
#  definition.  The anchor, the lepton ledgers (masses.py) and the
#  common mode all read it.
LEPTON_GAPS = {name: 1.0 + math.sqrt(d10) * math.cos(float(h10) + 2*math.pi*k/3)
               for k, name in ((1, "e"), (2, "mu"), (0, "tau"))}
assert LEPTON_GAPS["e"] == _Delta_e


def lepton_carrier_term(state, gap):
    """The depth 4 term on a lepton rate: the carrier read at W(fund)·dim J₃(O) = 27/2, ∝ Δₖ."""
    return -(d11**3 / d10) * state["eq_loop"] * gap


def _F_e_of(state):
    a = _alpha_of(state)
    qed = 1.0 - a / (2*math.pi)
    # the electron's depth 4 term reads the carrier node, as its ledger does
    vent = 1.0 + lepton_carrier_term(state, _Delta_e)
    return C_EFF * math.exp(-S_lepton) * _Delta_e**2 * qed * vent

WEB["M_Pl_MeV"] = Ledger("M_Pl_MeV",
                         lambda s: M_E_ANCHOR_MEV / _F_e_of(s),
                         "mul", "MeV", web_depth=4)
WEB.solve()
assert abs(WEB.state["M_Pl_MeV"] / M_Pl_MeV - 1.0) < 1e-14

#  THE COMMON MODE: gravity's strength G M_b² as a node of the web.
#  The background is not a separate object.  It is the field's common
#  mode, the component that carries no label.  With every interference
#  accounted for, the objects that close appear in it as a bias.  That
#  bias, measured, is gravity.  The node is a COUNT, not a search.  It
#  carries no weight of its own.  On the depth four web it reads the
#  closure its edges belong to.  G = common_mode·ħc/M_b².
#    base 0     the vacuum is the empty field, with no amplitude of its
#               own
#    edges, 4   P v P = 0 leaves no first moment. Each relative channel
#               q = 1, 2 of the lepton circulant vents its power |v_q|²
#               on the depth four rates.  On the searchlight's circulant
#               each is E[v²] = 1/2.
#  So G M_b² = P v² P = 3 Q_l − 1, the lepton circulant's Koide moment
#  at depth four: 1 − 7.50 (α/2π)², with Q_l 5.06 ppm below 2/3.  The
#  rival readings of the node fail (tests/probes/common_mode_count.py).
#  IDENTIFIED: that gravity's strength is this node's value.  The count
#  is fixed without using any measured G.
#  WHAT ENDS IT.  A G at 1e-5 that excludes 6.673936e-11.

def _lepton_z3_power(state, q):
    """|v_q|² of the lepton circulant on its depth four rates: harmonic power per channel over the q = 0 power."""
    rates = [g * g * (1.0 + lepton_carrier_term(state, g)) for g in LEPTON_GAPS.values()]
    amp = [math.sqrt(x) for x in rates]
    vq = sum(a * cmath.exp(-2j * math.pi * q * k / 3) for k, a in enumerate(amp)) / 3
    v0 = sum(amp) / 3
    return abs(vq) ** 2 / v0 ** 2


WEB["common_mode"] = (
    Ledger("common_mode(G M_b²)", 0.0, "add", web_depth=4)
    .echo(["Z3(q=1)"], lambda s: _lepton_z3_power(s, 1), 4,
          "FORCED", "relative channel q = 1, its power |v_1|² on the "
                    "depth four circulant", kind="vent", parent="eq_loop")
    .echo(["Z3(q=2)"], lambda s: _lepton_z3_power(s, 2), 4,
          "FORCED", "relative channel q = 2, its power |v_2|² on the "
                    "depth four circulant", kind="vent", parent="eq_loop"))
WEB.solve()
COMMON_MODE = WEB.state["common_mode"]
_L_now = WEB.state["eq_loop"]
_Q_l = (sum(g * g * (1 + lepton_carrier_term(WEB.state, g)) for g in LEPTON_GAPS.values())
        / sum(g * math.sqrt(1 + lepton_carrier_term(WEB.state, g)) for g in LEPTON_GAPS.values()) ** 2)
assert abs(_lepton_z3_power({"eq_loop": 0.0}, 1) - float(E_v2)) < 1e-14
assert abs(COMMON_MODE - (3 * _Q_l - 1)) < 1e-14
assert abs((COMMON_MODE - 1.0) / _L_now + 7.502) < 0.001

#  Newton's constant reads the unit through the common mode:
#  G = common_mode·ħc/M_b².  G_DEPTH_ONE is ħc/M_b² in CODATA's units,
#  the factor the common mode multiplies.
G_DEPTH_ONE = G_CODATA * (M_PL_CODATA_GEV / M_Pl_GeV)**2
G_PRED = COMMON_MODE * G_DEPTH_ONE
#  sanity: the derived Planck mass sits within G's experimental band
assert abs(M_Pl_GeV / M_PL_CODATA_GEV - 1.0) < 3.3e-5  # 3σ_G/2
assert abs(G_PRED - 6.673936e-11) < 1e-15  # frozen prediction value
assert abs(S_lepton - Fraction(9, 4) * 2 * math.pi**2) < 1e-12

# Friedmann critical density scale for the unreduced Planck mass:
# ρ_crit = 3 H₀²/(8πG) = (3/8π) M_b² H₀²/common_mode.
# Calling this ρ_Λ would silently assume Omega_Lambda = 1. The current
# dark energy density is Omega_Lambda times this value (dark_sector.py).
RHO_CRIT_MAPPED = (3.0 * H_0_GEV**2 * M_Pl_GeV**2
                   / (8.0 * math.pi * COMMON_MODE))

# # Electroweak scale (v_EW derivation chain) #########
#
# v_EW = M_Pl · exp(−S_eff),  where S_eff = S₀ − C₂(26) + δS
#
#   S₀ = 9π²/2:  lepton instanton action (above)
#   −C₂(26) = −d₁₀d₁₁ = −6:  functional determinant shift from quarks
#       living in the 26 of F₄ (Casimir of 26-dimensional representation)
#   δS = Tr_+(Q²) · α_EM / (2π): boundary current backreaction
#       Q E_α = (α·q)E_α for a unit Cartan current q. wiring.py sums one
#       representative of each E8 root pair and derives
#       Tr_+(Q²) = Σ_{α>0}(α·q)² = 30 = h∨(E₈), independent of q.
#       α_EM = π/512 (derived above)
#       δS = 30 · (π/512) / (2π) = 15/512
#
N_vertex = vertex_charge_trace           # not a count: = Tr_+(Q²) = 30
delta_S = N_vertex * alpha_EM / (2 * math.pi)  # = 15/512

#  - quark layer action ledger: ALL departures from the lepton action
#  are echoes.  The F4 knot vents its Casimir into the instanton
#  background (−C₂(26)).  The E8 current trace vents through the EM
#  channel at its emergence reading (+Tr_+(Q²)α/2π).  One law, two
#  channels.
WEB["S_quark"] = (
    Ledger("S_quark(action)", S_lepton, "add", web_depth=4)
    .echo(["F4-knot(26)"], -float(C2_26), 1,
          "FORCED", "Casimir vent of the 26 into the instanton background",
          kind="vent")
    .echo(["EM(alg)"], delta_S, 1,
          "FORCED", f"E8 current trace {N_vertex} = 15/512",
          kind="vent")
    # Depth-4 response: the action's own
    # depth 1 vertices compose across the e↔q loop,
    # charge_trace × C₂(26) = (8/3)·6 = 16 (= the Singh ratio = d₁₀⁴).
    # The endpoint product is exactly 16. The reverse path composes with
    # the Albert transport as g^-1 g, fixing junction normalization to 1.
    # Closes v_EW to G_F at 1e-7: 246.219645 vs 246.219640 from
    # G_F = 1.1663788(6)e-5 GeV⁻² (+0.02 ppm, 0.08σ_GF).
    .echo(["read:eq_loop(Singh)"],
          lambda s: -float(charge_trace * C2_26)
          * s["eq_loop"], 4,
          "FORCED", "depth 4 response: charge_trace × C₂(26) = 16",
          parent="eq_loop"))

#  - Higgs boundary ledger: the bridge echoes into the scalar channel -
#  Depth-2 vent: the fundamental channel share of the bridge vents
#  into the condensate and does not reach the scalar.  This reduces the
#  edge by the factor (1 − h₁₀).  The 7 of G₂ decomposes as 3 ⊕ 3̄ ⊕ 1
#  under the confined SU(3).  The triality charged share cannot cross a
#  confining interface (N-ality superselection, words.py).  A selection
#  rule carries no altitude factor, so this is a vent, like the Casimir
#  vent above.  The blocked share in the weight metric is h₁₀.  No Higgs
#  mass is read from the edge.  That reading would need a mixed RGE
#  convention, which is not canonical.
WEB["lambda_MPl"] = (
    Ledger("lambda(M_Pl)", 0.0, "add", web_depth=2)
    .echo(["bridge", "Higgs"], -(n7 * n26) * alpha_G2_Pl**2 * float(E_v2), 1,
          "FORCED", "N_bridge·α_G₂²·E[v²] bridge→Higgs echo")
    .echo(["WZW(fund)-vent"],
          (n7 * n26) * alpha_G2_Pl**2 * float(E_v2) * float(h10), 2,
          "FORCED", "fundamental share vent (1−h₁₀)",
          kind="vent"))

#  rerun the kernel so the new nodes join the solved state
WEB.solve()
S_quark = WEB.state["S_quark"]           # = 9π²/2 − 6 + 15/512

# # Koide Z₃ amplitude ratio (octonionic/MTC channel map) #
#
# The Z₃ circulant structure gives √mₖ = A[1 + (B/A)cos(θ + 2πk/3)].
# The Koide parameter Q = 1/3 + |B/A|²/6  (Z₃ circulant identity).
#
# The B/A channel map uses two structures:
#   (a) Octonionic CG from G₂ = Aut(O), 7 → 3⊕3̄⊕1:
#       |C₃̄|/|C₁| = √2 = √d₁₀  (Fano cross product, Schur's lemma)
#   (b) Quantum Schur suppression in SU(3)₃ MTC:
#       1/d(fund) = 1/d₁₀  (off diagonal channel, quantum dimension)
#   Assembly (two cycle orientations × CG × 1/d):
#       B/A = 2|C₃̄|/(|C₁|·d₁₀) = 2√2/2 = √2 = √d₁₀
#
# Then Q₀ = 1/3 + d₁₀/6 = 2/3 = d₁₀/d₁₁.
#
# The common mode reads the same number: on the searchlight's circulant
# its count is P v² P = (B/A)²/2 = 3Q − 1, so B/A = √d₁₀, Q₀ = 2/3 and
# the whole second moment G M_b² = 1 are one identification
# (tests/probes/boundary_whole.py B8).
#
# Credit: Koide (1982) for the empirical relation.
#         Brannen (2006) for the parameterisation.
# Derivation: Z₃ circulant algebra + octonionic CG + quantum Schur lemma.
#
BA_ratio = math.sqrt(d10)               # |B/A| = √d₁₀ = √2


# ═══════════════════════════════════════════════════════════════════════
#  PDG REFERENCE VALUES  (comparison only, not inputs)
# ═══════════════════════════════════════════════════════════════════════
#  Kept in reference.py (its rules and declared entry list are there).
#  Reexported here so consumers can read the references through
#  `interference.root`.  Both names refer to the same dict objects.

from .reference import (PDG_MASSES, PDG_CKM, NUFIT_PMNS,       # noqa: F401
                        NUFIT_OSCILLATION, PDG_EW, PDG_COSMO,
                        PLANCK_2018, ACT_DR6, SUPER_K,
                        PDG_LIGHT_RATIOS, ALPHA_INV_REF)


def rk4_run(f, y, t0, t1, n_steps):
    """Generic RK4 integrator dy/dt = f(y), float or ndarray state.
    Single shared implementation, used by couplings and gravity."""
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
    print(f"    δ = 1/(d₁₀d₁₁²) = {delta}  [since d₁₁²−d₁₀³ = 1]")
    print(f"    h(7) = d₁₀/(1+d₁₀²) = {h_7},  h(26) = d₁₀d₁₁/(1+d₁₁²) = {h_26}")
    print(f"    h_bridge = {h_bridge},  c_coset = {c_coset}")

    print(f"\n  Weak mixing angle:")
    print(f"    sin²θ_W = d₁₁/(d₁₀²+d₁₁²) = {d11}/{d10**2+d11**2} = {float(sin2W):.5f}")

    print(f"\n  Couplings:")
    print(f"    α_G₂(μ* = 253.5 GeV) = π/32,  α_EM = π/512 (zero edge base)")
    print(f"    1/α(0) = (512/π)(1−1/(2π)−2(α/2π)²) = {inv_alpha_phys:.9f}")
    print(f"            [depth 1 truncation: (512/π)(1−1/(2π)) = 137.036439]")

    print(f"\n  Action values (layer resolved web readings):")
    print(f"    S_lepton = d₁₁²π²/2 = {S_lepton:.4f}")
    print(f"    S_baryo  = (4/9)S_lepton = {S_baryo:.4f}  [depth 4 G₂ face]")

    return {
        'd10': d10, 'd11': d11, 'n7': n7, 'n26': n26,
        'wiring': WIRING_DATA.copy(),
        'h_J': WIRING_DATA['h_J'],
        'pointed_weight': WIRING_DATA['pointed_weight'],
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
        'b0_G2': float(b0_G2),
        'S_quark': S_quark,
        'vertex_charge_trace': vertex_charge_trace,
        'N_vertex': N_vertex,
        'BA_ratio': BA_ratio,
        'h_dual': {'G2': hv_G2, 'F4': hv_F4, 'E8': hv_E8},
        'C2_fund_G2': C2_7, 'C2_fund_F4': C2_26,
        'dim_G2_fund': n7, 'dim_F4_fund': n26,
    }
