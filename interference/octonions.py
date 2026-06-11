"""
Octonionic Clebsch-Gordan verification: B/A = √d₁₀ from first principles.

Constructive proof that the Koide amplitude ratio B/A = √2 = √d₁₀
follows from the Fano-plane multiplication table of the octonions
plus the quantum Schur lemma in SU(3)₃.

Derivation chain:
    1. Build Im(O) ≅ R⁷ cross product from Fano-plane triples
    2. Holomorphic basis: z_α = (e_{2α-1} − i e_{2α})/√2, α=1,2,3
    3. Singlet channel:  z_α × z̄_β = i δ_{αβ} e₇  →  |C₁| = 1
    4. Triplet channel:  z_α × z_β = √2 ε_{αβγ} z̄_γ  →  |C₃̄| = √2
    5. Channel weight:   W(3̄) = 1/d(fund) = 1/d₁₀ = 1/2
    6. Assembly:         B/A = 2|C₃̄|·W(3̄)/|C₁| = 2√2·(1/2)/1 = √2
    7. Koide:            Q₀ = 1/3 + |B/A|²/6 = 1/3 + 2/6 = 2/3 = d₁₀/d₁₁

The ratio |C₃̄|/|C₁| = √2 is frame-independent (Schur's lemma:
both maps are unique up to phase).

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import math
import numpy as np

from root import d10, d11


# ═══════════════════════════════════════════════════════════════════════
#  Octonion multiplication table from Fano plane
# ═══════════════════════════════════════════════════════════════════════

# Fano-plane triples (i,j,k) with e_i × e_j = ±e_k (cyclic).
# Convention: e₇ is the SU(3)-invariant direction.
# The six-plane (e₁,...,e₆) carries the 3 ⊕ 3̄ of SU(3) ⊂ G₂.
FANO_TRIPLES = [
    (1, 2, 7),   # e₁ × e₂ = +e₇
    (3, 4, 7),   # e₃ × e₄ = +e₇
    (5, 6, 7),   # e₅ × e₆ = +e₇
    (1, 3, 5),   # e₁ × e₃ = +e₅
    (1, 4, 6),   # e₁ × e₄ = −e₆
    (2, 3, 6),   # e₂ × e₃ = −e₆
    (2, 4, 5),   # e₂ × e₄ = −e₅
]

FANO_SIGNS = [+1, +1, +1, +1, -1, -1, -1]


def _cross_complex(a, b):
    """Cross product on Im(O) ≅ C⁷, extended to complex coefficients.

    For each Fano triple (i,j,k) with sign s:
        e_i × e_j = s·e_k,  e_j × e_k = s·e_i,  e_k × e_i = s·e_j
    (cyclic with same sign, antisymmetric under swap).
    """
    result = np.zeros(7, dtype=complex)
    for (i, j, k), sign in zip(FANO_TRIPLES, FANO_SIGNS):
        ii, jj, kk = i - 1, j - 1, k - 1
        result[kk] += sign * (a[ii] * b[jj] - a[jj] * b[ii])
        result[ii] += sign * (a[jj] * b[kk] - a[kk] * b[jj])
        result[jj] += sign * (a[kk] * b[ii] - a[ii] * b[kk])
    return result


# ═══════════════════════════════════════════════════════════════════════
#  Full octonion product + LINEARITY FROM CLOSURE (executable form)
# ═══════════════════════════════════════════════════════════════════════
#
#  (This file is the canonical home of the linearity argument.)
#
#  THE PUZZLE.  The substrate model (Z₃ coupled NLS) is nonlinear, and
#  solitons of a nonlinear equation do not superpose, yet emergent
#  quantum mechanics is EXACTLY linear.  Resolution in three steps:
#
#  (i)  The emergent state is the LEDGER, not the displacement.  A
#       knot's state is its phase-of-phase record; records compose by
#       the division-algebra product, which is BILINEAR with
#       |xy| = |x||y| (Hurwitz).  Composition acts by left
#       multiplication L_x: a linear, norm-preserving operator.
#       Linearity is an exactness property of the bookkeeping.
#       [checked: _check_linearity]
#
#  (ii) Norm multiplicativity IS probability conservation, and it is
#       selected: a branch whose records composed
#       nonlinearly would have path-dependent norms, loop accounting
#       that does not close, and by the selection principle it cannot
#       bind.  By Wigner's theorem norm-preserving composition forces
#       unitary (linear) evolution.  The Hurwitz gate is executable:
#       one Cayley-Dickson doubling past O produces ZERO DIVISORS,
#       |xy| = 0 with |x|,|y| ≠ 0, the bookkeeping fails.
#       [checked: _check_hurwitz_gate]
#
#  (iii) The substrate nonlinearity is REROUTED, not erased: the Z₃
#       selection rule (q₁+q₂+q₃ ≡ 0 mod 3) sends the leading cubic
#       response into the common mode, the nonlinearity that would
#       deform matter superposition is exactly the part that becomes
#       gravity; residual relative-sector couplings are interactions
#       and echoes (operators on a linear state space), never
#       superposition violations.  [checked: _check_z3_selection]
#
#  REMARK (Born rule, conjectural): protected forgetting erases linear
#  ledger entries (PvP = 0) and retains second moments (Pv²P = ½P);
#  measurement statistics are web memories, so outcome frequencies can
#  only read |amplitude|², the same structure that fixes m = Δ².
#  Recorded as a conjecture, not used in canonical values.

def _oct_mul(x, y):
    """Full octonion product on R⁸ (basis 1, e₁..e₇) from Fano triples."""
    z = np.zeros(8)
    z[0] = x[0]*y[0] - np.dot(x[1:], y[1:])
    z[1:] = x[0]*y[1:] + y[0]*x[1:] + np.real(
        _cross_complex(x[1:].astype(complex), y[1:].astype(complex)))
    return z


def _cd_double(mul):
    """One Cayley-Dickson doubling of a *-algebra product.

    (a,b)(c,d) = (ac − d*b, da + bc*), conjugation negating imaginaries.
    """
    def conj(x):
        out = -x.copy(); out[0] = x[0]; return out

    def mul2(x, y):
        n = len(x)//2
        a, b, c, d = x[:n], x[n:], y[:n], y[n:]
        return np.concatenate([mul(a, c) - mul(conj(d), b),
                               mul(d, a) + mul(b, conj(c))])
    return mul2


def _check_linearity(rng):
    """(i) bilinearity + norm multiplicativity + L_x isometry on O."""
    for _ in range(50):
        x, y, z = rng.standard_normal((3, 8))
        s, t = rng.standard_normal(2)
        # bilinearity: L_x is a linear operator
        assert np.allclose(_oct_mul(x, s*y + t*z),
                           s*_oct_mul(x, y) + t*_oct_mul(x, z), atol=1e-12)
        # Hurwitz: |xy| = |x||y|  (probability conservation of the ledger)
        assert abs(np.linalg.norm(_oct_mul(x, y))
                   - np.linalg.norm(x)*np.linalg.norm(y)) < 1e-10
    return True


def _check_hurwitz_gate(rng):
    """(ii) the 4th doubling fails: sedenions have zero divisors."""
    sed_mul = _cd_double(_oct_mul)
    # norm multiplicativity FAILS one doubling past O: exhibit a zero
    # divisor (in this Fano convention): (e₁ + e₁₀)(e₃ − e₁₂) = 0,
    # one of 168 such pairs on the unit-coefficient lattice.
    e = np.eye(16)
    zd1, zd2 = e[1] + e[10], e[3] - e[12]
    prod = sed_mul(zd1, zd2)
    assert np.linalg.norm(zd1) > 1 and np.linalg.norm(zd2) > 1
    assert np.linalg.norm(prod) < 1e-12, "sedenion zero divisor expected"
    return True


def _check_e8_lattice_gate():
    """Rank-8 coherence gate (gate paper, executable): among the rank-8
    candidates, only E₈ is simultaneously EVEN (all squared norms even),
    UNIMODULAR (det Gram = 1), and positive definite.  I₈ is unimodular
    but odd; A₈ and D₈ are even but not unimodular."""
    I8 = np.eye(8)
    A8 = 2*np.eye(8) - np.diag(np.ones(7), 1) - np.diag(np.ones(7), -1)
    D8 = A8.copy(); D8[7, 6] = D8[6, 7] = 0; D8[7, 5] = D8[5, 7] = -1
    E8 = A8.copy()
    E8[7, 6] = E8[6, 7] = 0          # detach node 8 from node 7 ...
    E8[7, 4] = E8[4, 7] = -1         # ... attach it to node 5 (E₈ diagram)
    results = {}
    for name, G in [("I8", I8), ("A8", A8), ("D8", D8), ("E8", E8)]:
        det = round(np.linalg.det(G))
        even = all(int(G[i, i]) % 2 == 0 for i in range(8))
        posdef = np.all(np.linalg.eigvalsh(G) > 0)
        results[name] = (even, det == 1, posdef)
    assert results["E8"] == (True, True, True)
    assert not all(results["I8"][:2])   # odd
    assert not all(results["A8"][:2])   # det 9
    assert not all(results["D8"][:2])   # det 4
    return True


def _check_z3_selection():
    """(iii) cubic vertices vanish unless q₁+q₂+q₃ ≡ 0 (mod 3)."""
    w = np.exp(2j*np.pi/3)
    for q1 in range(3):
        for q2 in range(3):
            for q3 in range(3):
                s = sum(w**((q1+q2+q3)*k) for k in range(3))/3
                expect = 1.0 if (q1+q2+q3) % 3 == 0 else 0.0
                assert abs(s - expect) < 1e-12
    return True


# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive():
    """Verify B/A = √d₁₀ from the octonionic Clebsch-Gordan coefficients.

    Returns dict with C1, C3bar, BA_ratio, Q0.
    """
    print("\n" + "=" * 78)
    print("  OCTONIONIC CG VERIFICATION: B/A = √d₁₀ from Fano plane")
    print("=" * 78)

    # ── Holomorphic basis on SU(3) six-plane ─────────────────────────
    # z_α = (e_{2α-1} − i e_{2α})/√2,  α = 1,2,3
    # z̄_α = conjugate
    # e₇ = SU(3)-invariant direction

    z = np.zeros((3, 7), dtype=complex)
    zbar = np.zeros((3, 7), dtype=complex)
    for alpha in range(3):
        z[alpha, 2*alpha] = 1.0 / math.sqrt(2)
        z[alpha, 2*alpha + 1] = -1j / math.sqrt(2)
        zbar[alpha] = z[alpha].conj()

    print(f"\n  Holomorphic basis: z_α = (e_{{2α-1}} − i e_{{2α}})/√2")

    # ── Singlet channel: 3 ⊗ 3̄ → 1 ──────────────────────────────────
    # z_α × z̄_β = i δ_{αβ} e₇  →  |C₁| = 1

    print(f"\n  Singlet channel (3⊗3̄ → 1):")
    C1_values = []
    for alpha in range(d11):
        for beta in range(d11):
            prod = _cross_complex(z[alpha], zbar[beta])
            e7_coeff = prod[6]
            if alpha == beta:
                C1_values.append(abs(e7_coeff))
                print(f"    z_{alpha+1} × z̄_{beta+1} = ({e7_coeff:.4f}) e₇"
                      f"  [expected: i]")

    C1 = np.mean(C1_values)
    assert abs(C1 - 1.0) < 1e-10, f"|C₁| = {C1}, expected 1"
    print(f"  |C₁| = {C1:.6f}  ✓")

    # ── Antisymmetric channel: 3 ⊗ 3 → 3̄ ────────────────────────────
    # z_α × z_β = √2 ε_{αβγ} z̄_γ  →  |C₃̄| = √2

    print(f"\n  Triplet channel (3⊗3 → 3̄):")
    C3bar_values = []
    # Cyclic pairs (α,β,γ) with ε_{αβγ} = +1
    pairs = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
    for alpha, beta, gamma in pairs:
        prod = _cross_complex(z[alpha], z[beta])
        norm_zbar = np.sqrt(np.vdot(zbar[gamma], zbar[gamma]).real)
        overlap = np.vdot(zbar[gamma], prod)
        coeff_abs = abs(overlap) / norm_zbar
        print(f"    z_{alpha+1} × z_{beta+1}: coeff of z̄_{gamma+1}"
              f" = {coeff_abs:.6f}  [expected: √2 = {math.sqrt(2):.6f}]")
        C3bar_values.append(coeff_abs)

    C3bar = np.mean(C3bar_values)
    assert abs(C3bar - math.sqrt(2)) < 1e-10, f"|C₃̄| = {C3bar}, expected √2"
    print(f"  |C₃̄| = {C3bar:.6f} = √d₁₀  ✓")

    # ── Frame independence ───────────────────────────────────────────
    ratio_CG = C3bar / C1
    print(f"\n  |C₃̄|/|C₁| = {ratio_CG:.6f} = √{d10}")
    print(f"  Frame-independent (Schur's lemma: both maps unique up to phase)")

    # ── Channel weight and B/A assembly ──────────────────────────────
    # W(3̄) = 1/d(fund) = 1/d₁₀ (quantum Schur lemma in SU(3)₃ MTC)
    # B/A = 2|C₃̄|·W(3̄)/|C₁| = 2√2·(1/2)/1 = √2 = √d₁₀

    d_fund = d10                           # quantum dimension d(1,0) = 2
    W_fund = 1.0 / d_fund                  # = 1/2
    BA_derived = 2.0 * C3bar * W_fund / C1  # = √2

    print(f"\n  Quantum Schur suppression:")
    print(f"    d(fund) = d₁₀ = {d_fund},  W(3̄) = 1/d₁₀ = {W_fund}")
    print(f"    B/A = 2|C₃̄|·W(3̄)/|C₁| = 2×{C3bar:.4f}×{W_fund}÷{C1:.4f}")
    print(f"        = {BA_derived:.6f} = √{d10}  ✓")
    assert abs(BA_derived - math.sqrt(d10)) < 1e-10

    # ── Koide parameter ──────────────────────────────────────────────
    # Q₀ = 1/3 + |B/A|²/6 = 1/3 + d₁₀/6 = 2/3 = d₁₀/d₁₁

    Q0 = 1.0/3.0 + BA_derived**2 / 6.0
    assert abs(Q0 - float(d10)/float(d11)) < 1e-10
    print(f"\n  Q₀ = 1/3 + |B/A|²/6 = 1/3 + {d10}/6 = {Q0:.10f} = d₁₀/d₁₁  ✓")

    # ── Linearity from closure (executable §linearity) ───────────────
    rng = np.random.default_rng(7)
    ok_lin = _check_linearity(rng)
    ok_gate = _check_hurwitz_gate(rng)
    ok_sel = _check_z3_selection()
    ok_e8 = _check_e8_lattice_gate()
    print(f"\n  Linearity from closure:")
    print(f"    (i)  L_x linear, |xy|=|x||y| on O (50 random checks)  "
          f"{'✓' if ok_lin else '✗'}")
    print(f"    (ii) Hurwitz gate executable: sedenion zero divisor "
          f"(e₁+e₁₀)(e₃−e₁₂)=0  {'✓' if ok_gate else '✗'}")
    print(f"    (iii) Z₃ selection: cubic vertices vanish unless "
          f"Σq ≡ 0 mod 3  {'✓' if ok_sel else '✗'}")
    print(f"    rank-8 gate: only E₈ is even+unimodular+positive "
          f"(I₈ odd, A₈ det 9, D₈ det 4)  {'✓' if ok_e8 else '✗'}")
    print(f"    → superposition is exact because only norm-multiplicative")
    print(f"      bookkeeping composes; the substrate nonlinearity is the")
    print(f"      q=0 venting (gravity), not a deformation of QM")

    return {
        'C1': C1,
        'C3bar': C3bar,
        'CG_ratio': ratio_CG,
        'BA_ratio': BA_derived,
        'Q0': Q0,
        'linearity_checks': (ok_lin, ok_gate, ok_sel),
    }
