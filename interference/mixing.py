"""
Mixing: CKM + PMNS from (d₁₀, d₁₁, n₇, n₂₆, π).

The parameter map evaluates all four Wolfenstein coordinates from
(d₁₀, d₁₁, π).  The modular phase π/9 is computed; assigning it
to the physical Wolfenstein coordinate η̄ is an identification hypothesis,
as are the assignments of λ, A, and ρ̄.
CKM scorecard: diagonal pull score = 3.13, max marginal pull 0.70σ.
The comparison block mixes correlated global-fit outputs with separately
reported phase constraints and supplies no joint covariance, so this score
is descriptive and is not a chi-square statistic.
The map's −arg(U_e3) = 77.0 deg is convention-dependent and is not the
physical δ_CP.

═══════════════════════════════════════════════════════════════════════
CKM parameter map
═══════════════════════════════════════════════════════════════════════

The engine uses:

  λ = h₁₀·(1 + h₁₁/K²) = (2/9)(73/72) = 73/324 , METRIC face
      Sugawara weight h₁₀ = d₁₀/d₁₁² carrying the adjoint backreaction
      h₁₁/K² = 1/72, the same 73/72 root.py applies to bridge².
  A = √(n¹¹_off/d₁₁) = √(d₁₀/d₁₁) = √(2/3), ALGEBRAIC face
      The D⁽⁶⁾ adjoint nimrep derives a two-of-three off-diagonal
      fraction inside each boundary pair. It has zero matrix elements
      between the three pairs.

  η̄ = π/d₁₁² = π/9                            , MODULAR face (S-matrix)
      The Kac-Peterson S-matrix gives arg(S_ff) = π/d₁₁² exactly.

  ρ̄ = h₁₀|X| = h₁₀/√d₁₀ = √d₁₀/d₁₁² = √2/9, TOPOLOGICAL face
      The exact D⁽⁶⁾* cell system gives a Yang-Baxter connection
      with |X| = 1/√d₁₀.

  γ = arctan(η̄/ρ̄) = arctan(π/√d₁₀)           , downstream

  The apex is one complex number:
      apex = ρ̄ + iη̄ = (h₁₀/√d₁₀)·(1 + iπ/√d₁₀)
  γ and J follow from this apex.

═══════════════════════════════════════════════════════════════════════
PMNS, democratic mixing + charged-lepton backreaction
═══════════════════════════════════════════════════════════════════════

  Conjugation-orbit compression of the bulk fusion matrices gives a
  democratic adjoint block with TBM eigenvectors.

  U_PMNS = U_e† · U_TBM, with U_e = R₂₃(−φ) · R₁₃(θ_C, δ):
    φ = √(m_e/m_μ)/√d₁₀
    θ_C = h₁₀ = 2/9  (bare Sugawara weight)
    δ = arccos(−λ)    (echoed λ = 73/324)

Reference:
    Kahsay, Kibrom Kidane (2026). One-field repository.
    Papers and archived releases are listed in README.md.
"""

from fractions import Fraction
from functools import lru_cache
import cmath
import itertools
import math
import numpy as np

from .root import (d10, d11, h10, h11, K, Q0,
                   PDG_CKM, NUFIT_PMNS)


# ═══════════════════════════════════════════════════════════════════════
#  WZW S-matrix infrastructure (SU(3) at level k = d₁₁ = 3)
# ═══════════════════════════════════════════════════════════════════════

def _su3_inner(a, b):
    """SU(3) Killing-form inner product on weight space."""
    return (2*a[0]*b[0] + a[0]*b[1] + a[1]*b[0] + 2*a[1]*b[1]) / 3.0


_WEYL = [
    (lambda a, b: (a, b), +1),
    (lambda a, b: (-a, a+b), -1),
    (lambda a, b: (a+b, -b), -1),
    (lambda a, b: (-(a+b), a), +1),
    (lambda a, b: (b, -(a+b)), +1),
    (lambda a, b: (-b, -a), -1),
]


def _compute_wzw():
    """Compute the SU(3)₃ S-matrix and Verlinde fusion coefficients."""
    K_LEVEL = d11                         # = 3
    K_ALT   = K_LEVEL + d11              # altitude = 6 = d₁₀·d₁₁

    PRIMARIES = [(l1, l2) for l1 in range(K_LEVEL + 1)
                 for l2 in range(K_LEVEL + 1) if l1 + l2 <= K_LEVEL]
    N_PRIM = len(PRIMARIES)
    IDX = {p: i for i, p in enumerate(PRIMARIES)}

    # Kac-Peterson S-matrix
    S_raw = np.zeros((N_PRIM, N_PRIM), dtype=complex)
    for a, (l1a, l2a) in enumerate(PRIMARIES):
        pa = (l1a + 1, l2a + 1)
        for b, (l1b, l2b) in enumerate(PRIMARIES):
            pb = (l1b + 1, l2b + 1)
            val = sum(d * cmath.exp(-2j * math.pi * _su3_inner(w(*pa), pb) / K_ALT)
                      for w, d in _WEYL)
            S_raw[a, b] = val

    norm = math.sqrt(np.sum(np.abs(S_raw[0, :])**2))
    S_mat = S_raw / norm
    S_mat /= S_mat[0, 0] / abs(S_mat[0, 0])

    # Verlinde fusion
    N_fus = np.zeros((N_PRIM, N_PRIM, N_PRIM), dtype=int)
    for a in range(N_PRIM):
        for b in range(N_PRIM):
            for cc in range(N_PRIM):
                val = sum(S_mat[a, s] * S_mat[b, s] * S_mat[cc, s].conj() / S_mat[0, s]
                          for s in range(N_PRIM))
                N_fus[a, b, cc] = round(val.real)

    return {
        'PRIMARIES': PRIMARIES, 'N_PRIM': N_PRIM, 'IDX': IDX,
        'S_mat': S_mat, 'N_fus': N_fus,
    }


def _compute_zc_boundary_nimrep(wzw):
    """Reconstruct the actual two-boundary nimrep of the Z_C invariant.

    The exponents are the two self-conjugate primaries (0,0) and (1,1).
    Their two-character spectral representation, transformed by the
    Hadamard boundary basis, gives non-negative integer matrices and
    satisfies the full SU(3)_3 fusion algebra.
    """
    primaries = wzw['PRIMARIES']
    idx = wzw['IDX']
    s_mat = wzw['S_mat']
    n_fus = wzw['N_fus']
    exponents = [idx[p] for p in primaries if p == (p[1], p[0])]
    assert tuple(primaries[i] for i in exponents) == ((0, 0), (1, 1))

    boundary_basis = np.array([[1.0, 1.0], [1.0, -1.0]]) / math.sqrt(2.0)
    matrices = []
    for a in range(len(primaries)):
        eigenvalues = np.array([s_mat[a, m] / s_mat[0, m]
                                for m in exponents])
        raw = boundary_basis @ np.diag(eigenvalues) @ boundary_basis.T
        integer = np.rint(raw.real).astype(int)
        assert np.allclose(raw, integer, atol=1e-12)
        assert np.all(integer >= 0)
        matrices.append(integer)

    fusion_residual = 0.0
    for a in range(len(primaries)):
        for b in range(len(primaries)):
            rhs = sum(n_fus[a, b, c] * matrices[c]
                      for c in range(len(primaries)))
            fusion_residual = max(
                fusion_residual,
                float(np.linalg.norm(matrices[a] @ matrices[b] - rhs)))
    assert fusion_residual == 0.0
    return {
        'exponents': tuple(primaries[i] for i in exponents),
        'matrices': tuple(matrices),
        'fusion_residual': fusion_residual,
    }


@lru_cache(maxsize=1)
def _compute_d6_boundary():
    """Compute the six-boundary graph and its exact braid connection.

    The fundamental nimrep is the canonical graph used by the generation
    word lemma. Its adjoint action is fixed by 3 x 3bar = 1 + 8.  The
    connected components of that adjoint action are three boundary
    pairs, and every pair has one diagonal and two off-diagonal channels.
    There are no adjoint channels between distinct pairs.

    In cell-system notation this regular graph is D⁽⁶⁾*.  The two
    algebraic cell magnitudes are fixed by the type-I equations and the
    relative sign by the type-II equation.  Substitution in the Ocneanu
    formula gives a unitary Yang--Baxter connection whose entries all
    have magnitude 1/sqrt(2).  Thus 1/sqrt(2) is derived here from the
    solved connection, not assigned from graph valency.  This is the
    n=1 specialization of Evans--Pugh, Theorem 8.2 (arXiv:0906.4307).
    """
    from .words import N_FUND

    n10 = np.asarray(N_FUND, dtype=float)
    identity = np.eye(n10.shape[0])
    n11 = n10 @ n10.T - identity

    # Obtain the generation pairs from the graph itself, not from labels.
    unseen = set(range(n11.shape[0]))
    boundary_pairs = []
    while unseen:
        seed = min(unseen)
        component = tuple(int(i) for i in np.flatnonzero(n11[seed]))
        assert len(component) == d10 and set(component) <= unseen
        boundary_pairs.append(component)
        unseen.difference_update(component)
    assert len(boundary_pairs) == d11

    blocks = tuple(n11[np.ix_(pair, pair)] for pair in boundary_pairs)
    expected = np.array([[1.0, 2.0], [2.0, 1.0]])
    assert all(np.array_equal(block, expected) for block in blocks)
    total_channels = float(blocks[0][0].sum())
    pair_offdiagonal_channels = float(blocks[0][0, 1])
    pair_offdiagonal_fraction = pair_offdiagonal_channels / total_channels

    degrees = n10.sum(axis=1)
    assert np.all(degrees == d10)

    # D^(6)* cell system W+ at q=exp(i*pi/6).  Each boundary pair is
    # (i_1,i_2), (j_1,j_2), or (k_1,k_2).  Cyclic rotations have the
    # same cell.  W- is obtained by swapping the two internal labels.
    sqrt3 = math.sqrt(3.0)
    cell_hi = math.sqrt((sqrt3 + 1.0) / 2.0)
    cell_lo = math.sqrt((sqrt3 - 1.0) / 2.0)
    internal_label = {
        vertex: label
        for pair in boundary_pairs
        for label, vertex in enumerate(pair, start=1)
    }
    high_cells = {(1, 1, 1), (1, 2, 2), (2, 1, 2), (2, 2, 1)}

    def cell(x, y, z):
        labels = (internal_label[x], internal_label[y], internal_label[z])
        if labels == (2, 2, 2):
            return -cell_lo
        return cell_hi if labels in high_cells else cell_lo

    triangles = tuple(
        (x, y, z)
        for x, y, z in itertools.product(range(6), repeat=3)
        if n10[x, y] and n10[y, z] and n10[z, x]
    )
    assert len(triangles) == 24

    # Type I: each oriented edge frame has norm [2]_q = sqrt(3).
    type_i_residual = 0.0
    for x, y in zip(*np.nonzero(n10)):
        norm = sum(abs(cell(x, y, z))**2
                   for z in range(6) if (x, y, z) in triangles)
        type_i_residual = max(type_i_residual, abs(norm - sqrt3))
    assert type_i_residual < 1e-12

    # Hecke projector U and braided connection X on every two-step path
    # space with fixed endpoints.  The graph is simple and has PF weight
    # one, so U is the outer product of the two closing triangle cells.
    q = cmath.exp(1j * math.pi / 6.0)
    local_connections = {}
    connection_magnitudes = []
    connection_unitarity_residual = 0.0
    for x, z in itertools.product(range(6), repeat=2):
        mids = tuple(y for y in range(6) if n10[x, y] and n10[y, z])
        if not mids:
            continue
        values = np.array([cell(x, y, z) for y in mids])
        hecke = np.outer(values, values.conjugate())
        connection = (q**(2.0/3.0) * np.eye(len(mids))
                      - q**(-1.0/3.0) * hecke)
        local_connections[(x, z)] = (mids, connection)
        connection_magnitudes.extend(abs(connection).ravel())
        connection_unitarity_residual = max(
            connection_unitarity_residual,
            float(np.linalg.norm(connection @ connection.conjugate().T
                                 - np.eye(len(mids)))))
    connection_magnitude = 1.0 / math.sqrt(float(d10))
    assert np.allclose(connection_magnitudes, connection_magnitude,
                       atol=1e-12)
    assert connection_unitarity_residual < 1e-12

    def path_operator(paths, position):
        index = {path: i for i, path in enumerate(paths)}
        operator = np.zeros((len(paths), len(paths)), dtype=complex)
        for path in paths:
            mids, connection = local_connections[
                (path[position - 1], path[position + 1])]
            source = mids.index(path[position])
            for target, middle in enumerate(mids):
                moved = list(path)
                moved[position] = middle
                operator[index[tuple(moved)], index[path]] = \
                    connection[target, source]
        return operator

    # Yang--Baxter is the type-II cell equation in operator form.
    paths3 = tuple(
        path for path in itertools.product(range(6), repeat=4)
        if all(n10[path[i], path[i + 1]] for i in range(3))
    )
    braid1 = path_operator(paths3, 1)
    braid2 = path_operator(paths3, 2)
    braid_residual = float(np.linalg.norm(
        braid1 @ braid2 @ braid1 - braid2 @ braid1 @ braid2))
    braid_commutator = float(np.linalg.norm(
        braid1 @ braid2 - braid2 @ braid1))
    assert braid_residual < 1e-12 and braid_commutator > 1e-12

    # At web depth four a fixed-endpoint path block is eight-dimensional.
    # The central full twist splits it canonically as 3 + 3 + 2.
    start = boundary_pairs[0][0]
    endpoint = boundary_pairs[1][0]
    paths4 = tuple(
        path for path in itertools.product(range(6), repeat=5)
        if path[0] == start and path[-1] == endpoint
        and all(n10[path[i], path[i + 1]] for i in range(4))
    )
    assert len(paths4) == 8
    depth4_braids = tuple(path_operator(paths4, i) for i in (1, 2, 3))
    full_twist = np.linalg.matrix_power(
        depth4_braids[0] @ depth4_braids[1] @ depth4_braids[2], 4)
    twist_spectrum = np.linalg.eigvals(full_twist)
    roots = (1.0 + 0.0j,
             cmath.exp(2j * math.pi / 3.0),
             cmath.exp(-2j * math.pi / 3.0))
    depth4_carrier_dimensions = tuple(
        int(np.count_nonzero(abs(twist_spectrum - root) < 1e-9))
        for root in roots
    )
    assert sorted(depth4_carrier_dimensions) == [2, 3, 3]

    return {
        'n10': n10,
        'n11': n11,
        'boundary_pairs': tuple(boundary_pairs),
        'adjoint_blocks': blocks,
        'pair_offdiagonal_fraction': pair_offdiagonal_fraction,
        'cell_magnitudes': (cell_hi, cell_lo),
        'type_i_residual': type_i_residual,
        'connection_magnitude': connection_magnitude,
        'connection_unitarity_residual': connection_unitarity_residual,
        'braid_residual': braid_residual,
        'braid_commutator': braid_commutator,
        'depth4_carrier_dimensions': depth4_carrier_dimensions,
    }


# ═══════════════════════════════════════════════════════════════════════
#  CKM derivation
# ═══════════════════════════════════════════════════════════════════════

def lambda_echo_menu():
    """Price the committed λ backreaction over the engine's echo menu.

    The λ assignment carries one backreaction echo x/Kⁿ on the base
    h₁₀, the inserted constant x from {h₁₀, h₁₁, Q₀, 1} and the
    altitude Kⁿ from {K, K², K³}.  All twelve candidates h₁₀(1 + x/Kⁿ) are
    evaluated exactly and compared with the measured Wolfenstein λ
    here, once, in the open: the uniqueness of the committed entry is
    a computed count.  Comparison output only; no prediction consumes
    it, and tests/test_interference.py freezes the count.
    """
    lam_pdg, lam_sig = PDG_CKM['lambda']
    inserted = (('h10', h10), ('h11', h11), ('Q0', Q0),
                ('1', Fraction(1)))
    entries = []
    for x_name, x in inserted:
        for n in (1, 2, 3):
            cand = h10 * (1 + x / K**n)
            pull = (float(cand) - lam_pdg) / lam_sig
            entries.append(((x_name, n), cand, pull))
    hits = [e for e in entries if abs(e[2]) <= 1.0]
    committed = ('h11', 2)
    committed_pull = next(p for key, _, p in entries if key == committed)
    return {
        'menu': {key: (str(cand), pull) for key, cand, pull in entries},
        'n_candidates': len(entries),
        'n_within_1sigma': len(hits),
        'hit_keys': [key for key, _, _ in hits],
        'committed': committed,
        'committed_is_unique_hit':
            len(hits) == 1 and hits[0][0] == committed,
        'committed_pull': committed_pull,
        'next_best_abs_pull': min(abs(p) for key, _, p in entries
                                  if key != committed),
    }


def _derive_ckm(wzw):
    """Derive the full CKM matrix from the SU(3)₃ boundary data.

    Four Wolfenstein parameters, each a face of the D⁽⁶⁾* graph:
      λ , metric (Sugawara weight + adjoint echo)
      A , algebraic (nimrep within-pair off-diagonal fraction)
      η̄ , modular (S-matrix phase)
      ρ̄ , topological (Ocneanu braid-connection amplitude)
    """

    # # Wolfenstein parameters (all from d₁₀/d₁₁) #
    d6 = _compute_d6_boundary()

    # λ = h₁₀·(1 + h₁₁/K²): the fundamental Sugawara weight h₁₀ = 2/9,
    # carrying the engine's quadratic-altitude echo.  The
    # echo factor is 73/72, the SAME one root.py already applies to
    # bridge² (root.py: .echo(["WZW(adj/K²)"], h11/K², 2)): a PAIRED
    # insertion is Z₃-invariant already at second order, so K² is the
    # admissible altitude, and the adjoint is the stepper the down-type
    # channel closes through.  Over the engine's full echo
    # menu, {h₁₀, h₁₁, Q₀, 1} × {K, K², K³}, this is the only entry
    # within 1σ of the measured Wolfenstein λ: the enumeration is
    # executable (lambda_echo_menu above), reported in the returned
    # dict, and frozen by tests/test_interference.py.
    lam_exact = h10 * (1 + h11 / K**2)          # exact rational ledger
    assert lam_exact == Fraction(73, 324)
    lam = float(lam_exact)

    # A² is the normalized off-diagonal fraction of each adjoint block.
    A = math.sqrt(d6['pair_offdiagonal_fraction'])
    assert abs(A**2 - float(d10) / float(d11)) < 1e-15

    # The oriented fundamental/fundamental modular phase is pi/9.
    # Identifying this mathematical phase with physical etabar is an IH.
    fund = wzw['IDX'][(1, 0)]
    eta = cmath.phase(wzw['S_mat'][fund, fund])
    assert abs(eta - math.pi / float(d11**2)) < 1e-12

    # The solved D^(6)* braid connection has |X|=1/√d₁₀.
    rho = float(h10) * d6['connection_magnitude']

    # Package the two apex coordinates once for the downstream CKM map.
    apex = complex(rho, eta)
    assert abs(apex - (float(h10)/math.sqrt(d10))
               * complex(1.0, math.pi/math.sqrt(d10))) < 1e-15

    # γ = arg(apex) = arctan(π/√d₁₀), downstream.
    gamma_rad = math.atan2(eta, rho)
    gamma_deg = math.degrees(gamma_rad)

    # # Full CKM matrix #
    s12 = lam
    s23 = A * lam**2
    c12 = math.sqrt(1 - s12**2)
    c23 = math.sqrt(1 - s23**2)
    z13 = apex * s23 * s12 * c23 / (c12 * (1.0 - apex * s23**2))
    s13 = abs(z13)
    delta_rad = cmath.phase(z13)
    c13 = math.sqrt(1 - s13**2)

    ed = cmath.exp(1j * delta_rad)
    V = np.array([
        [c12*c13,                        s12*c13,                       s13*cmath.exp(-1j*delta_rad)],
        [-s12*c23 - c12*s23*s13*ed,      c12*c23 - s12*s23*s13*ed,     s23*c13],
        [s12*s23 - c12*c23*s13*ed,      -c12*s23 - s12*c23*s13*ed,     c23*c13],
    ], dtype=complex)
    Vabs = np.abs(V)

    # Jarlskog
    J_ckm = c12 * c23 * c13**2 * s12 * s23 * s13 * math.sin(delta_rad)

    # UT angles
    beta = math.atan2(eta, 1 - rho)
    alpha_ut = math.pi - gamma_rad - beta
    sin2beta = math.sin(2 * beta)

    alpha_deg = math.degrees(alpha_ut)
    beta_deg = math.degrees(beta)

    # # Correlated comparison coordinates and diagonal pull score #
    labels_row = ['u', 'c', 't']
    labels_col = ['d', 's', 'b']
    lam_pdg, lam_sig = PDG_CKM['lambda']
    A_pdg, A_sig = PDG_CKM['A']
    rho_pdg, rho_sig = PDG_CKM['rhobar']
    eta_pdg, eta_sig = PDG_CKM['etabar']
    gamma_pdg, gamma_sig = PDG_CKM['gamma_deg']
    J_pdg, J_sig = PDG_CKM['J']

    ckm_obs = [
        ('lambda', lam, lam_pdg, lam_sig),
        ('A',      A,   A_pdg,   A_sig),
        ('rhobar', rho, rho_pdg, rho_sig),
        ('etabar', eta, eta_pdg, eta_sig),
    ]
    for i in range(3):
        for j in range(3):
            key = f"V{labels_row[i]}{labels_col[j]}"
            pdg_val, pdg_err = PDG_CKM[key]
            ckm_obs.append((f"|V{labels_row[i]}{labels_col[j]}|", Vabs[i, j], pdg_val, pdg_err))
    ckm_obs.append(('gamma_deg', gamma_deg, gamma_pdg, gamma_sig))
    ckm_obs.append(('J', J_ckm, J_pdg, J_sig))

    diagonal_pull_score_ckm = sum(
        ((p - v) / e)**2 for _, p, v, e in ckm_obs)
    max_pull_ckm = max(abs((p - v) / e) for _, p, v, e in ckm_obs)

    return {
        'lam': lam, 'A': A, 'eta': eta, 'rho': rho,
        'V': V, 'Vabs': Vabs,
        'gamma_deg': gamma_deg, 'J_ckm': J_ckm,
        'alpha_deg': alpha_deg, 'beta_deg': beta_deg,
        'sin2beta': sin2beta,
        'd6_boundary_pairs': d6['boundary_pairs'],
        'd6_adjoint_blocks': d6['adjoint_blocks'],
        'd6_cell_magnitudes': d6['cell_magnitudes'],
        'd6_connection_magnitude': d6['connection_magnitude'],
        'd6_braid_residual': d6['braid_residual'],
        'd6_depth4_carrier_dimensions': d6['depth4_carrier_dimensions'],
        'ckm_obs': ckm_obs,
        'diagonal_pull_score_ckm': diagonal_pull_score_ckm,
        'max_pull_ckm': max_pull_ckm,
        'lambda_menu': lambda_echo_menu(),
    }


# ═══════════════════════════════════════════════════════════════════════
#  PMNS derivation
# ═══════════════════════════════════════════════════════════════════════

def _derive_pmns(wzw, lam_ckm, m_e, m_mu):
    """Derive the PMNS candidate from the SU(3)₃ modular data.

    The conjugation modular invariant Z_C of SU(3)₃ pairs each primary
    (l₁,l₂) with its conjugate (l₂,l₁). Compressing bulk fusion to
    conjugation-even orbit sums produces a democratic principal block:

    1. The adjoint principal block and a further projected fundamental
       block give TBM magnitudes up to ordering and rephasing. They are
       an algebraic candidate, not a Z_C nimrep or a neutrino mass operator.
    2. Charged-lepton rotation U_e provides backreaction (called "NLO
       corrections" in standard convention):
         φ = √(m_e/m_μ)/√d₁₀  (lepton mass ratio, derived in masses.py)
         θ_C = h₁₀ = 2/9       (conformal weight = Cabibbo angle)
         δ = arccos(−λ_CKM)    (phase consistency)
    3. U_PMNS = U_e† · U_TBM gives three mixing angles.  The directly
       reported −arg(U_e3) is convention-dependent, not physical δ_CP.
    """

    PRIMARIES = wzw['PRIMARIES']
    IDX = wzw['IDX']
    N_fus = wzw['N_fus']
    S_mat = wzw['S_mat']
    zc_nimrep = _compute_zc_boundary_nimrep(wzw)

    # # Conjugation invariant orbit structure #
    CONJ_MAP = {(l1, l2): (l2, l1) for l1, l2 in PRIMARIES}
    orbits = []; seen = set()
    for p in PRIMARIES:
        if p not in seen:
            cp = CONJ_MAP[p]
            orbits.append((p,) if cp == p else (p, cp))
            seen.update([p, cp])

    N_orb = len(orbits)

    # Normalized conjugation-orbit compression Q^T N_lambda Q.  This is
    # bulk fusion data, not the Z_C Cardy nimrep (which has Tr Z_C = 2
    # boundary states).  The diagnostics below keep that distinction
    # executable rather than burying it in prose.
    idx_11 = IDX[(1, 1)]
    idx_10 = IDX[(1, 0)]
    orbit_basis = np.zeros((len(PRIMARIES), N_orb))
    for j, orbit in enumerate(orbits):
        for primary in orbit:
            orbit_basis[IDX[primary], j] = 1.0 / math.sqrt(len(orbit))

    folded = tuple(orbit_basis.T @ N_fus[a] @ orbit_basis
                   for a in range(len(PRIMARIES)))
    b11_fold = folded[idx_11]
    b10_fold = folded[idx_10]
    # The three non-invertible conjugate-pair orbits form the generation
    # block. This replaces the former hand-written label filter; the
    # physical identification of this block with generations: the
    # non-invertible conjugate pairs are the generation carriers.
    qdim = np.real(S_mat[0, :] / S_mat[0, 0])
    gen = [i for i, o in enumerate(orbits)
           if len(o) == 2 and qdim[IDX[o[0]]] > 1.0 + 1e-9]
    assert len(gen) == d11
    blk = b11_fold[np.ix_(gen, gen)]
    assert blk.shape == (3, 3) and np.allclose(blk, blk[0, 0]), \
        "conjugation-orbit adjoint principal block must be democratic"
    ev = sorted(np.linalg.eigvalsh(blk), reverse=True)
    assert abs(ev[0] - 3*blk[0, 0]) < 1e-9 and abs(ev[1]) < 1e-9

    # The projected fundamental principal block resolves that plane. Its eigenvalues
    # {-1/2, 1/2, 1} are non-degenerate, and the absolute eigenvector
    # matrix is TBM up to the unphysical ordering and rephasing of the
    # boundary states. Canonicalize those freedoms to the convention used
    # by the charged-lepton rotation below.
    complement = [i for i in range(N_orb) if i not in gen]
    fundamental_generation_leakage = float(np.linalg.norm(
        b10_fold[np.ix_(gen, complement)]))
    assert fundamental_generation_leakage > 1e-12

    # A representation would obey B_a B_b = sum_c N_ab^c B_c. Orbit
    # compression fails this test because non-self-conjugate fusion does
    # not preserve the conjugation-even subspace.
    fold_fusion_residual = 0.0
    for a in range(len(PRIMARIES)):
        for b in range(len(PRIMARIES)):
            rhs = sum(N_fus[a, b, c] * folded[c]
                      for c in range(len(PRIMARIES)))
            fold_fusion_residual = max(
                fold_fusion_residual,
                float(np.linalg.norm(folded[a] @ folded[b] - rhs)))
    assert fold_fusion_residual > 1e-12

    blk10 = b10_fold[np.ix_(gen, gen)]
    vals10, vecs10 = np.linalg.eigh(blk10)
    targets = (-0.5, 1.0, 0.5)  # mass-column convention: 1, 2, 3
    col_order = [int(np.argmin(abs(vals10 - target))) for target in targets]
    assert len(set(col_order)) == d11
    U_TBM = vecs10[:, col_order]
    row_order = np.argsort(-np.abs(U_TBM[:, 0]), kind='stable')
    U_TBM = U_TBM[row_order, :]
    if U_TBM[0, 0] < 0:
        U_TBM[:, 0] *= -1
    if U_TBM[2, 0] < 0:
        U_TBM[2, :] *= -1
    if U_TBM[0, 1] < 0:
        U_TBM[:, 1] *= -1
    if U_TBM[1, 2] < 0:
        U_TBM[:, 2] *= -1
    U_TBM = U_TBM.astype(complex)
    tbm_abs = np.array([
        [2/math.sqrt(6), 1/math.sqrt(3), 0],
        [1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
        [1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
    ])
    assert np.allclose(np.abs(U_TBM), tbm_abs)

    # # Charged-lepton backreaction: U_e = R₂₃(−φ) · R₁₃(θ_C, δ) #
    # (Standard convention: "NLO correction to TBM.")
    # The rotation form uses the CKM-derived coordinates.
    # φ = √(m_e/m_μ)/√d₁₀; m_e and m_μ come from masses.py.
    phi = math.sqrt(m_e / m_mu) / math.sqrt(float(d10))
    # θ_C = h₁₀ BARE. arctan(λ) forks from h₁₀ under the echoed λ;
    # both readings are measured in the header. The primary one is bare.
    theta = float(h10)
    # δ = arccos(−λ_CKM).
    delta = math.acos(-lam_ckm)

    def R13(a, d):
        c, s, e = math.cos(a), math.sin(a), cmath.exp(1j*d)
        return np.array([[c,0,s*e],[0,1,0],[-s*e.conjugate(),0,c]], dtype=complex)

    def R23(a):
        c, s = math.cos(a), math.sin(a)
        return np.array([[1,0,0],[0,c,s],[0,-s,c]], dtype=complex)

    U_e = R23(-phi) @ R13(theta, delta)
    U_PMNS = np.conj(U_e.T) @ U_TBM

    # Extract mixing angles
    s13 = abs(U_PMNS[0, 2])
    c13sq = 1 - s13**2
    sin2_13 = s13**2
    sin2_12 = abs(U_PMNS[0, 1])**2 / c13sq
    sin2_23 = abs(U_PMNS[1, 2])**2 / c13sq
    ue3_phase_deg = -math.degrees(cmath.phase(U_PMNS[0, 2]))

    J_lep = abs((U_PMNS[0,0]*U_PMNS[1,1]*np.conj(U_PMNS[0,1])*np.conj(U_PMNS[1,0])).imag)

    # Marginal pulls vs the frozen NuFit 6.0 snapshot. The published
    # intervals are asymmetric and the angles are correlated, so their
    # sum of squared symmetric pulls is descriptive, not a chi-square.
    pull_12 = (sin2_12 - NUFIT_PMNS['sin2_12'][0]) / NUFIT_PMNS['sin2_12'][1]
    pull_23 = (sin2_23 - NUFIT_PMNS['sin2_23'][0]) / NUFIT_PMNS['sin2_23'][1]
    pull_13 = (sin2_13 - NUFIT_PMNS['sin2_13'][0]) / NUFIT_PMNS['sin2_13'][1]
    diagonal_pull_score_pmns = pull_12**2 + pull_23**2 + pull_13**2

    return {
        'sin2_12': sin2_12, 'sin2_23': sin2_23, 'sin2_13': sin2_13,
        'ue3_phase_deg': ue3_phase_deg,
        'J_lep': J_lep,
        'diagonal_pull_score_pmns': diagonal_pull_score_pmns,
        'pulls': {'12': pull_12, '23': pull_23, '13': pull_13},
        'generation_orbits': tuple(orbits[i] for i in gen),
        'zc_boundary_states': len(zc_nimrep['exponents']),
        'zc_n10': zc_nimrep['matrices'][IDX[(1, 0)]],
        'zc_n11': zc_nimrep['matrices'][IDX[(1, 1)]],
        'zc_fusion_residual': zc_nimrep['fusion_residual'],
        'fundamental_generation_leakage': fundamental_generation_leakage,
        'fold_fusion_residual': fold_fusion_residual,
        'democratic_eigenvalues': tuple(float(x) for x in ev),
        'fundamental_eigenvalues': targets,
        'tbm_orbit_order': tuple(orbits[gen[int(i)]] for i in row_order),
        'U_TBM': U_TBM,
        'U_PMNS': U_PMNS,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive(R, masses):
    """
    Compute CKM + PMNS mixing matrices.

    R      : dict from root.derive()
    masses : dict from masses.derive()
    """
    from fractions import Fraction

    print("\n" + "=" * 78)
    print("  MIXING: CKM + PMNS from (d₁₀, d₁₁) = (2, 3)")
    print("=" * 78)

    # # WZW S-matrix #
    print(f"\n  SU(3) at level k = d₁₁ = {d11}, altitude K = d₁₀·d₁₁ = {d10*d11}")
    wzw = _compute_wzw()
    S_mat = wzw['S_mat']
    N_PRIM = wzw['N_PRIM']
    U_check = S_mat @ S_mat.conj().T
    unitarity_err = np.max(np.abs(U_check - np.eye(N_PRIM)))
    print(f"  S-matrix: {N_PRIM}x{N_PRIM}, unitarity |SS†−I| = {unitarity_err:.2e}")

    # # CKM #
    print(f"\n  CKM candidate map (physical coordinate assignments are IHs):")
    ckm = _derive_ckm(wzw)

    h10_frac = Fraction(d10, d11**2)
    print(f"    lambda [IH] = h₁₀(1+h₁₁/K²) = ({h10_frac})(73/72) = "
          f"{Fraction(73, 324)} = {ckm['lam']:.5f}"
          f"  (PDG: {PDG_CKM['lambda'][0]}, pull: {(ckm['lam']-PDG_CKM['lambda'][0])/PDG_CKM['lambda'][1]:+.2f}sigma)")
    _menu = ckm['lambda_menu']
    print(f"    lambda echo menu {{h₁₀,h₁₁,Q₀,1}}x{{K,K²,K³}}: "
          f"{_menu['n_within_1sigma']}/{_menu['n_candidates']} within "
          f"1sigma (committed {_menu['committed_pull']:+.2f}sigma, "
          f"next best {_menu['next_best_abs_pull']:.1f}sigma)")
    print(f"    A [IH] = sqrt(d₁₀/d₁₁) = sqrt({Fraction(d10,d11)}) = {ckm['A']:.5f}"
          f"  (PDG: {PDG_CKM['A'][0]}, pull: {(ckm['A']-PDG_CKM['A'][0])/PDG_CKM['A'][1]:+.2f}sigma)")
    print(f"    modular phase pi/d₁₁² = pi/{d11**2} -> etabar [IH] "
          f"= {ckm['eta']:.5f}"
          f"  (PDG: {PDG_CKM['etabar'][0]}, pull: {(ckm['eta']-PDG_CKM['etabar'][0])/PDG_CKM['etabar'][1]:+.2f}sigma)")
    print(f"    rhobar [IH] = h₁₀/sqrt(d₁₀) = sqrt({d10})/{d11**2} = {ckm['rho']:.5f}"
          f"  (PDG: {PDG_CKM['rhobar'][0]}, pull: {(ckm['rho']-PDG_CKM['rhobar'][0])/PDG_CKM['rhobar'][1]:+.2f}sigma)")
    print(f"    gamma [IH-map output] = arctan(pi/sqrt({d10})) = {ckm['gamma_deg']:.2f} deg"
          f"  (PDG: {PDG_CKM['gamma_deg'][0]})")

    labels_row = ['u', 'c', 't']
    labels_col = ['d', 's', 'b']
    Vabs = ckm['Vabs']
    for i in range(3):
        row = "    "
        for j in range(3):
            key = f"V{labels_row[i]}{labels_col[j]}"
            pdg_val, pdg_err = PDG_CKM[key]
            pull = (Vabs[i, j] - pdg_val) / pdg_err
            row += f"|V_{labels_row[i]}{labels_col[j]}|={Vabs[i,j]:.5f}({pull:+.1f}s) "
        print(row)

    print(f"    J_ckm = {ckm['J_ckm']:.2e}  (PDG: {PDG_CKM['J'][0]:.2e})")
    print(f"    diagonal pull score (CKM) = "
          f"{ckm['diagonal_pull_score_ckm']:.2f} / "
          f"{len(ckm['ckm_obs'])} comparison coordinates = "
          f"{ckm['diagonal_pull_score_ckm']/len(ckm['ckm_obs']):.2f}")
    print("      global-fit outputs plus a separate phase constraint; no joint")
    print("      covariance model is supplied, so this is descriptive, not chi2")

    # # PMNS #
    m_e = masses['m_e']
    m_mu = masses['m_mu']
    pmns = _derive_pmns(wzw, ckm['lam'], m_e, m_mu)

    print(f"\n  PMNS candidate (conjugation-orbit projection + charged-lepton backreaction):")
    print(f"    Z_C Cardy boundaries = {pmns['zc_boundary_states']}; "
          f"fundamental projection leakage = "
          f"{pmns['fundamental_generation_leakage']:.3f}")
    print(f"    phi = sqrt(m_e/m_mu)/sqrt({d10}) = {math.sqrt(m_e/m_mu)/math.sqrt(float(d10)):.6f} rad")
    print(f"    theta_C = h₁₀ = {h10_frac} = {math.degrees(float(h10)):.3f} deg")
    print(f"    delta = arccos(-lambda) = {math.degrees(math.acos(-ckm['lam'])):.3f} deg")

    for label, key in [('12', 'sin2_12'), ('23', 'sin2_23'), ('13', 'sin2_13')]:
        val = pmns[key]
        nf_val, nf_err = NUFIT_PMNS[f'sin2_{label}']
        pull = (val - nf_val) / nf_err
        print(f"    sin2_theta_{label} = {val:.5f}  (NuFit: {nf_val}, pull: {pull:+.2f}sigma)")

    print(f"    -arg(U_e3) = {pmns['ue3_phase_deg']:.1f} deg  "
          "(current phase convention; not physical delta_CP)")
    print(f"    J_lep = {pmns['J_lep']:.6f}")
    print(f"    diagonal pull score (PMNS) = "
          f"{pmns['diagonal_pull_score_pmns']:.2f} / 3 correlated angles = "
          f"{pmns['diagonal_pull_score_pmns']/3:.2f}")
    print("      approximate symmetric marginals; not a covariance-aware chi2")

    # # Merge and return #
    result = {}
    result.update(ckm)
    result.update(pmns)
    return result
