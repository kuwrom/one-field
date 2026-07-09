"""
circulant_normalization.py: three circulants, one slot, and which of
them can be the mass matrix.

THE SLOT.  Lepton masses come from a Z₃-equivariant insertion, which is
a real circulant Δ = c₀I + c₁S + c₁*S†.  Its eigenvalues are

    Δ_k = A + |B| cos(θ + 2πk/3),   A = c₀,  |B| = 2|c₁|,  θ = arg c₁

and on the positive-gap branch (all Δ_k > 0) the framework reads
√m_k = Δ_k.  Three different objects in this repository claim that
slot.  This probe puts them side by side and computes which can be it.

  (1) CANON        A = 1, |B| = √d₁₀ = √2, θ = h₁₀ = 2/9
                   root.py:574, masses.py:277, octonions.py:329
  (2) HARMONIC     v = (1, −½, −½) at ‖v‖² = 3/2, i.e. (B, θ) = (1, 0)
                   the retired normalisation argument in
                   protected_forgetting.py
  (3) BdG UNIFORM  the (B/A, θ) of the substrate's own circulant
                   (nls_soliton.interaction_matrix) linearised about
                   the uniform condensate; at α = 0,
                   B/A = 2(1 + 3r/2)/(1 − r) with r = g₁/g₀

WHAT IS PROVED HERE (computed, not cited):

  PART A: the master identity.  Q ≡ Σm_k/(Σ√m_k)² = 1/3 + (B/A)²/6,
      exactly and for every θ.  Hence Q is blind to the phase, and
      Q = 2/3 pins the AMPLITUDE alone: B/A = √2.  A normalisation
      condition can therefore never fix θ.

  PART B: the moment bookkeeping.  Σv_k²/3 = B²/2 (total, both Z₃
      channels) while |v̂₁|² = B²/4 (per channel), so Pv²P = 2·E[v²].
      At the canonical B = √2 the total is 1 and the per-mode is ½.
      The web's `E_v2` is the per-mode power; the projector notation
      "Pv²P = ½P" names the total.  The two are distinct, and
      conflating them is the error this part rules out.

  PART C: the θ = 0 degeneracy theorem.  For θ = 0 the spectrum is
      (1+B, 1−B/2, 1−B/2) for EVERY B: two masses exactly equal.  Any
      circulant with a real off-diagonal is θ = 0 and therefore cannot
      carry three generations, at any coupling.

  PART D: the branch-exclusivity theorem for the substrate.  For the
      BdG uniform circulant of the hermitian equation at α = 0,
          B/A = 2(1 + 3r/2)/(1 − r),   r = g₁/g₀,
      and B/A = √2 requires r = −0.1327 < 0, i.e. no modulational
      instability (max growth 0).  The positive-gap branch (which the
      Koide reading requires) and the unstable branch (which makes
      knots, i.e. matter) are disjoint.

  PART E: the positivity window.  At the Koide-forced B = √2, the
      requirement that all three gaps stay positive admits exactly a
      quarter of θ ∈ [0, 2π).  θ = 2/9 lies inside it with a thin
      margin: the electron sits near the boundary, and the lepton
      hierarchy is that thinness.

STATUS:

  · PARTS F-G construct and dynamically verify a stable UNIFORM state
    whose effective circulant has (B/A, θ) = (√2, 2/9).
    knot_circulant.py reads the anomalous
    channel on a certified stationary knot (|sin β| = 0.70 at the
    core, against exactly 0 for every uniform state) and lifts the
    q1/q2 degeneracy at the eigenbasis level.
    So the knot escapes the uniform-background obstruction
    QUALITATIVELY.  What the internal spectrum of an isolated knot
    carries is not measured here: no probe in this repository
    constructs one reproducibly.  The interpretation and its
    falsifier are LAYER SCOPING below.
  · The substrate runs at α = 0, where θ = 0 identically; the blind
    clock measures Q ≈ 1/3 (emergent_clock.py).  A complex Z₃ phase
    α ≠ 0 is hermitian-admissible (nls_soliton.interaction_matrix)
    and reaches (√2, 2/9) at exactly two inverse images (r, α), both
    of which are modulationally STABLE (full 6x6 BdG, PART F;
    confirmed by direct dynamics, PART G).  The canonical point is a
    stable uniform phase carrying the Koide circulant; knots cannot
    form there.  The effective-coupling flow and confined-phase route
    have no implementation.  A dynamical
    origin for B/A = sqrt(2) via instanton tunnelling was tested and
    REJECTED by rigorous computation.

LAYER SCOPING, stated with its falsifier.  The classical field is the
first folding (amplitude and phase) and would provide the SLOTS: a
stable loop, three internal modes, channel-lifting.  The RATIOS are
written at the deeper foldings: theta = 2/9 is a Sugawara conformal
weight, a datum of the quantized bookkeeping, and the Brannen
hierarchy is a near-cancellation (one gap almost closing), an O(1)
structural fact, not an alpha/2pi response.  The classical layer is
structurally blind to it: Q = 1/3 classically is the exact image of
the spectrum at theta = 0.  Q in [1/3, 1) holds identically for any
three positive lines, so the classical readout confirms nothing on
its own; tests/test_emergence.py reports Q and declines to assert it.
FALSIFIER: a classical measurement of Q = 2/3 would falsify the layer
scoping.

Usage: python3 circulant_normalization.py
"""

import cmath
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
import _bootstrap  # noqa: E402,F401
import nls_soliton as _nls  # noqa: E402

B_CANON = math.sqrt(2.0)
THETA_CANON = 2.0 / 9.0


def spectrum(A, B, theta):
    """Δ_k = A + B cos(θ + 2πk/3)."""
    return np.array([A + B * math.cos(theta + 2 * math.pi * k / 3)
                     for k in range(3)])


def koide_Q(delta):
    """Q = Σm/(Σ√m)² with √m = Δ.  None if off the positive-gap branch."""
    if (delta <= 0).any():
        return None
    return float((delta ** 2).sum() / delta.sum() ** 2)


def moments(B, theta):
    """(total Σv²/3, per-mode |v̂₁|²) for v_k = B cos(θ + 2πk/3)."""
    v = spectrum(0.0, B, theta)
    total = float((v ** 2).sum() / 3.0)
    omega = np.exp(2j * np.pi / 3)
    v_hat1 = sum(v[k] * omega ** (-k) for k in range(3)) / 3.0
    return total, float(abs(v_hat1) ** 2)


def bdg_circulant(g0, g1, rho0=1.0, alpha=0.0):
    """The uniform-condensate normal self-energy circulant of the
    hermitian substrate, read in the (A, |B|, θ) convention.

        c₀ = (g₀ − g₁cos α) ρ₀
        c₁ = [g₀ + g₁cos α + (g₁/2) e^{−iα}] ρ₀
        A = c₀,   |B| = 2|c₁|,   θ = arg c₁

    At α = 0 this reduces to A = (g₀−g₁)ρ₀, |B| = 2(g₀+3g₁/2)ρ₀, θ = 0,
    i.e. B/A = 2(1 + 3r/2)/(1 − r) with r = g₁/g₀, not the 2 + r of a
    non-hermitian equation (see nls_soliton.py).  Verified against a
    finite-difference Wirtinger Jacobian of nls_soliton's own
    interaction_matrix in `test_matches_live_substrate` below.

    Returns (A, |B|, theta, eigenvalues-of-the-normal-circulant).
    """
    ca = math.cos(alpha)
    c0 = (g0 - g1 * ca) * rho0
    c1 = (g0 + g1 * ca + (g1 / 2.0) * cmath.exp(-1j * alpha)) * rho0
    theta = math.atan2(c1.imag, c1.real)
    B = 2.0 * abs(c1)
    eigs = np.array([c0 + B * math.cos(theta + 2 * math.pi * k / 3)
                     for k in range(3)])
    return c0, B, theta, eigs


def bdg_max_growth(g0, g1, rho0=1.0, alpha=0.0, n=1200, p_max=4.0):
    """Maximum modulational growth rate over p from the FULL 6x6 BdG
    L = [[eps I + A, B], [-B*, -(eps I + A)*]], no channel-pairing
    assumption.  A and B are finite-difference Wirtinger derivatives
    of nls_soliton's own interaction matrix, so this tracks the live
    code rather than a transcription.  Returns max Im omega >= 0.

    CORRECTION: the per-channel
    form (eps + a_q)^2 - |b_q|^2 pairs a_q with b_q of the SAME q,
    which is valid only at alpha = 0 where a_1 = a_2.  At alpha != 0
    the anomalous term couples channel q with channel -q, and the
    per-channel form reported spurious instability at the canonical
    point.  The full 6x6 reproduces the alpha = 0 rates exactly
    (g1 rho0 at both r values) and is confirmed by direct dynamics
    (PART G): the seeded 'unstable' mode does not grow.
    """
    s = math.sqrt(rho0)
    p0 = np.full(3, s, dtype=complex)

    def F(p):
        M = _nls.interaction_matrix(p.reshape(3, 1), g0, g1, alpha)[0]
        return M @ p

    mu = (F(p0)[0] / p0[0]).real
    h = 1e-7
    A = np.zeros((3, 3), dtype=complex)
    Bm = np.zeros((3, 3), dtype=complex)
    for j in range(3):
        e = np.zeros(3, dtype=complex)
        e[j] = 1
        dr = (F(p0 + h * e) - F(p0 - h * e)) / (2 * h)
        di = (F(p0 + 1j * h * e) - F(p0 - 1j * h * e)) / (2 * h)
        A[:, j] = 0.5 * (dr - 1j * di)
        Bm[:, j] = 0.5 * (dr + 1j * di)
    A = A - mu * np.eye(3)
    worst = 0.0
    for p in np.linspace(1e-3, p_max, n):
        H0 = (p * p / 2.0) * np.eye(3) + A
        L = np.block([[H0, Bm], [-Bm.conj(), -H0.conj()]])
        worst = max(worst, float(np.linalg.eigvals(L).imag.max()))
    return worst


def run(report=print):
    out = {}

    report("=" * 78)
    report("  CIRCULANT NORMALIZATION: three claimants, one slot")
    report("=" * 78)

    ## PART A: the master identity ##
    report("\n  PART A: Q = 1/3 + (B/A)²/6, exactly, for every θ")
    report(f"    {'B/A':>8} {'θ':>8} {'Q measured':>16} "
           f"{'Q identity':>16} {'branch':>10}")
    max_dev = 0.0
    for BA in (0.0, 0.5, B_CANON, 1.9):
        for th in (0.0, 2.0 / 9.0, 0.7, 2.2):
            d = spectrum(1.0, BA, th)
            Q = koide_Q(d)
            Q_id = 1.0 / 3.0 + BA ** 2 / 6.0
            tag = "positive" if Q is not None else "OFF-BRANCH"
            if Q is not None:
                max_dev = max(max_dev, abs(Q - Q_id))
                report(f"    {BA:8.4f} {th:8.4f} {Q:16.12f} "
                       f"{Q_id:16.12f} {tag:>10}")
            else:
                report(f"    {BA:8.4f} {th:8.4f} {'--':>16} "
                       f"{Q_id:16.12f} {tag:>10}")
    report(f"    max deviation = {max_dev:.2e}   → Q is θ-BLIND")
    report(f"    therefore  Q = 2/3  ⟺  B/A = √2 = {B_CANON:.9f}, "
           f"and θ is NOT fixed by Q")
    assert max_dev < 1e-12
    out["identity_max_dev"] = max_dev

    ## PART B: moment bookkeeping ##
    report("\n  PART B: total vs per-mode second moment")
    report(f"    {'B':>10} {'Σv²/3 (total)':>16} {'|v̂₁|² (mode)':>16} "
           f"{'ratio':>8}")
    for B in (1.0, B_CANON, 1.7):
        tot, mode = moments(B, THETA_CANON)
        report(f"    {B:10.6f} {tot:16.9f} {mode:16.9f} "
               f"{tot / mode:8.4f}")
        assert abs(tot - 2.0 * mode) < 1e-12
    tot_c, mode_c = moments(B_CANON, THETA_CANON)
    report(f"    at B = √2:  total = {tot_c:.9f} (= 1),  "
           f"per-mode = {mode_c:.9f} (= ½)")
    report(f"    the web's E_v2 is the PER-MODE power; "
           f"'Pv²P = ½P' named the total")
    assert abs(tot_c - 1.0) < 1e-12 and abs(mode_c - 0.5) < 1e-12
    out["Pv2P_total"] = tot_c
    out["E_v2_mode"] = mode_c

    ## PART C: the θ = 0 degeneracy theorem ##
    report("\n  PART C: θ = 0 ⇒ exact two-fold degeneracy, for every B")
    report(f"    {'B':>10} {'Δ₀':>12} {'Δ₁':>12} {'Δ₂':>12} "
           f"{'|Δ₁−Δ₂|':>10}")
    theta0_deg = True
    for B in (0.2, 1.0, B_CANON, 1.9):
        d = spectrum(1.0, B, 0.0)
        split = abs(d[1] - d[2])
        theta0_deg &= split < 1e-12
        report(f"    {B:10.6f} {d[0]:12.6f} {d[1]:12.6f} {d[2]:12.6f} "
               f"{split:10.2e}")
    assert theta0_deg
    report("    → any circulant with a REAL off-diagonal (θ = 0) is "
           "structurally")
    report("      incapable of three generations, at any coupling")
    out["theta0_degenerate"] = bool(theta0_deg)

    ## the three claimants side by side ##
    report("\n  THE THREE CLAIMANTS")
    report(f"    {'name':>14} {'A':>8} {'B':>8} {'B/A':>8} {'θ':>7} "
           f"{'Q':>12} {'branch':>11}")

    d_canon = spectrum(1.0, B_CANON, THETA_CANON)
    Q_canon = koide_Q(d_canon)
    report(f"    {'CANON':>14} {1.0:8.4f} {B_CANON:8.4f} "
           f"{B_CANON:8.4f} {THETA_CANON:7.4f} {Q_canon:12.9f} "
           f"{'positive':>11}")

    d_harm = spectrum(1.0, 1.0, 0.0)
    Q_harm = koide_Q(d_harm)
    report(f"    {'HARMONIC':>14} {1.0:8.4f} {1.0:8.4f} {1.0:8.4f} "
           f"{0.0:7.4f} {Q_harm:12.9f} {'degenerate':>11}")

    g1_over_g0 = 1.0 / math.sqrt(2.0)
    A_b, B_b, th_b, eigs_b = bdg_circulant(1.0, g1_over_g0)
    Q_b = koide_Q(np.array(sorted(eigs_b)))
    report(f"    {'BdG UNIFORM':>14} {A_b:8.4f} {B_b:8.4f} "
           f"{B_b / A_b:8.4f} {th_b:7.4f} "
           f"{'--':>12} {'UNSTABLE':>11}")
    report(f"    normal-circulant eigenvalues = {np.round(eigs_b, 6)}  "
           f"→ off the positive-gap branch")
    assert Q_canon is not None and abs(Q_canon - 2 / 3) < 1e-12
    assert abs(Q_harm - 0.5) < 1e-12
    assert Q_b is None
    out["Q_canon"] = Q_canon
    out["Q_harmonic"] = Q_harm
    out["BdG_BA"] = B_b / A_b

    # # PART D: branch exclusivity for the substrate ##
    report("\n  PART D: the substrate's branch exclusivity (hermitian "
           "dynamics)")
    report("    At α = 0 the normal circulant of the hermitian equation "
           "gives")
    report("        B/A = 2(1 + 3r/2)/(1 − r),   r = g₁/g₀")
    report("    (a non-hermitian equation reads 2 + r instead; see "
           "nls_soliton.py)")
    report(f"    {'g₁/g₀':>10} {'B/A':>10} {'max growth':>12} "
           f"{'gaps>0':>8} {'unstable':>9}")
    for r in (0.0, 0.3, g1_over_g0, -0.1327, -0.5, -1.0):
        A_r, B_r, _, eigs_r = bdg_circulant(1.0, r)
        g_r = bdg_max_growth(1.0, r)
        report(f"    {r:10.4f} {B_r / A_r:10.6f} {g_r:12.6f} "
               f"{str(bool((eigs_r > 0).all())):>8} "
               f"{str(bool(g_r > 1e-8)):>9}")
        # |B| = 2|c₁| is a magnitude, so compare magnitudes: when
        # c₁ changes sign the phase jumps to θ = π rather than B/A
        # going negative.
        assert abs((B_r / A_r)
                   - abs(2.0 * (1 + 1.5 * r) / (1 - r))) < 1e-12
    r_needed = (B_CANON - 2.0) / (3.0 + B_CANON)
    A_n, B_n, _, _ = bdg_circulant(1.0, r_needed)
    g_n = bdg_max_growth(1.0, r_needed)
    report(f"    B/A = √2 requires g₁/g₀ = {r_needed:.6f} < 0 "
           f"(check: B/A = {B_n / A_n:.9f})")
    report(f"    at that coupling max growth = {g_n:.6f} → "
           f"NO instability, no knots")
    report("    → the Koide positive-gap branch and the knot-forming "
           "branch are")
    report("      DISJOINT at α = 0.  Not a tuning problem.")
    out["r_needed_for_sqrt2"] = r_needed
    assert r_needed < 0
    assert abs(B_n / A_n - B_CANON) < 1e-9
    assert g_n < 1e-8
    # validation of the 6x6 form: it reproduces the alpha = 0 theory
    # rate max growth = |b_1| = g1 rho0 at the knot-forming coupling
    g_check = bdg_max_growth(1.0, 1.0 / math.sqrt(2.0))
    assert abs(g_check - 1.0 / math.sqrt(2.0)) < 1e-3

    ## PART E: positivity window at the forced amplitude ##
    report("\n  PART E: positivity window at B/A = √2")
    n_scan = 36000
    n_ok = sum(1 for j in range(n_scan)
               if (spectrum(1.0, B_CANON,
                            2 * math.pi * j / n_scan) > 0).all())
    frac = n_ok / n_scan
    min_cos = min(math.cos(THETA_CANON + 2 * math.pi * k / 3)
                  for k in range(3))
    report(f"    admissible fraction of θ ∈ [0, 2π): {frac:.6f}")
    report(f"    boundary: min_k cos(θ + 2πk/3) > −1/√2 = "
           f"{-1 / math.sqrt(2):.6f}")
    report(f"    at θ = 2/9: min_k cos = {min_cos:.6f}  "
           f"(margin {min_cos + 1 / math.sqrt(2):+.6f})")
    report(f"    smallest gap Δ_min = {d_canon.min():.6f}  → the "
           f"electron sits")
    report(f"    near the positivity boundary")
    assert d_canon.min() > 0
    out["positivity_fraction"] = frac
    out["theta_margin"] = min_cos + 1 / math.sqrt(2)
    out["gap_min"] = float(d_canon.min())

    # # PART F: the complex Z₃ phase and what it does and doesn't buy #
    report("\n  PART F: a complex Z₃ coupling (α ≠ 0)")
    report("    The real projection that forces θ = 0 is NOT required by "
           "hermiticity.")
    report("    E_α = ∫{½Σ|∇ψ|² + (g₀/2)ρ_tot² + g₁Σρ_k Re(e^{iα}ψ̄ψ)} is "
           "manifestly")
    report("    real for every α, so its gradient is hermitian; α is the "
           "synthetic")
    report("    flux through the Z₃ triangle, gauge-invariant mod 2π/3, "
           "and α ≠ 0")
    report("    breaks time reversal.  Solving the two conditions "
           "B/A = √2, θ = 2/9:")
    r_c, a_c = 0.383029723, 4.333441893
    A_c, B_c, th_c, d_c = bdg_circulant(1.0, r_c, alpha=a_c)
    Q_c = koide_Q(d_c)
    g_c = bdg_max_growth(1.0, r_c, alpha=a_c)
    # the "second inverse image" is the SAME theory: (g1, alpha) ->
    # (-g1, alpha - pi) maps the interaction term g1 rho_k
    # Re(e^{i alpha} psibar psi) to itself identically, so computing it
    # exercises an exact redundancy, not a second point (asserted below)
    g_c2 = bdg_max_growth(1.0, -r_c, alpha=a_c - math.pi)
    redundancy_dev = 0.0
    for r_probe, a_probe in ((r_c, a_c), (0.61, 1.3), (-0.27, 5.1)):
        one = bdg_circulant(1.0, r_probe, alpha=a_probe)
        two = bdg_circulant(1.0, -r_probe, alpha=a_probe - math.pi)
        redundancy_dev = max(
            redundancy_dev,
            abs(one[0] - two[0]), abs(one[1] - two[1]),
            abs(one[2] - two[2]),
            float(np.max(np.abs(np.asarray(one[3]) - np.asarray(two[3])))))
    ms = np.sort(d_c ** 2)
    report(f"      g₁/g₀ = {r_c:.9f},  α = {a_c:.9f} rad "
           f"({math.degrees(a_c):.3f}°)")
    report(f"      B/A   = {B_c / A_c:.12f}   (√2, err "
           f"{B_c / A_c - B_CANON:+.1e})")
    report(f"      θ     = {th_c:.12f}   (2/9, err "
           f"{th_c - THETA_CANON:+.1e})")
    report(f"      gaps  = {np.round(d_c, 7)}   Q = {Q_c:.12f}")
    report(f"      m₂/m₁ = {ms[1] / ms[0]:.4f}   "
           f"m₃/m₁ = {ms[2] / ms[0]:.4f}")
    assert abs(B_c / A_c - B_CANON) < 1e-9
    assert abs(th_c - THETA_CANON) < 1e-9
    assert Q_c is not None and abs(Q_c - 2.0 / 3.0) < 1e-9

    report("    WHAT THIS IS: the canonical point is REACHABLE by the "
           "substrate.")
    report("      The anomalous channel is what reaches it: at α = 0, "
           "θ = 0 forces")
    report("      degeneracy and B/A ≥ 2, and no coupling reaches it "
           "at all.")
    report("    WHAT THIS IS NOT: a prediction.  Two knobs (g₁/g₀, α) "
           "are solved")
    report("      against two targets (B/A, θ), so the system is exactly "
           "determined")
    report("      and has ZERO predictive content.  m₂/m₁ and m₃/m₁ "
           "follow from")
    report("      B/A and θ; they are not independent confirmations.  "
           "Nothing here")
    report("      derives α, and Q is θ-blind (PART A), so a whole "
           "one-parameter")
    report("      family reproduces Q = 2/3 with θ free.")
    report(f"    AND THE BRANCH EXCLUSIVITY SURVIVES: max growth = "
           f"{g_c:.6f} here and")
    report(f"      {g_c2:.6f} at the coordinate copy (r, α) = "
           f"({-r_c:.6f}, {a_c - math.pi:.6f}).")
    report("    REDUNDANCY (correction): (g₁, α) → (−g₁, α−π) maps the "
           "interaction")
    report("      functional to itself identically (max deviation "
           f"{redundancy_dev:.1e}), so the two")
    report("      'inverse images' are ONE theory in two coordinate "
           "patches: 'both")
    report("      stable' is one fact stated twice, and NO dynamical "
           "selection between")
    report("      images exists.  The canonical point is modulationally "
           "STABLE under")
    report("      the full 6x6 BdG, and direct dynamics agrees "
           "(PART G).")
    report("      A per-channel-pairing reading reports instability "
           "here; that was an")
    report("      artifact of per-channel BdG pairing, valid only at "
           "α = 0.  The full")
    report("      6x6 statement: the Koide-"
           "readable branch and")
    report("      the knot-forming branch are DISJOINT for EVERY "
           "admissible (r, α).")
    report("      The canonical point is a STABLE uniform phase whose "
           "linearised")
    report("      spectrum is the Koide circulant; knots cannot form "
           "there from noise.")
    report("    STATUS: the canonical uniform state is stable.")
    report("      The effective-coupling flow and confined-phase route")
    report("      have no implementation.")
    out["alpha_canon"] = a_c
    out["r_alpha_canon"] = r_c
    out["BA_alpha"] = B_c / A_c
    out["theta_alpha"] = th_c
    out["alpha_max_growth"] = g_c
    out["alpha_point_unstable"] = bool(g_c > 1e-8)
    out["alpha_images_all_stable"] = bool(g_c < 1e-8 and g_c2 < 1e-8)
    out["alpha_point_positive_gaps"] = bool((d_c > 0).all())
    out["inverse_image_redundancy_dev"] = redundancy_dev
    # the corrected headline: the canonical point is STABLE with gaps
    # positive, and its 'two inverse images' are one theory under the
    # exact redundancy (g1, alpha) -> (-g1, alpha - pi).  The
    # exclusivity survives the complex phase.
    assert g_c < 1e-8 and g_c2 < 1e-8 and (d_c > 0).all()
    assert redundancy_dev < 1e-9

    ## the closed form of α ##
    # Eliminating (r, sin α) from the hermitian circulant equations at
    # the canonical targets leaves a BIQUADRATIC for c = cos α with
    # coefficients rational in t = tan θ alone:
    #   (729t⁴+154t²+1)c⁴ − (108t⁴+158t²+2)c² + (4t⁴+4t²+1) = 0.
    # α is an algebraic function of the Sugawara angle from two targets
    # and two variables.  Its nearest committed angle-menu entry is 0.20
    # away.  The transformation (-r, alpha - pi) is an exact coordinate
    # redundancy.
    t_ = math.tan(THETA_CANON)
    _a = 729 * t_ ** 4 + 154 * t_ ** 2 + 1
    _b = 108 * t_ ** 4 + 158 * t_ ** 2 + 2
    _d = _b * _b - 4 * _a * (4 * t_ ** 4 + 4 * t_ ** 2 + 1)
    assert _d > 0
    root_lo = (_b - math.sqrt(_d)) / (2 * _a)
    assert abs(math.cos(a_c) ** 2 - root_lo) < 1e-9
    out["alpha_algebraic_in_theta"] = True
    report("    closed form: cos²α solves the biquadratic in tan θ "
           "(no new constant);")
    report("    α is off the committed menu: the verification target is "
           "the dynamical route")
    report("    to the canonical point (not a selection: the images "
           "are one theory).")

    ## PART G: direct dynamics as arbiter ##
    # The 6x6 BdG and the per-channel form disagree at the
    # canonical point (stable vs sigma = 0.134).  The integrator
    # decides: seed the mode the per-channel form calls unstable
    # (q = 1 pattern, p inside its claimed band) at the canonical
    # point and measure its growth; positive control at the alpha = 0
    # knot-forming coupling, where both forms agree (rate = g1 rho0).
    report("\n  PART G: direct dynamics arbitration")
    Lbox, Nx = 40.0, 256
    xg = np.linspace(0.0, Lbox, Nx, endpoint=False)
    dxg = Lbox / Nx
    wph = np.exp(2j * np.pi / 3)

    def seeded_rate(r_run, alpha_run, mode, t_chunk=5.0, n_chunks=5):
        p = 2 * np.pi * mode / Lbox
        psi = np.ones((3, Nx), dtype=complex)
        for k in range(3):
            psi[k] += 1e-6 * (wph ** k) * np.cos(p * xg)
        amps = []
        for _ in range(n_chunks):
            psi = _nls._split_step(psi, dxg, 0.005, 1.0, r_run,
                                   n_steps=int(t_chunk / 0.005),
                                   alpha=alpha_run)
            q1 = (psi[0] + np.conj(wph) * psi[1]
                  + np.conj(wph ** 2) * psi[2]) / 3.0
            amps.append(abs(np.fft.fft(q1)[mode]) / Nx)
        rates = [math.log(amps[i + 1] / amps[i]) / t_chunk
                 for i in range(len(amps) - 1)]
        return float(np.mean(rates))

    rate_canonical = seeded_rate(r_c, a_c, mode=2)      # p = 0.314
    rate_control = seeded_rate(1.0 / math.sqrt(2.0), 0.0,
                               mode=12, t_chunk=2.0)    # p = 1.885,
    #                            the band's peak (band p in 1.46-2.22)
    report(f"    canonical point, seeded mode p = 0.314: measured "
           f"rate {rate_canonical:+.4f} (per-channel form claims "
           f"+0.13)")
    report(f"    control (r = 1/sqrt2, alpha = 0), p = 1.885: measured "
           f"rate {rate_control:+.4f} (theory g1 rho0 = 0.707)")
    assert abs(rate_canonical) < 0.02      # stable, as the 6x6 says
    assert 0.5 < rate_control < 0.85       # the instability is real
    out["dynamics_rate_canonical"] = rate_canonical
    out["dynamics_rate_control"] = rate_control

    ####### verdict #######
    report("\n  VERDICT")
    report("    CANON       : the only claimant on the positive-gap "
           "branch with Q = 2/3")
    report("    HARMONIC    : (B, θ) = (1, 0) → Q = ½, two masses "
           "degenerate. RETIRED")
    report("    BdG UNIFORM : at α = 0, θ = 0 forces degeneracy and "
           "B/A = √2 needs")
    report("                  g₁ < 0, which kills the instability. "
           "Disjoint branches.")
    report("    the three are NOT reconcilable at α = 0; the amplitude "
           "is pinned by Q")
    report("    alone, and the phase must come from elsewhere "
           "(Sugawara).")
    report("  AND WITH THE COMPLEX PHASE: a complex Z₃ coupling IS "
           "hermitian-admissible")
    report("    and reaches the canonical point, but BOTH inverse "
           "images are STABLE")
    report("    (PART F, 6x6 BdG; PART G, direct dynamics), so the "
           "disjointness holds")
    report("    for every admissible (r, α): the canonical point is a "
           "stable uniform")
    report("    phase carrying the Koide circulant, and the "
           "knot-forming phase cannot")
    report("    carry it at uniform background.")
    report("    The effective-coupling route from the α = 0 phase to the")
    report("    canonical circulant has no implementation.")

    # reconciliation would require instability at a coupling that
    # reaches B/A = √2; none exists at α = 0 (PART D, growth g_n) or
    # at either α ≠ 0 inverse image (PART F); computed
    out["reconcilable"] = bool(g_n > 1e-8 or g_c > 1e-8 or g_c2 > 1e-8)
    return out


def test_matches_live_substrate():
    """The closed-form circulant must track nls_soliton's own
    interaction matrix, not a transcription of it."""
    for alpha, r in ((0.0, 1 / math.sqrt(2)), (2 / 9, 0.5),
                     (4.333441893, 0.383029723), (1.0, 0.9)):
        c0, B, th, _ = bdg_circulant(1.0, r, alpha=alpha)
        p0 = np.ones(3, dtype=complex)

        def F(p):
            M = _nls.interaction_matrix(p.reshape(3, 1), 1.0, r, alpha)[0]
            return M @ p

        mu = (F(p0)[0] / p0[0]).real
        h = 1e-7
        A = np.zeros((3, 3), dtype=complex)
        for j in range(3):
            e = np.zeros(3, dtype=complex)
            e[j] = 1
            dr = (F(p0 + h * e) - F(p0 - h * e)) / (2 * h)
            di = (F(p0 + 1j * h * e) - F(p0 - 1j * h * e)) / (2 * h)
            A[:, j] = 0.5 * (dr - 1j * di)
        A = A - mu * np.eye(3)
        assert abs(A[0, 0].real - c0) < 1e-6, (alpha, r)
        assert abs(2 * abs(A[0, 1]) - B) < 1e-6, (alpha, r)
    return True


if __name__ == "__main__":
    run()
