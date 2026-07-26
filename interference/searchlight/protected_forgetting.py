"""
Protected forgetting: the Z₃ harmonic, its two moments, and what
actually fixes it.

THE HARMONIC.  The Z₃-symmetric substrate has three components ψ_k
(k=0,1,2).  Any Z₃-equivariant insertion is a real circulant, so its
eigenvalues take the one-parameter-pair form

    v_k = B cos(θ + 2πk/3),      k = 0, 1, 2

with amplitude B and phase θ.  Every such v is traceless automatically
(Σ cos(θ + 2πk/3) = 0), so PvP = 0 holds for the whole family.  The
physics is in (B, θ), and the canon fixes them SEPARATELY:

    B = √d₁₀ = √2      octonionic CG ratio |C₃̄|/|C₁| (octonions.py),
                       equivalently the Koide value Q₀ = 2/3
    θ = h₁₀ = 2/9      Sugawara conformal weight of the SU(3)₃
                       fundamental (root.py)

so v = (+1.379438, −0.959650, −0.419788),  ‖v‖² = 3.

═══════════════════════════════════════════════════════════════════════
THE TWO MOMENTS, AND WHY THEY DIFFER BY A FACTOR OF TWO
═══════════════════════════════════════════════════════════════════════

Let P = |u⟩⟨u| project onto the uniform vacuum u = (1,1,1)/√3, and let
v act as a diagonal operator (elementwise multiplication).

  PvP = 0        FORGETTING.  u·(v∘u) = Σv_k/3 = 0, for every (B, θ).
                 Forces the minimal visible closure to be QUADRATIC,
                 which is the |w₁| = 2 length-rule seed.

  Pv²P = (B²/2)P TOTAL second moment.  u·(v²∘u) = Σv_k²/3 = B²/2,
                 because Σcos²(θ + 2πk/3) = 3/2 for every θ.
                 At the canonical B = √2 this equals 1, NOT ½.

  E[v²] ≡ |v̂₁|² = B²/4    PER-MODE second moment.  The Z₃ Fourier
                 component v̂₁ = (1/3)Σ_k v_k ω^(−k) = (B/2)e^(iθ), so
                 the q=1 channel carries |v̂₁|² = B²/4 = ½ at B = √2.
                 The q=2 channel carries the same, and
                     Pv²P = |v̂₁|² + |v̂₂|² = 2·E[v²].

  The quantity consumed by the gravity/Higgs ledger (root.py `E_v2`,
  gravity.py ξ_bridge = α_G₂·E[v²]·h_bridge = 1/(48π)) is the PER-MODE
  power, because the bridge couples one Z₃ channel at a time.  It is ½
  at the canonical amplitude.  The TOTAL is 1.  Both are correct; they
  are different contractions; the notation "Pv²P = ½P" conflates
  them.  No numerical value in the web changes: E[v²] = ½ throughout.

═══════════════════════════════════════════════════════════════════════
THE KOIDE IDENTITY, AND WHAT IT PINS
═══════════════════════════════════════════════════════════════════════

On the positive-gap branch (√m_k = 1 + v_k > 0 for all k):

    Q ≡ Σm_k / (Σ√m_k)²  =  1/3 + B²/6  =  1/3 + (Pv²P)/3
                          =  1/3 + 2·E[v²]/3

Q is INDEPENDENT of θ.  Hence:

    Q = 2/3  ⟺  B = √2  ⟺  Pv²P = 1  ⟺  E[v²] = ½  ⟺  ‖v‖² = 3

A normalisation therefore fixes the AMPLITUDE and can never fix the
phase.  θ is fixed by Sugawara, not by any moment condition.

═══════════════════════════════════════════════════════════════════════
RETIRED: the long-root normalisation argument
═══════════════════════════════════════════════════════════════════════

The long-root normalisation argument selects v = (1, −½, −½), the
G₂ long root at ‖v‖² = 3/2, by requiring the TOTAL moment to equal
½.  That argument is retired, and the scan below shows why, rather
than asserting it:

  • v = (1, −½, −½) is (B, θ) = (1, 0), giving √m = (2, ½, ½):
    two masses EXACTLY degenerate, ratios 1 : 1 : 16, and Q = ½.
  • θ = 0 forces √m = (1+B, 1−B/2, 1−B/2) for EVERY B, so the
    θ = 0 family can never carry three generations.
  • the short-root direction, scaled to the same norm, gives three
    DISTINCT positive masses (1 : 55.7 : 194); a discriminant
    printed without consulting the computed ratios reads the
    opposite.

What survives from that argument is PvP = 0, which holds for the whole
family and needs no normalisation.

═══════════════════════════════════════════════════════════════════════
WHAT THIS FILE PROVES / WHAT REMAINS OPEN
═══════════════════════════════════════════════════════════════════════

PROVES (computed here):
  • PvP = 0 for every (B, θ): forgetting is generic, not selected.
  • Pv²P = B²/2 and E[v²] = B²/4 exactly; the factor 2 is the
    two-channel Z₃ split, verified by explicit Fourier transform.
  • Q = 1/3 + B²/6, θ-independent, to machine precision.
  • θ = 0 ⟹ exact two-fold degeneracy, for every B.
  • at the Koide-forced B = √2, positivity of all three masses
    restricts θ to a computed window; θ = 2/9 lies inside it, and the
    margin is thin: the electron sits near the positivity boundary.

OPEN (labelled, not claimed):
  • B = √2 comes from octonions.py (Fano CG); θ = 2/9 comes from
    Sugawara.  Nothing in this file derives either; it verifies their
    consequences and the moment bookkeeping around them.
  • no dynamical state of the Z₃-NLS is known to realise (B, θ) =
    (√2, 2/9); the substrate's contact interaction is real-projected,
    hence θ = 0, hence degenerate.  See tests/probes/
    circulant_normalization.py for the three-circulant comparison.

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import math

import numpy as np

# The canonical harmonic parameters (canonical values live in
# root.py; duplicated here as literals so this file stays runnable
# standalone).  DRIFT GUARD: when the engine is importable (every
# test and run context) the literals are asserted equal to root's
# values, so the two copies cannot diverge silently
# (tests/test_probes.py::test_harmonic_single_source freezes the
# guard's presence).
B_CANON = math.sqrt(2.0)        # = √d₁₀, octonionic CG ratio
THETA_CANON = 2.0 / 9.0         # = h₁₀, Sugawara weight of the fundamental

try:                            # standalone runs lack the engine path
    import root as _root
except ImportError:
    _root = None
if _root is not None:
    assert abs(B_CANON - _root.BA_ratio) < 1e-15, \
        "harmonic drift: B_CANON != root.BA_ratio"
    assert abs(THETA_CANON - float(_root.h10)) < 1e-15, \
        "harmonic drift: THETA_CANON != root.h10"


def harmonic(B=B_CANON, theta=THETA_CANON):
    """The Z₃ harmonic v_k = B cos(theta + 2πk/3) as a 3-vector."""
    return np.array([B * math.cos(theta + 2 * math.pi * k / 3)
                     for k in range(3)])


def moments(v):
    """Return (PvP, Pv2P_total, E_v2_mode) for a 3-vector harmonic.

    P = |u><u| with u = (1,1,1)/√3; v acts by elementwise multiplication.

      PvP         = Σv_k / 3
      Pv2P_total  = Σv_k² / 3                    (both Z₃ channels)
      E_v2_mode   = |v̂₁|², v̂₁ = (1/3)Σ_k v_k ω^(−k)   (one channel)

    The identity Pv2P_total = 2·E_v2_mode is asserted.
    """
    u = np.ones(3) / math.sqrt(3.0)
    PvP = float(np.dot(u, v * u))
    Pv2P_total = float(np.dot(u, v ** 2 * u))

    omega = np.exp(2j * np.pi / 3)
    v_hat1 = sum(v[k] * omega ** (-k) for k in range(3)) / 3.0
    E_v2_mode = float(abs(v_hat1) ** 2)

    assert abs(Pv2P_total - 2.0 * E_v2_mode) < 1e-12, (
        f"channel split broken: {Pv2P_total} != 2*{E_v2_mode}")
    return PvP, Pv2P_total, E_v2_mode


def koide_spectrum(v):
    """Return (sqrt_m, masses, Q) for √m_k = 1 + v_k, or None if any
    gap is non-positive (off the positive-gap branch)."""
    sq_m = 1.0 + v
    if any(s <= 0 for s in sq_m):
        return None
    masses = sq_m ** 2
    Q = float(masses.sum() / sq_m.sum() ** 2)
    return sq_m, masses, Q


def mass_ratios(v):
    """Sorted mass ratios normalised to the lightest, or None."""
    spec = koide_spectrum(v)
    if spec is None:
        return None
    _, masses, _ = spec
    m_min = float(masses.min())
    return sorted(float(m) / m_min for m in masses)


def positivity_window(B, n_scan=36000):
    """Measure of θ ∈ [0, 2π) for which all three gaps are positive,
    at fixed amplitude B, plus the margin at θ = THETA_CANON.

    All gaps positive ⟺ min_k cos(θ + 2πk/3) > −1/B.
    (The default n_scan matches the derive() call site.)
    """
    n_ok = 0
    for j in range(n_scan):
        th = 2 * math.pi * j / n_scan
        if (1.0 + harmonic(B, th)).min() > 0:
            n_ok += 1
    frac = n_ok / n_scan
    v_can = harmonic(B, THETA_CANON)
    min_cos = float(min(math.cos(THETA_CANON + 2 * math.pi * k / 3)
                        for k in range(3)))
    return frac, min_cos, float((1.0 + v_can).min())


def derive(verbose=True):
    """Verify the harmonic's moments, the Koide identity, and the
    positivity window.  Returns a dict of computed (not literal) values."""
    print("\n" + "=" * 78)
    print("  PROTECTED FORGETTING: Z₃ HARMONIC, MOMENTS, AND KOIDE")
    print("=" * 78)

    # ═══════════════════════════════════════════════════════════════════
    #  1. The canonical harmonic and its two moments
    # ═══════════════════════════════════════════════════════════════════

    v = harmonic()
    PvP, Pv2P_total, E_v2_mode = moments(v)
    sq_m, masses, Q = koide_spectrum(v)

    print(f"\n  canonical harmonic  B = √d₁₀ = {B_CANON:.9f},  "
          f"θ = h₁₀ = {THETA_CANON:.9f}")
    print(f"    v      = ({v[0]:+.6f}, {v[1]:+.6f}, {v[2]:+.6f})")
    print(f"    ‖v‖²   = {float(np.dot(v, v)):.9f}   (= 3, not 3/2)")
    print(f"    PvP        = {PvP:+.15f}   (0, for every B and θ)")
    print(f"    Pv²P total = {Pv2P_total:.15f}   (= B²/2 = 1)")
    print(f"    E[v²] mode = {E_v2_mode:.15f}   (= B²/4 = ½, the "
          f"quantity the web consumes)")
    print(f"    Pv²P = 2·E[v²]  ✓  (two Z₃ channels, q = 1 and q = 2)")

    assert abs(PvP) < 1e-14
    assert abs(Pv2P_total - 1.0) < 1e-12
    assert abs(E_v2_mode - 0.5) < 1e-12

    print(f"\n    √m_k = 1 + v_k = ({sq_m[0]:.9f}, {sq_m[1]:.9f}, "
          f"{sq_m[2]:.9f})")
    print(f"    Koide Q = {Q:.15f}   (target 2/3, "
          f"|ΔQ| = {abs(Q - 2/3):.2e})")
    assert abs(Q - 2.0 / 3.0) < 1e-12

    # ═══════════════════════════════════════════════════════════════════
    #  2. Q = 1/3 + B²/6 = 1/3 + (Pv²P)/3, and it is θ-independent
    # ═══════════════════════════════════════════════════════════════════

    print(f"\n  identity  Q = 1/3 + B²/6 = 1/3 + Pv²P/3 = 1/3 + 2E[v²]/3")
    print(f"    {'B':>8} {'θ':>8} {'Pv²P':>10} {'E[v²]':>10} "
          f"{'Q(measured)':>14} {'Q(identity)':>14}")
    max_dev = 0.0
    for B_t in (0.3, 1.0, B_CANON, 1.7):
        for th_t in (0.0, 0.7, 2.2):
            v_t = harmonic(B_t, th_t)
            _, tot_t, mode_t = moments(v_t)
            spec_t = koide_spectrum(v_t)
            if spec_t is None:
                continue
            Q_t = spec_t[2]
            Q_id = 1.0 / 3.0 + B_t ** 2 / 6.0
            max_dev = max(max_dev, abs(Q_t - Q_id))
            print(f"    {B_t:8.4f} {th_t:8.4f} {tot_t:10.6f} "
                  f"{mode_t:10.6f} {Q_t:14.10f} {Q_id:14.10f}")
    print(f"    max |Q_measured − Q_identity| = {max_dev:.2e}  "
          f"(θ-independence verified)")
    assert max_dev < 1e-12

    # ═══════════════════════════════════════════════════════════════════
    #  3. θ = 0 is degenerate for every B  (kills the long-root harmonic)
    # ═══════════════════════════════════════════════════════════════════

    print(f"\n  θ = 0 degeneracy theorem: √m = (1+B, 1−B/2, 1−B/2)")
    theta0_degenerate = True
    for B_t in (0.2, 0.5, 1.0, B_CANON, 1.9):
        v0 = harmonic(B_t, 0.0)
        split = abs(v0[1] - v0[2])
        theta0_degenerate &= (split < 1e-12)
        print(f"    B = {B_t:6.4f}:  v = ({v0[0]:+.6f}, {v0[1]:+.6f}, "
              f"{v0[2]:+.6f})   |v₁−v₂| = {split:.2e}")
    assert theta0_degenerate
    print(f"    → the θ = 0 family cannot carry three generations, "
          f"for ANY amplitude")

    # The retired harmonic, computed rather than asserted.
    v_old = np.array([1.0, -0.5, -0.5])
    PvP_old, tot_old, mode_old = moments(v_old)
    r_old = mass_ratios(v_old)
    Q_old = koide_spectrum(v_old)[2]
    print(f"\n  retired harmonic v = (1, −½, −½) = (B, θ) = (1, 0):")
    print(f"    Pv²P total = {tot_old:.6f}   E[v²] = {mode_old:.6f}")
    print(f"    mass ratios 1 : {r_old[1]:.3f} : {r_old[2]:.3f}"
          f"   → two masses degenerate")
    print(f"    Q = {Q_old:.6f}  (= ½, not 2/3)")
    assert abs(Q_old - 0.5) < 1e-12
    assert abs(r_old[1] - 1.0) < 1e-9, "retired harmonic must be degenerate"

    # ═══════════════════════════════════════════════════════════════════
    #  4. The short-root claim, computed
    # ═══════════════════════════════════════════════════════════════════

    print(f"\n  G₂ short-root direction (1, −1, 0) scaled to ‖v‖² = 3/2:")
    v_short_raw = np.array([1.0, -1.0, 0.0])
    v_short = v_short_raw * math.sqrt(1.5 / float(np.dot(v_short_raw,
                                                         v_short_raw)))
    PvP_s, tot_s, mode_s = moments(v_short)
    r_s = mass_ratios(v_short)
    print(f"    v = ({v_short[0]:+.6f}, {v_short[1]:+.6f}, "
          f"{v_short[2]:+.6f})")
    print(f"    PvP = {PvP_s:+.2e}   Pv²P total = {tot_s:.6f}")
    if r_s is None:
        short_three_distinct = False
        print(f"    at least one mass is non-positive")
    else:
        short_three_distinct = (r_s[1] / r_s[0] > 1.01
                                and r_s[2] / r_s[1] > 1.01)
        print(f"    mass ratios 1 : {r_s[1]:.3f} : {r_s[2]:.3f}")
        print(f"    three distinct positive masses: "
              f"{'YES' if short_three_distinct else 'NO'}")
    print(f"    NOTE: the short root gives three DISTINCT masses; a "
          f"degeneracy")
    print(f"    discriminant would select against the working "
          f"direction.  Both")
    print(f"    directions are ruled out for the real reason instead: "
          f"neither has")
    print(f"    B = √2, so neither gives Q = 2/3.")

    # ═══════════════════════════════════════════════════════════════════
    #  5. Positivity window at the Koide-forced amplitude
    # ═══════════════════════════════════════════════════════════════════

    print(f"\n  positivity window at B = √2 (all three gaps > 0):")
    frac, min_cos, gap_min = positivity_window(B_CANON, n_scan=36000)
    print(f"    admissible fraction of θ ∈ [0, 2π):  {frac:.4f}")
    print(f"    boundary condition: min_k cos(θ + 2πk/3) > −1/√2 = "
          f"{-1/math.sqrt(2):.6f}")
    print(f"    at θ = 2/9:  min_k cos = {min_cos:.6f}  "
          f"(margin {min_cos + 1/math.sqrt(2):+.6f})")
    print(f"    smallest gap √m_min = {gap_min:.6f}  → the electron sits "
          f"near the")
    print(f"    positivity boundary; the lepton hierarchy is that "
          f"thinness.")
    assert gap_min > 0, "canonical harmonic must be on the positive-gap branch"

    # ═══════════════════════════════════════════════════════════════════
    #  6. Weyl orbit (Z₃ relabelling leaves the spectrum alone)
    # ═══════════════════════════════════════════════════════════════════

    print(f"\n  Z₃ relabelling orbit of the canonical harmonic:")
    orbit_Q = []
    for shift in range(3):
        v_w = np.roll(v, shift)
        Q_w = koide_spectrum(v_w)[2]
        orbit_Q.append(Q_w)
        print(f"    roll {shift}: ({v_w[0]:+.4f}, {v_w[1]:+.4f}, "
              f"{v_w[2]:+.4f})   Q = {Q_w:.12f}")
    assert max(orbit_Q) - min(orbit_Q) < 1e-12
    print(f"    all rolls give identical Q (Z₃ symmetry)")

    # ═══════════════════════════════════════════════════════════════════
    #  7. Summary
    # ═══════════════════════════════════════════════════════════════════

    unique_amplitude = abs(Q - 2.0 / 3.0) < 1e-12 and abs(
        B_CANON ** 2 - 2.0) < 1e-12

    print(f"\n  VERIFIED:")
    print(f"    ✓ PvP = 0 for the whole (B, θ) family: forgetting is "
          f"generic")
    print(f"    ✓ Pv²P = B²/2 = 1 (total), E[v²] = B²/4 = ½ (per mode)")
    print(f"    ✓ Q = 1/3 + B²/6, θ-independent → Q = 2/3 pins B = √2 "
          f"ALONE")
    print(f"    ✓ θ = 0 degenerate for every B → long-root harmonic excluded")
    print(f"    ✓ short-root claim computed")
    print(f"    ✓ θ = 2/9 lies inside the positivity window, thin margin")
    print(f"  OPEN:")
    print(f"    · B = √2 from octonions.py, θ = 2/9 from Sugawara; "
          f"neither derived here")
    print(f"    · no dynamical Z₃-NLS state realises (√2, 2/9); the "
          f"contact")
    print(f"      interaction is real-projected, hence θ = 0, hence "
          f"degenerate")

    return {
        'B': B_CANON,
        'theta': THETA_CANON,
        'PvP': PvP,
        'Pv2P_total': Pv2P_total,
        'E_v2_mode': E_v2_mode,
        'Q': Q,
        'sqrt_masses': [float(s) for s in sq_m],
        'mass_ratios': mass_ratios(v),
        'theta0_degenerate': bool(theta0_degenerate),
        'retired_Q': Q_old,
        'short_root_three_distinct': bool(short_three_distinct),
        'positivity_fraction': frac,
        'gap_min': gap_min,
        'amplitude_pinned_by_koide': bool(unique_amplitude),
    }


if __name__ == "__main__":
    derive()
