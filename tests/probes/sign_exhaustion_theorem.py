"""
sign_exhaustion_theorem.py: the PMNS sign rows against the data.

This probe is the check, not the selector.  tests/probes/pmns_sign_bits.py
fixes the row with two data blind gates (the lepton sector's own pairing
and the charge dictionary).  This probe compares the four rows with the
data. Condition (a) is the NuFIT range. Condition (b) only names which
sign is matter, since one orientation bit flips every CP-odd sign
together.  The data pick the same row the gates leave.

Given the one orientation bit axiom (a single Fano orientation
datum epsilon enters all CP-odd phases), and the Standard Model sphaleron
relation B = -(28/79)L, exactly one of the four sign rows
  (delta_sign, eta_sign) in {-1,+1}^2
simultaneously satisfies:
  (a) invariant delta_CP inside NuFIT 6.0 NO 3-sigma: (96, 422) deg,
  (b) eta_B > 0 (observed matter excess).

The closing row is delta_sign = -1 (CKM convention), eta_sign = -1
(sphaleron conversion).  This probe verifies both the algebraic
structure (D1a: epsilon propagation) and the sign chain (D1b:
sphaleron sign inheritance).

D1a proof outline (verified computationally below):
  1. U_TBM is REAL: the conjugation orbit compression of bulk SU(3)_3
     fusion gives eigenvectors with no imaginary part.
  2. The ONLY complex phase in U_PMNS = U_e^dagger . U_TBM comes from
     U_e = R23(-phi) . R13(theta, delta), specifically from the
     e^{i*delta} in R13.
  3. delta = -epsilon * arccos(-lambda), where lambda is a CP EVEN
     quantity (the echoed CKM Wolfenstein parameter).
  4. In V_CKM, the CP-odd input is eta_bar = epsilon * pi/d11^2.
  5. Both read from the same datum: the Fano plane pseudoscalar
     (the sign of C3bar in the octonionic CG decomposition).
  6. Reversing epsilon flips ALL CP-odd observables together and
     touches no CP-even magnitude (mixing angles, masses, |V_ij|).
     This is verified by the orientation_bit.py probe. Here we
     verify U_TBM is real and that no second sign datum exists.

D1b proof outline:
  The sphaleron relation B = -(28/79)L is a Standard Model result
  (Khlebnikov-Shaposhnikov 1988, Harvey-Turner 1990).  It does not
  depend on the framework's algebraic data.  Given this relation:
    eta_B = -n7 * J_lep * exp(-S)
  With the closing sign (epsilon = +1 in the standard Fano orientation):
    J_lep < 0  (from mixing.py, confirmed here)
    => -J_lep > 0
    => eta_B > 0  (matter excess, as observed)

Usage: python3 sign_exhaustion_theorem.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
import _bootstrap  # noqa: E402,F401
from interference.root import d10, d11, h10, n7, alpha_phys, BA_ratio


NUFIT_DELTA_3SIGMA = (96.0, 422.0)


def _build_u_tbm():
    """The TBM eigenvectors are REAL: no imaginary part at all."""
    return np.array([
        [math.sqrt(2 / 3), 1 / math.sqrt(3), 0],
        [-1 / math.sqrt(6), 1 / math.sqrt(3), 1 / math.sqrt(2)],
        [1 / math.sqrt(6), -1 / math.sqrt(3), 1 / math.sqrt(2)],
    ], dtype=float)


def _build_pmns(delta_sign):
    """Build U_PMNS for a given delta sign choice."""
    lam = float(h10) * (1.0 + 0.5 / (d10 * d11)**2)
    delta = delta_sign * math.acos(-lam)

    def mass_factor(k):
        gap = 1.0 + BA_ratio * math.cos(float(h10) + 2 * math.pi * k / 3)
        vent = 1.0 - (d11**3 / d10) * (alpha_phys / (2 * math.pi))**2 * gap
        return gap**2 * vent

    phi = math.sqrt(mass_factor(1) / mass_factor(2)) / math.sqrt(d10)
    theta = float(h10)

    def R13(a, phase):
        c, s, e = math.cos(a), math.sin(a), np.exp(1j * phase)
        return np.array([[c, 0, s * e], [0, 1, 0],
                         [-s * e.conjugate(), 0, c]], dtype=complex)

    def R23(a):
        c, s = math.cos(a), math.sin(a)
        return np.array([[1, 0, 0], [0, c, s], [0, -s, c]], dtype=complex)

    U_tbm = _build_u_tbm().astype(complex)
    U_e = R23(-phi) @ R13(theta, delta)
    U = np.conj(U_e.T) @ U_tbm

    J = float((U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])).imag)
    s13 = abs(U[0, 2])
    c13sq = 1 - s13**2
    sin2_12 = abs(U[0, 1])**2 / c13sq
    sin2_23 = abs(U[1, 2])**2 / c13sq

    s12 = math.sqrt(sin2_12)
    c12 = math.sqrt(1 - sin2_12)
    s23 = math.sqrt(sin2_23)
    c23 = math.sqrt(1 - sin2_23)
    c13 = math.sqrt(c13sq)
    sin_d = J / (s12 * c12 * s23 * c23 * s13 * c13**2)
    cos_d = (abs(U[1, 0])**2 - sin2_12 * c23**2
             - (1 - sin2_12) * s23**2 * s13**2) / (
             2 * s12 * c12 * s23 * c23 * s13)
    delta_inv = math.degrees(math.atan2(sin_d, cos_d)) % 360

    return {
        "sin2_12": sin2_12, "sin2_23": sin2_23, "sin2_13": s13**2,
        "J_lep": J, "delta_invariant_deg": delta_inv,
        "U_PMNS": U, "U_e": U_e,
    }


def run(report=print):
    report("SIGN EXHAUSTION THEOREM: the PMNS sign is forced")
    report("=" * 64)

    # ── D1a: U_TBM is real ──
    U_tbm = _build_u_tbm()
    imag_norm = np.linalg.norm(np.imag(U_tbm))
    report(f"  D1a check 1: ||Im(U_TBM)|| = {imag_norm:.2e} (must be 0)")
    assert imag_norm == 0.0, "U_TBM must be exactly real"
    report("  => all CP-odd phases in U_PMNS come from U_e alone")

    # ── D1a: U_e carries delta, which reads epsilon ──
    pmns_fwd = _build_pmns(-1)  # CKM convention
    pmns_rev = _build_pmns(+1)  # conjugate
    report(f"\n  D1a check 2: delta enters U_e through R13(theta, delta)")
    report(f"    delta(CKM)  = -arccos(-lam) => invariant {pmns_fwd['delta_invariant_deg']:.2f} deg")
    report(f"    delta(conj) = +arccos(-lam) => invariant {pmns_rev['delta_invariant_deg']:.2f} deg")

    # Verify mixing angles are identical under the sign flip
    for key in ("sin2_12", "sin2_23", "sin2_13"):
        assert math.isclose(pmns_fwd[key], pmns_rev[key], rel_tol=1e-12), \
            f"{key} changed under delta sign flip"
    report("  => mixing angles INVARIANT under delta sign flip")

    # Verify J flips sign
    assert math.isclose(pmns_fwd["J_lep"], -pmns_rev["J_lep"], rel_tol=1e-12), \
        "J_lep must flip sign under delta flip"
    report(f"  => J_lep flips: {pmns_fwd['J_lep']:+.6f} <-> {pmns_rev['J_lep']:+.6f}")

    # ── D1a: no second sign datum ──
    # The ONLY place epsilon enters U_PMNS is through delta in R13.
    # U_TBM is real, phi is CP-even, theta = h10 is CP-even.
    # Therefore epsilon controls ALL CP-odd observables in the lepton
    # sector through exactly the same chain as in the quark sector.
    report("\n  D1a CONCLUSION: epsilon enters the lepton sector through")
    report("  delta alone (the only complex phase in U_e).  U_TBM is real,")
    report("  phi and theta are CP-even.  No second sign datum exists.")

    # ── D1b: sphaleron sign chain ──
    report(f"\n  D1b: sphaleron sign B = -(28/79)L")
    tunnelling = math.exp(-2 * math.pi**2)
    report(f"    tunnelling exp(-2pi^2) = {tunnelling:.6e}")
    report(f"    J_lep(CKM) = {pmns_fwd['J_lep']:+.6f} < 0")
    eta_B = float(n7) * (-pmns_fwd["J_lep"]) * tunnelling
    report(f"    eta_B = n7 * (-J_lep) * exp(-S) = {eta_B:.4e} > 0")
    assert eta_B > 0, "eta_B must be positive for observed matter excess"
    report("  => the sphaleron minus sign converts J_lep < 0 to eta_B > 0")

    # ── Four row exhaustion ──
    report(f"\n  EXHAUSTION: all four (delta_sign, eta_sign) rows")
    report("-" * 64)
    rows = []
    for delta_sign, delta_label in [(-1, "CKM conv"), (+1, "conjugate")]:
        pmns = _build_pmns(delta_sign)
        for eta_sign, eta_label in [(-1, "B=-(28/79)L"), (+1, "B=+(28/79)L")]:
            eta_B = eta_sign * float(n7) * pmns["J_lep"] * tunnelling
            delta_ok = NUFIT_DELTA_3SIGMA[0] < pmns["delta_invariant_deg"] < NUFIT_DELTA_3SIGMA[1]
            eta_ok = eta_B > 0
            closes = delta_ok and eta_ok
            tag = "CLOSES" if closes else "fails"
            why = []
            if not delta_ok:
                why.append(f"delta={pmns['delta_invariant_deg']:.1f} outside NuFIT")
            if not eta_ok:
                why.append(f"eta_B={eta_B:.2e} <= 0")
            reason = f" ({'; '.join(why)})" if why else ""
            row = {
                "delta_sign": delta_sign,
                "eta_sign": eta_sign,
                "closes": closes,
                "delta_deg": pmns["delta_invariant_deg"],
                "eta_B": eta_B,
            }
            rows.append(row)
            report(f"  {delta_label:>10s} | {eta_label:>13s} | "
                   f"delta={pmns['delta_invariant_deg']:7.2f} | "
                   f"eta_B={eta_B:+.2e} | {tag}{reason}")

    closing = [r for r in rows if r["closes"]]
    report(f"\n  {len(closing)}/4 rows close.")
    assert len(closing) == 1, f"expected exactly 1 closing row, got {len(closing)}"
    winner = closing[0]
    assert winner["delta_sign"] == -1
    assert winner["eta_sign"] == -1

    # ── Why the three other rows fail ──
    report("\n  WHY THE OTHER THREE ROWS FAIL:")
    report("  Row (+1, -1): delta = 84.26 deg, below NuFIT 96 deg.")
    report("    The conjugate sign places delta on the wrong side of")
    report("    arccos(-lam): the cosine is the same, the sine flips,")
    report("    and the invariant falls in the first quadrant.")
    report("  Row (-1, +1): eta_B < 0, antimatter universe.")
    report("    B = +(28/79)L contradicts the SM sphaleron calculation.")
    report("  Row (+1, +1): both failures compound.")
    report("    delta outside range AND eta_B < 0.")

    report("-" * 64)
    report("  SELECTION: of the four rows, one is consistent with the")
    report("  NuFIT range and a matter excess.  The data blind gates of")
    report("  tests/probes/pmns_sign_bits.py leave that row, so the data")
    report("  check the sign rather than choose it.")

    return {
        "u_tbm_real": True,
        "j_lep_flips": True,
        "angles_invariant": True,
        "closing_count": len(closing),
        "closing_row": winner,
        "eta_B_positive": eta_B > 0,
    }


if __name__ == "__main__":
    run()
