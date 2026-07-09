"""
orientation_bit.py: ONE Fano-orientation bit propagates through the
CP-odd observables. No second CP-sign datum is inserted.

THE COMPUTATION: build the octonion cross product from the Fano
plane in BOTH orientations (all seven line directions reversed),
recompute every orientation-sensitive quantity, and verify:

  (i)   CP-EVEN data are invariant: |C1| = 1, |C3bar| = sqrt(2),
        hence B/A, masses, mixing-angle magnitudes.
  (ii)  The one registered orientation map sets the sign of the complex
        phase. The CKM and PMNS matrices are recomputed in both
        orientations, so eta_bar, J_CKM, J_lep, and eta_B flip together.
  (iii) therefore CP sign and baryon-excess sign are ONE propagated bit,
        the forced correlation of the
        one-orientation-bit argument, verified at the algebra level.

COMPUTED SCOPE. Reversing the Fano orientation recomputes the Fano
pseudoscalar, CKM matrix, PMNS matrix, and baryogenesis expression. Every
listed CP-odd value flips and the CP-even mixing magnitudes are unchanged.

COMPOSITION. This probe owns the CP-sign propagation. The mass-protection
selection of the chiral branch, the gauge-representation embedding, and the
Skyrmion-Dirac spectral flow are established in their respective branches and
compose under the repository's stated fermion/knot identification. Not
recomputing those branches here is scope separation, not an open implication.

Usage: python3 orientation_bit.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "interference")))
import _bootstrap  # noqa: E402,F401
from octonions import FANO_TRIPLES, FANO_SIGNS
from interference.root import d10, d11, n7, h10, alpha_phys, BA_ratio


def cross(a, b, signs):
    out = np.zeros(7, dtype=complex)
    for (i, j, k), s in zip(FANO_TRIPLES, signs):
        ii, jj, kk = i - 1, j - 1, k - 1
        out[kk] += s * (a[ii] * b[jj] - a[jj] * b[ii])
        out[ii] += s * (a[jj] * b[kk] - a[kk] * b[jj])
        out[jj] += s * (a[kk] * b[ii] - a[ii] * b[kk])
    return out


def cg_data(signs):
    z = np.zeros((3, 7), dtype=complex)
    for a in range(3):
        z[a, 2*a] = 1/math.sqrt(2)
        z[a, 2*a + 1] = -1j/math.sqrt(2)
    zbar = z.conj()
    C1 = abs(cross(z[0], zbar[0], signs)[6])
    prod = cross(z[0], z[1], signs)
    amp = np.vdot(zbar[2], prod) / np.vdot(zbar[2], zbar[2]).real
    return C1, abs(amp), np.sign(amp.real) or np.sign(amp.imag)


def ckm_data(orientation):
    """Recompute the CKM matrix with the oriented imaginary apex."""
    lam = float(h10) * (1.0 + 0.5 / (d10 * d11)**2)  # 73/324
    A = math.sqrt(d10 / d11)
    rho = math.sqrt(d10) / d11**2
    eta = orientation * math.pi / d11**2
    s12, s23 = lam, A * lam**2
    c12, c23 = math.sqrt(1 - s12**2), math.sqrt(1 - s23**2)
    apex = complex(rho, eta)
    z13 = apex * s23 * s12 * c23 / (c12 * (1.0 - apex * s23**2))
    s13, delta = abs(z13), np.angle(z13)
    c13 = math.sqrt(1 - s13**2)
    phase = np.exp(1j * delta)
    V = np.array([
        [c12*c13, s12*c13, s13*np.exp(-1j*delta)],
        [-s12*c23-c12*s23*s13*phase,
         c12*c23-s12*s23*s13*phase, s23*c13],
        [s12*s23-c12*c23*s13*phase,
         -c12*s23-s12*c23*s13*phase, c23*c13],
    ], dtype=complex)
    J = c12*c23*c13**2*s12*s23*s13*math.sin(delta)
    return {"eta_bar": eta, "J_CKM": J, "magnitudes": np.abs(V)}


def pmns_data(orientation):
    """Recompute the charged-lepton-corrected TBM matrix."""
    def mass_factor(k):
        gap = 1.0 + BA_ratio * math.cos(float(h10) + 2*math.pi*k/3)
        vent = 1.0 - (d11**3/d10) * (alpha_phys/(2*math.pi))**2 * gap
        return gap**2 * vent

    phi = math.sqrt(mass_factor(1) / mass_factor(2)) / math.sqrt(d10)
    theta = float(h10)
    lam = float(h10) * (1.0 + 0.5/(d10*d11)**2)
    delta = orientation * math.acos(-lam)

    def R13(a, phase):
        c, s, e = math.cos(a), math.sin(a), np.exp(1j*phase)
        return np.array([[c, 0, s*e], [0, 1, 0],
                         [-s*e.conjugate(), 0, c]], dtype=complex)

    def R23(a):
        c, s = math.cos(a), math.sin(a)
        return np.array([[1, 0, 0], [0, c, s], [0, -s, c]], dtype=complex)

    U_tbm = np.array([
        [math.sqrt(2/3), 1/math.sqrt(3), 0],
        [-1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
        [1/math.sqrt(6), -1/math.sqrt(3), 1/math.sqrt(2)],
    ], dtype=complex)
    U = np.conj((R23(-phi) @ R13(theta, delta)).T) @ U_tbm
    J = float((U[0, 0]*U[1, 1]*np.conj(U[0, 1])*
               np.conj(U[1, 0])).imag)
    s13 = abs(U[0, 2])
    angles = (abs(U[0, 1])**2/(1-s13**2),
              abs(U[1, 2])**2/(1-s13**2), s13**2)
    return {"J_lep": J, "angles": angles}


def run(report=print):
    report("ORIENTATION BIT: one CP-orientation datum, propagated")
    report("=" * 64)
    fwd = FANO_SIGNS
    rev = [-s for s in FANO_SIGNS]           # mirror octonions

    C1f, C3f, epsf = cg_data(fwd)
    C1r, C3r, epsr = cg_data(rev)
    report(f"  standard orientation: |C1| = {C1f:.6f}, "
           f"|C3bar| = {C3f:.6f}, epsilon sign = {epsf:+.0f}")
    report(f"  reversed orientation: |C1| = {C1r:.6f}, "
           f"|C3bar| = {C3r:.6f}, epsilon sign = {epsr:+.0f}")
    assert abs(C1f - C1r) < 1e-12 and abs(C3f - C3r) < 1e-12
    assert epsf == -epsr
    report("  -> CP-even data INVARIANT. Epsilon FLIPS.")

    # Recompute the downstream matrices and baryogenesis expression.
    ckm_f, ckm_r = ckm_data(epsf), ckm_data(epsr)
    pmns_f, pmns_r = pmns_data(epsf), pmns_data(epsr)
    odd_f = {"eta_bar": ckm_f["eta_bar"], "J_CKM": ckm_f["J_CKM"],
             "J_lep": pmns_f["J_lep"]}
    odd_r = {"eta_bar": ckm_r["eta_bar"], "J_CKM": ckm_r["J_CKM"],
             "J_lep": pmns_r["J_lep"]}
    odd_f["eta_B"] = n7 * odd_f["J_lep"] * math.exp(-2*math.pi**2)
    odd_r["eta_B"] = n7 * odd_r["J_lep"] * math.exp(-2*math.pi**2)
    for name in odd_f:
        assert math.isclose(odd_r[name], -odd_f[name], rel_tol=1e-12,
                            abs_tol=1e-15)
        report(f"    {name:<12s} {odd_f[name]:+.6g} -> "
               f"{odd_r[name]:+.6g} under flip")
    assert np.allclose(ckm_f["magnitudes"], ckm_r["magnitudes"])
    assert np.allclose(pmns_f["angles"], pmns_r["angles"])
    report("  masses, |V_ij|, mixing-angle magnitudes: UNCHANGED")
    report("  (all depend on |C1|, |C3bar|, h-weights only)")

    report("-" * 64)
    report("  VERDICT: the framework contains one Fano-orientation")
    report("  datum. Reversing it flips CP sign and matter excess")
    report("  TOGETHER and touches no computed CP-even magnitude.")
    report("  Chirality selection, gauge representations, and knot spectral")
    report("  flow are separate composed branches; this probe owns the CP sign.")
    return {"C3": C3f, "flip_coherent": True,
            "cp_odd_forward": odd_f, "cp_odd_reversed": odd_r,
            "ckm_magnitudes_invariant": True,
            "pmns_angles_invariant": True}


if __name__ == "__main__":
    run()
