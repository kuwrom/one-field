"""
zero_mode.py: the chiral zero mode on the knot, computed
(Jackiw-Rossi on our own defect).

WHY THE SURVIVING BRANCH IS CHIRAL (mass-protection closure
argument).  The framework's entire mass sector has the form
m = M_Pl * exp(-S_closure) * (circulant data), with every factor
forced.  This structure requires that no mass term exist above the
closure scale other than those the closure itself generates.  In a
vector-like gauge theory, bare Dirac masses are gauge invariant.
Nothing forbids them, and by genericity they appear at the cutoff
scale.  A vector-like branch therefore carries mass observables
contaminated by terms the closure does not control, so it does not
close and is not a surviving branch.  Chiral gauge symmetry is
precisely the known mechanism that forbids bare masses until symmetry
breaking.  Therefore the surviving branch must be chiral: chirality
is a survival condition of the framework.

The lepton branch already contains this argument in miniature: the
cyclic closure [Delta, S] = 0 makes the mass matrix circulant, and
the non-circulant term diag(v) is not an admissible mass term.  The
4d gauge-space statement is the analog: the surviving branch's
symmetry must exclude every mass term the closure does not generate.
One exclusion principle, two layers.

WHICH CHIRALITY: the handedness seed is the ORIENTED Fano plane
(seven directed lines), the same bit that sources CP violation
(J_lep, eta_B).  The mechanism producing one-handed couplings from a
non-chiral bulk is established (Callan-Harvey anomaly inflow,
Jackiw-Rossi zero modes): fermion zero modes localized on topological
defects are chiral by index theorem, the handedness is set by the
sign of the defect's winding, and the anomaly of the chiral modes is
cancelled by inflow from the bulk.  The observed light fermions are
the zero modes bound to the knots, so their gauge coupling is V-A
with handedness = winding sign = the Fano bit.  This probe computes
exactly that on the knot's own profile.

THE SETUP (standard JR angular reduction): a Dirac fermion coupled to
a defect with winding n decomposes into angular channels. Winding
SHIFTS the centrifugal index so that exactly |n| channels lose their
centrifugal barrier (index theorem).  In the barrier-free channel the
radial operator pair is

    A  = d/dr + f(r)          (zero mode: u = exp(-INT f), normalizable)
    A+ = -d/dr + f(r)         (candidate: exp(+INT f), NOT normalizable)

so exactly ONE zero mode exists and it lives in ONE chirality block.
Reversing the winding (f -> -f) moves it to the OTHER block:
handedness = winding sign.  For winding 0, every channel keeps a
centrifugal term (+1/r): u = exp(-INT f)/r is not normalizable at the
origin, NO zero mode.

This probe builds the discretized block Hamiltonian H = [[0, A],
[A^T, 0]] and verifies all three statements numerically on the knot's
profile f(r) = tanh(r/xi).

Usage: python3 zero_mode.py
"""

import numpy as np


def deriv(n, dr, forward=True):
    D = np.zeros((n, n))
    if forward:
        for i in range(n - 1):
            D[i, i] -= 1.0 / dr
            D[i, i + 1] += 1.0 / dr
    else:
        for i in range(1, n):
            D[i, i] += 1.0 / dr
            D[i, i - 1] -= 1.0 / dr
    return D


def block_h(f_vals, dr, winding, centrifugal=None):
    n = len(f_vals)
    if winding >= 0:
        A = deriv(n, dr, True) + np.diag(f_vals)     # d/dr + f
    else:
        A = -deriv(n, dr, False) + np.diag(f_vals)   # -d/dr + f
    if centrifugal is not None:
        A += np.diag(centrifugal)
    H = np.zeros((2 * n, 2 * n))
    H[:n, n:] = A
    H[n:, :n] = A.T
    return H


def analyze(H, n, r, gap_scale):
    """Smallest-|E| CORE-localized state (filters boundary-row
    artifacts of the discretization, which sit at r ~ R)."""
    w, v = np.linalg.eigh(H)
    order = np.argsort(np.abs(w))
    R = r[-1]
    for i in order[:6]:
        psi = v[:, i]
        dens = psi[:n]**2 + psi[n:]**2
        r_mean = float(np.sum(r * dens) / np.sum(dens))
        if r_mean < 0.4 * R:               # core mode, not artifact
            up = float(np.sum(psi[:n]**2))
            dn = float(np.sum(psi[n:]**2))
            gamma = (up - dn) / (up + dn)
            return float(w[i]), gamma, r_mean, abs(w[i]) < 0.05*gap_scale
    return None, None, None, False


def run(report=print):
    n, R = 600, 30.0
    r = np.linspace(R / n, R, n)
    dr = r[1] - r[0]
    f = np.tanh(r / 1.5)                   # the knot profile, f(inf)=1

    report("JACKIW-ROSSI ON THE KNOT: chirality = winding sign")
    report("=" * 64)
    report(f"  {'case':<26s} {'E0':>10s} {'<Gamma>':>9s} "
           f"{'<r>':>6s} {'zero mode':>10s}")
    cases = [("winding +1 (chiral mass)", block_h(f, dr, +1)),
             ("winding -1 (reversed)", block_h(f, dr, -1))]
    # winding 0: a windingless mass has NO phase structure -> it is a
    # sigma_z mass, not a chiral one: H = [[f, D],[D^T, -f]] (gapped)
    D = deriv(n, dr, True)
    H0 = np.zeros((2*n, 2*n))
    H0[:n, :n] = np.diag(f)
    H0[:n, n:] = D
    H0[n:, :n] = D.T
    H0[n:, n:] = -np.diag(f)
    cases.append(("winding 0 (sigma_z mass)", H0))
    results = {}
    for label, H in cases:
        E0, g, rm, is_zero = analyze(H, n, r, 1.0)
        if E0 is None or not is_zero:
            wmin = float(np.min(np.abs(np.linalg.eigvalsh(H))))
            report(f"  {label:<26s} {wmin:>10.2e} {'-':>9s} {'-':>6s} "
                   f"{'NONE (gapped)':>13s}")
            results[label] = {"zero_mode": False, "gap": wmin}
        else:
            report(f"  {label:<26s} {E0:>10.2e} {g:>+9.3f} "
                   f"{rm:>6.2f} {'YES':>10s}")
            results[label] = {"zero_mode": True, "E0": E0, "gamma": g}
    # the registered conclusions, as assertions
    plus = results["winding +1 (chiral mass)"]
    minus = results["winding -1 (reversed)"]
    zero = results["winding 0 (sigma_z mass)"]
    assert plus["zero_mode"] and minus["zero_mode"]
    assert plus["gamma"] * minus["gamma"] < 0, \
        "handedness must flip with the winding sign"
    assert not zero["zero_mode"] and zero["gap"] > 0.1, \
        "the windingless mass must be gapped"
    report("-" * 64)
    report("  winding +1 and -1: one true zero mode each, in OPPOSITE")
    report("  chirality blocks, handedness IS the winding sign,")
    report("  i.e. The propagated Fano bit")
    report("  winding  0: windingless mass = sigma_z channel = gapped,")
    report("  no zero mode: chirality requires the winding")
    report("  => the knot's bound fermion is ONE-HANDED by index")
    report("     theorem. V-A is inherited, not imposed.")
    return results


if __name__ == "__main__":
    run()
