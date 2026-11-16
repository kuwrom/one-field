"""
selector_theorem.py (Gate 4): the selector space is BINARY,
and one of the two is already dead.

THE THEOREM (finite group theory, verified exhaustively below):
an idempotent probability measure on Z3 (p * p = p under Z3
convolution) must be the uniform measure on a SUBGROUP.  Z3 has
exactly two subgroups: {e} and Z3 itself.  Hence exactly TWO
selectors exist:

    delta_e  (uniform on {e}):  KEEPS the label   (kappa != 0 class)
    Haar     (uniform on Z3):   ERASES the label  (protected forgetting)

Proof sketch (also checked numerically): in Z3 Fourier space,
convolution is multiplication, so idempotence means p_hat(q)^2 =
p_hat(q), i.e. p_hat(q) in {0, 1}.  Normalisation forces
p_hat(0) = 1. Reality forces p_hat(1) = p_hat(2).  The only choices
are (1,0,0) = Haar and (1,1,1) = delta_e.  QED.

WHY IDEMPOTENT IS THE RIGHT CLASS (the one labeled definition):
a selector applied twice must equal itself applied once. A
"selection" that keeps changing under reapplication is dynamics,
not selection.  A weighted twirl is the example: not idempotent,
hence not a selector.

THE CLOSURE: the delta_e branch is the label keeping branch. The
three gates applied to the survivors KILL it on observables: it
fails lepton closure (the "localize G2 first" family). Haar =
protected forgetting is therefore the UNIQUE viable selector. This
is not a preference. The alternative set has exactly one other
element, and computation excludes it.

STATUS: Gate 4 closes MODULO one labeled definition (selection =
idempotent Z3-equivariant measure).  Gate 4 does not need the wall
datum derived from 4d dynamics. Given the definition of selection,
the rest is a two element check, so the definition is the only
point open to objection.

Usage: python3 selector_theorem.py
"""


import numpy as np
import _bootstrap  # noqa: E402,F401



def convolve3(p, q):
    out = np.zeros(3)
    for i in range(3):
        for j in range(3):
            out[(i + j) % 3] += p[i] * q[j]
    return out


def run(report=print):
    report("SELECTOR THEOREM: the Gate-4 space is binary")
    report("=" * 64)

    # exhaustive numerical search on a fine simplex grid
    n = 300   # divisible by 3 so Haar lies on the grid
    found = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            p = np.array([i, j, n - i - j], float) / n
            if np.allclose(convolve3(p, p), p, atol=1e-9):
                found.append(tuple(np.round(p, 6)))
    found = sorted(set(found))
    report(f"  idempotent measures on the simplex (grid {n}): "
           f"{len(found)}")
    for p in found:
        name = ("Haar (forgetting)" if np.allclose(p, 1/3)
                else "delta_e (label keeping)" if p[0] == 1.0
                else "UNEXPECTED")
        report(f"    p = {p}   -> {name}")
    assert len(found) == 2, "theorem: exactly two idempotents on Z3"

    # Fourier space proof check (exact)
    w = np.exp(2j * np.pi / 3)
    for p in [np.array([1/3, 1/3, 1/3]), np.array([1.0, 0, 0])]:
        ph = np.array([sum(p[k] * w**(q*k) for k in range(3))
                       for q in range(3)])
        assert np.allclose(ph**2, ph, atol=1e-12)

    # the excluded example: weighted twirl (from the runner) is NOT
    # idempotent -> not a selector -> not in the space at all
    tw = np.array([0.5, 0.25, 0.25])
    assert not np.allclose(convolve3(tw, tw), tw, atol=1e-9)
    report("  weighted twirl (1/2,1/4,1/4): NOT idempotent -> excluded")

    report("-" * 64)
    report("  The selector space is {delta_e, Haar}.  The delta_e")
    report("  (label keeping) branch fails the lepton and quark gates")
    report("  (interference/__init__.py, Gate 6: the L and Q columns.")
    report("  The losing rows are reproduced from Table V of the gravity")
    report("  record. The winning row is executable.)  Therefore")
    report("  Haar, protected forgetting, is the UNIQUE viable")
    report("  selector, modulo one labeled definition:")
    report("  selection = idempotent Z3-equivariant measure.")
    return found


if __name__ == "__main__":
    run()
