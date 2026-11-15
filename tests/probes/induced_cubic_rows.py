"""The cubic term for the Weyl and vector rows.

THE QUESTION. induced_cubic.py finds the scalar row's two derivative
term to be (1/6 - xi) sqrt g R beyond quadratic order. The fermion and
vector rows enter the ledger with weights 1/12 and 1/6 - c/2 per state,
read from their heat traces at quadratic order (ledger_kernel.py). On
backgrounds that are not small, do their heat traces carry the Einstein
density with those weights?

PREMISES. P1 the backgrounds of induced_cubic.py P2,
ds^2 = dx^2 + exp(2 b1) dy1^2 + exp(2 b2)(dy2^2 + dy3^2), with the
orthonormal frame e^i = exp(b_i) dy_i: the only connection components
are the rotations in the (i, x) planes, at rate b_i' exp(b_i) along
y_i. P2 for a field in a representation with generators T_i of those
rotations, the Laplace type operator on a Fourier mode k along y is
    D = [-(1/sqrt g) d sqrt g d + sum k_i^2 exp(-2 b_i)]
        - 2i sum k_i b_i' exp(-b_i) T_i - sum b_i'^2 T_i^2 + E,
with T_i = gamma_i gamma_x / 2 and E = R/4 for the squared Dirac
operator (Lichnerowicz), and T_i the vector rotation and E the Ricci
tensor, R_xx = -sum (b_i'' + b_i'^2), R_ii = -b_i'' - b_i' sigma', for
the vector's Hodge operator. P3 the spectral method of
induced_cubic.py, on a coarser grid.

STATEMENTS (the numbers below check each one).

  Y1 (the squared Dirac operator). (4 pi s)^2 Tr exp(-s D) is
     4 c0 + c1 s + ... with c1 = (4/6 - 1) int sqrt g R: the fitted
     value agrees to 1e-6 at lam = +1 and -1, the odd part included.
     With the fermion loop's sign and the half of the square root the
     Dirac field weighs 1/3, four states at 1/12.
  Y2 (the vector). The Hodge operator's trace has
     c1 = (4/6 - 1) int sqrt g R at the same accuracy. Less two ghost
     scalars the gauge field weighs -2/3, two states at 1/6 - 1/2, the
     covariant contact.

RESULT. The Weyl and vector rows induce the Einstein density beyond
quadratic order with the weights the ledger uses.  The cubic term of
the whole list is Einstein's with the stiffness Sigma.

Scope: the term with two derivatives. Backgrounds depending on one
coordinate. The covariant contact for the vector.

Usage: python3 tests/probes/induced_cubic_rows.py     (about fifteen seconds)
"""

import numpy as np

N = 48                                           # spectral grid in x
NK, KMAX = 40, 38.0                              # Gauss-Legendre nodes in each momentum
S_LO, S_HI, S_PTS, DEGREE = 0.05, 0.17, 14, 7    # proper times the grid resolves
X = 2 * np.pi * np.arange(N) / N
FREQ = np.fft.fftfreq(N, 1 / N)
D2 = np.real(np.fft.ifft(-(FREQ**2)[:, None] * np.fft.fft(np.eye(N), axis=0), axis=0))
DX = 2 * np.pi / N


def _d(f):
    return np.real(np.fft.ifft(1j * FREQ * np.fft.fft(f)))


def background(lam):
    b1 = lam * (0.30 * np.cos(X) + 0.20 * np.sin(2 * X))
    b2 = lam * (0.25 * np.sin(X) - 0.15 * np.cos(2 * X) + 0.10 * np.cos(X))
    p1, p2 = _d(b1), _d(b2)
    sp = p1 + 2 * p2
    return dict(b1=b1, b2=b2, p1=p1, p2=p2, q1=_d(p1), q2=_d(p2), sig=b1 + 2 * b2, sp=sp, spp=_d(sp))


def einstein(lam):
    """int sqrt g dx and int sqrt g R dx."""
    g = background(lam)
    w = np.exp(g["sig"]) * DX
    return float(np.sum(w)), float(np.sum(w * (g["sp"]**2 - g["p1"]**2 - 2 * g["p2"]**2)))


def generators(rep):
    """T_1, T_2, T_3^2 and the multiplicity: the momentum lies in the (y1, y2) plane."""
    if rep == "spinor":                          # gamma_i gamma_x / 2 on a two component block, twice
        sz, sx = np.diag([1.0 + 0j, -1.0]), np.array([[0, 1], [1, 0]], dtype=complex)
        return 0.5j * sz, 0.5j * sx, -0.25 * np.eye(2), 2.0
    T = []
    for i in (1, 2, 3):                          # frame order (x, 1, 2, 3)
        t = np.zeros((4, 4), dtype=complex)
        t[i, 0], t[0, i] = 1.0, -1.0
        T.append(t)
    return T[0], T[1], T[2] @ T[2], 1.0


def trace(lam, s, rep):
    """Tr exp(-s D) per unit coordinate volume in y for the row's operator of P2."""
    g = background(lam)
    T1, T2, T3sq, mult = generators(rep)
    dim = T1.shape[0]
    I = np.eye(dim)
    W = g["spp"] / 2 + g["sp"]**2 / 4
    base = np.kron(-D2 + np.diag(W), I).astype(complex)
    base -= (np.kron(np.diag(g["p1"]**2), T1 @ T1) + np.kron(np.diag(g["p2"]**2), T2 @ T2)
             + np.kron(np.diag(g["p2"]**2), T3sq))
    if rep == "spinor":
        R = -2 * g["spp"] - (g["p1"]**2 + 2 * g["p2"]**2) - g["sp"]**2
        base += np.kron(np.diag(R / 4), I)
    else:
        ric = np.zeros((N, 4))
        ric[:, 0] = -(g["q1"] + 2 * g["q2"]) - (g["p1"]**2 + 2 * g["p2"]**2)
        ric[:, 1] = -g["q1"] - g["p1"] * g["sp"]
        ric[:, 2] = ric[:, 3] = -g["q2"] - g["p2"] * g["sp"]
        base += np.diag(ric.ravel())
    lin1 = np.kron(np.diag(g["p1"] * np.exp(-g["b1"])), -2j * T1)
    lin2 = np.kron(np.diag(g["p2"] * np.exp(-g["b2"])), -2j * T2)
    K1, K2 = np.kron(np.diag(np.exp(-2 * g["b1"])), I), np.kron(np.diag(np.exp(-2 * g["b2"])), I)
    node, wt = np.polynomial.legendre.leggauss(NK)
    k, wk = 0.5 * KMAX * (node + 1), 0.5 * KMAX * wt
    out = np.zeros(len(s))
    for k1, w1 in zip(k, wk):
        H = (base + k1 * k1 * K1 + k1 * lin1)[None] + (k**2)[:, None, None] * K2[None] + k[:, None, None] * lin2[None]
        ev = np.linalg.eigvalsh(H)
        out += w1 * np.array([np.sum(wk * k * np.exp(-si * ev).sum(axis=1)) for si in s])
    return mult * out / (2 * np.pi**2)


def coefficients(lam, rep):
    cheb = np.cos(np.pi * (np.arange(S_PTS) + 0.5) / S_PTS)
    s = 0.5 * (S_HI + S_LO) + 0.5 * (S_HI - S_LO) * cheb
    c = np.polynomial.polynomial.polyfit(s, (4 * np.pi * s)**2 * trace(lam, s, rep), DEGREE)
    return float(c[0]), float(c[1])


def run(report=print):
    out = {}
    exact = {lam: einstein(lam) for lam in (1.0, -1.0)}
    for rep, label in (("spinor", "Y1 squared Dirac operator"), ("vector", "Y2 vector Hodge operator")):
        fit = {lam: coefficients(lam, rep) for lam in (1.0, -1.0)}
        err0 = max(abs(fit[lam][0] / (4 * exact[lam][0]) - 1) for lam in fit)
        err1 = max(abs(fit[lam][1] / ((4 / 6 - 1) * exact[lam][1]) - 1) for lam in fit)
        odd = 0.5 * (fit[1.0][1] - fit[-1.0][1]) / ((4 / 6 - 1) * 0.5 * (exact[1.0][1] - exact[-1.0][1])) - 1
        assert err0 < 1e-8 and err1 < 1e-6 and abs(odd) < 1e-5
        out[rep] = dict(c1=fit[1.0][1], error=err1, odd_error=abs(odd))
        report(f"{label}: c1 = {fit[1.0][1]:.7f} against (4/6 - 1) int sqrt g R = "
               f"{(4 / 6 - 1) * exact[1.0][1]:.7f} ({err1:.0e}. Odd part {abs(odd):.0e})")
    dirac = -(4 / 6 - 1)                              # fermion loop sign on the squared operator
    gauge = (4 / 6 - 1) - 2 / 6                       # less two ghost scalars
    assert abs(dirac - 4 / 12) < 1e-12 and abs(gauge - 2 * (1 / 6 - 1 / 2)) < 1e-12
    out.update(dirac_weight=dirac, gauge_weight=gauge)
    report(f"weights: Dirac field {dirac:.4f} = 4 x 1/12. Gauge field {gauge:.4f} = 2 x (1/6 - 1/2)")
    report("RESULT: the Weyl and vector rows induce the Einstein density beyond quadratic order with the ledger's "
           "weights.")
    return out


if __name__ == "__main__":
    run()
