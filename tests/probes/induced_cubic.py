"""The induced action beyond quadratic order: the cubic term is Einstein's.

THE QUESTION. composite_pole.py and ledger_kernel.py compute the induced
kernel, the part of the action quadratic in the metric's departure from
flat. Newton's law and gamma = 1 depend on that part alone. The next
order, the three graviton vertex and with it the post-Newtonian beta,
depends on the cubic part. Is the two derivative term of the induced action the
Einstein-Hilbert density at that order, with the stiffness of the
quadratic part?

PREMISES. P1 a scalar row on the metric g, its operator
-(1/sqrt g) d(sqrt g g^-1 d) + xi R, under the proper time weight of
FACT 1: the action is -(1/2) int ds/s rho(s) Tr exp(-s D). P2 the
backgrounds
    ds^2 = dx^2 + exp(2 b1(x)) dy1^2 + exp(2 b2(x)) (dy2^2 + dy3^2),
x periodic, b1 and b2 trigonometric polynomials scaled by an amplitude
lam. They contain a trace part sigma = b1 + 2 b2 and a shear b1 - b2,
the polarization transverse to x.  At these amplitudes the metric's
components vary by a factor of five, so no expansion in the metric is
made. P3 the trace is computed from the operator's spectrum: Fourier
modes along y, a spectral grid in x, Gauss-Legendre nodes in the two
momenta.

STATEMENTS (the numbers below check each one).

  X1 (the expansion). (4 pi s)^2 Tr exp(-s D) per unit coordinate
     volume in y is c0 + c1 s + ..., with c0 = int sqrt g dx and
     c1 = (1/6) int sqrt g R dx = (1/6) int exp(sigma) (sigma'^2 -
     b1'^2 - 2 b2'^2) dx: the fitted coefficients agree with these to
     1e-6 at lam = 1 and at lam = -1.
  X2 (the cubic term). c1 is not even in lam: its odd part, the terms
     of third and higher odd order in the metric, is 11 percent of it
     at lam = 1.  The trace reproduces that part to 1e-5. The
     quadratic kernel alone leaves this part free. The trace fixes it
     at the Einstein value, with the 1/6 of the quadratic part.
  X3 (the coupling xi). With xi R in the operator c1 is
     (1/6 - xi) int sqrt g R dx at the same accuracy: xi multiplies
     the whole Einstein density, cubic term included.
  X4 (a row with a speed of its own). A row whose kinetic form carries
     a frame term, its x-derivative scaled by c^2 = 1.3
     (saddle_stabilizer.py Z7), induces the Einstein density of its
     own metric: c1/c0 is c^2 times the common metric's value. Two such
     rows give two metrics, and no single cubic vertex.

RESULT. For a row that couples to the one metric, the two derivative
term of the induced action is (1/6 - xi) sqrt g R beyond quadratic
order, on backgrounds that are not small.  The cubic term is Einstein's
with the stiffness of the kernel. With gamma = 1 (composite_pole.py K5)
the action is Einstein's up to the order on which the post-Newtonian
beta depends.
The same vertex is found from the diagrams, with one Planck mass for
species of different mass (Carone, Claringbold, Vaman
arXiv:1710.09367).

Scope: scalar rows. Backgrounds depending on one coordinate. The term
with two derivatives. The fermion and vector rows are computed on the
same backgrounds in induced_cubic_rows.py.

Usage: python3 tests/probes/induced_cubic.py     (about ten seconds)
"""

import math

import numpy as np

N = 64                                           # spectral grid in x
NK, KMAX = 120, 46.0                             # Gauss-Legendre nodes in each momentum
S_LO, S_HI, S_PTS, DEGREE = 0.03, 0.15, 14, 7    # proper times of the fit
X = 2 * np.pi * np.arange(N) / N
FREQ = np.fft.fftfreq(N, 1 / N)
D2 = np.real(np.fft.ifft(-(FREQ**2)[:, None] * np.fft.fft(np.eye(N), axis=0), axis=0))
DX = 2 * np.pi / N


def _d(f):
    return np.real(np.fft.ifft(1j * FREQ * np.fft.fft(f)))


def background(lam):
    """b1, b2 and the geometry they fix: sigma, its derivatives, the scalar curvature."""
    b1 = lam * (0.30 * np.cos(X) + 0.20 * np.sin(2 * X))
    b2 = lam * (0.25 * np.sin(X) - 0.15 * np.cos(2 * X) + 0.10 * np.cos(X))
    p1, p2 = _d(b1), _d(b2)
    sig, sp = b1 + 2 * b2, p1 + 2 * p2
    spp = _d(sp)
    R = -2 * spp - (p1**2 + 2 * p2**2) - sp**2
    return dict(b1=b1, b2=b2, sig=sig, sp=sp, spp=spp, shear2=p1**2 + 2 * p2**2, R=R)


def einstein(lam, xi=0.0):
    """c0 = int sqrt g and c1 = (1/6 - xi) int sqrt g R, the second after one integration by parts."""
    g = background(lam)
    c0 = float(np.sum(np.exp(g["sig"])) * DX)
    direct = float(np.sum(np.exp(g["sig"]) * g["R"]) * DX)
    parts = float(np.sum(np.exp(g["sig"]) * (g["sp"]**2 - g["shear2"])) * DX)
    assert abs(direct - parts) < 1e-10 * max(1.0, abs(parts))
    return c0, (1 / 6 - xi) * parts


def trace(lam, s, xi=0.0, c2=1.0):
    """Tr exp(-s D) per unit coordinate volume in y, D symmetrized by u = g^(1/4) f; c2 scales the x-derivative."""
    g = background(lam)
    W = g["spp"] / 2 + g["sp"]**2 / 4                      # (g^(1/4))'' / g^(1/4)
    A = c2 * (-D2 + np.diag(W)) + np.diag(xi * g["R"])
    node, wt = np.polynomial.legendre.leggauss(NK)
    k, wk = 0.5 * KMAX * (node + 1), 0.5 * KMAX * wt
    P1, P2 = np.exp(-2 * g["b1"]), np.exp(-2 * g["b2"])
    eye = np.eye(N)[None]
    out = np.zeros(len(s))
    for k1, w1 in zip(k, wk):
        H = A[None] + (k1 * k1 * P1 + k[:, None]**2 * P2[None])[:, :, None] * eye
        ev = np.linalg.eigvalsh(H)                         # one spectrum per transverse momentum
        out += w1 * np.array([np.sum(wk * k * np.exp(-si * ev).sum(axis=1)) for si in s])
    return out / (2 * np.pi**2)                            # int d^3k/(2 pi)^3, axial symmetry


def coefficients(lam, xi=0.0, c2=1.0):
    """c0 and c1 of (4 pi s)^2 Tr exp(-s D), from a polynomial fit on Chebyshev proper times."""
    cheb = np.cos(np.pi * (np.arange(S_PTS) + 0.5) / S_PTS)
    s = 0.5 * (S_HI + S_LO) + 0.5 * (S_HI - S_LO) * cheb
    F = (4 * np.pi * s)**2 * trace(lam, s, xi, c2)
    c = np.polynomial.polynomial.polyfit(s, F, DEGREE)
    return float(c[0]), float(c[1])


def run(report=print):
    out = {}
    spread = float(np.exp(2 * np.ptp(background(1.0)["b1"])))
    assert spread > 2.5                                    # the metric is not near flat

    # X1: the first two coefficients on both signs of the amplitude.
    fit = {lam: coefficients(lam) for lam in (1.0, -1.0)}
    exact = {lam: einstein(lam) for lam in (1.0, -1.0)}
    err = max(abs(fit[lam][i] / exact[lam][i] - 1) for lam in fit for i in (0, 1))
    assert err < 1e-6
    out.update(c1={lam: fit[lam][1] for lam in fit}, einstein_c1={lam: exact[lam][1] for lam in exact},
               expansion_error=err, metric_spread=spread)
    report(f"X1 expansion: c1 = {fit[1.0][1]:.7f} and {fit[-1.0][1]:.7f} at lam = +1, -1 against "
           f"(1/6) int sqrt g R = {exact[1.0][1]:.7f} and {exact[-1.0][1]:.7f}; c0 and c1 agree to {err:.0e}; "
           f"g_11 varies by a factor {spread:.1f}")

    # X2: the part odd in the amplitude.
    odd_fit = 0.5 * (fit[1.0][1] - fit[-1.0][1])
    odd_exact = 0.5 * (exact[1.0][1] - exact[-1.0][1])
    share = odd_exact / exact[1.0][1]
    assert share > 0.1 and abs(odd_fit / odd_exact - 1) < 1e-5
    out.update(odd_share=share, odd_error=abs(odd_fit / odd_exact - 1))
    report(f"X2 cubic term: the part of c1 odd in the amplitude is {100 * share:.1f}% of it, and the trace gives "
           f"{odd_fit:.8f} against Einstein's {odd_exact:.8f} ({out['odd_error']:.0e})")

    # X3: the coupling xi multiplies the whole density.
    xi = 0.1
    fx, ex = coefficients(1.0, xi=xi), einstein(1.0, xi=xi)
    assert abs(fx[1] / ex[1] - 1) < 1e-6 and abs(ex[1] / exact[1.0][1] - (1 - 6 * xi)) < 1e-12
    out.update(xi_error=abs(fx[1] / ex[1] - 1))
    report(f"X3 with xi = {xi}: c1 = {fx[1]:.7f} against (1/6 - xi) int sqrt g R = {ex[1]:.7f}")

    # X4: a row with a speed of its own induces the Einstein density of its own metric.
    c2 = 1.3
    fc = coefficients(1.0, c2=c2)
    own = (fc[1] / fc[0]) / (exact[1.0][1] / exact[1.0][0])
    assert abs(own / c2 - 1) < 1e-6 and abs(fc[0] * math.sqrt(c2) / exact[1.0][0] - 1) < 1e-6
    out.update(frame_ratio=own)
    report(f"X4 a row with c^2 = {c2}: c1/c0 is {own:.6f} times the common metric's value")
    report("RESULT: the two derivative term of the induced action is (1/6 - xi) sqrt g R beyond quadratic order, "
           "the cubic term with the kernel's stiffness. A row with a speed of its own induces the Einstein density "
           "of its own metric.")
    return out


if __name__ == "__main__":
    run()
