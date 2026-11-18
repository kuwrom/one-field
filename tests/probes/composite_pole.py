"""The massless pole of the composite metric, computed for the carrier.

composite_metric.py shows that the carrier written without a background
metric is the theory of an actionless metric g = h/U. This probe
computes what that metric does at one loop: the full tensor kernel of
the induced quadratic action, in all ten components, at every momentum,
for the carrier's real scalar components in the empty field.

PREMISES. P1 the empty field: free components of mass m (the carrier's
interactions enter at higher order in the fields). P2 Euclidean D = 4,
one loop of the components, no loop of the composite metric itself:
leading order in the number of components. P3 a generally covariant
regulator, written as a weight rho(t) on the proper time of the heat
trace, W = -(1/2) int dt/t rho(t) Tr exp(-t D_g). Two weights are run:
Pauli-Villars, rho = sum_j c_j exp(-t s_j), and the proper time cutoff
rho = theta(t - t0), the regulator of sakharov_normalisation.py FACT 1.

THE KERNEL. With g = 1 + h, S = int sqrt g [(1/2) g^{mu nu} d phi d phi
+ (1/2) M^2 phi^2], the quadratic action is
    Gamma2[h; q] = <S2> - (1/2) <S1 S1>
      = (I1/4)(h.h/2 - (tr h)^2/4) - (1/4) int_p (h:V)^2 G(p) G(p+q),
    h:V = h_{mu nu} p_mu (p+q)_nu - (1/2) tr h [p.(p+q) + M^2].
A Feynman parameter x and the shift l = p + x q reduce it to
    A = h - (1/2) tr h,  c0 = -x(1-x) qAq,  X = x(1-x) q^2,
    int_0^1 dx [ 6 T4 J(2,0) + 2 (b^2/4 + c0 trA/2) J(1,0)
                 - (1/2) tr h trA J(1,1) + c0^2 J(0,0)
                 - c0 tr h J(0,1) + (1/4)(tr h)^2 J(0,2) ] / (16 pi^2),
    T4 = ((trA)^2 + 2 tr A^2)/24,  b^2 = (1-2x)^2 |Aq|^2,
    J(a,n;X) = int dt rho(t) e^{-t m^2} d_t^n [t^{-a-1} e^{-tX}],
and I1 = 2 J(2,0;0)/(16 pi^2). The form is linear in the weight, so it
holds for every regulator of P3.

STATEMENTS (the numbers below check each one, for both weights).

  K1 (vacuum term). Gamma2[h; 0] = rho_vac (-h.h/4 + (tr h)^2/8), the
     second order part of rho_vac sqrt g. With the potential constant
     V0 = -rho_vac the q = 0 term cancels in every component. Off that
     point it is a mass term of the form h.h - (tr h)^2/2, not
     Pauli-Fierz (Carone arXiv:1911.13166 Eq. 3.5).
  K2 (vacuum metric). <d_mu phi d_nu phi> = (I1/4) delta and
     I1/4 = -rho_vac, so <g> = <d phi d phi>/V0 = delta at that same
     point. The flat metric is the fluctuation term. With the proper
     time weight, I1 > 0 and V0 > 0: a positive regulator gives a
     positive Euclidean metric. The fluctuations' vacuum term and the
     flat metric are one number, I1/4 = -rho_vac. A state with no
     fluctuation carries no metric. With V0 = 0 the elimination g = h/U
     has no vacuum at all. So in the full form flatness fixes the
     potential's constant against that number, and a nonzero
     Lambda = 8 pi G (V0 + rho_vac) is a constant of the action. In the
     traceless form Lambda is a constant of the state, and the form
     keeps a fixed volume element (|det h|^(1/4) has density weight one
     half).
  K3 (protection). K(q) = Gamma2(q) - Gamma2(0) annihilates every
     h = q xi + xi q, at every q.
  K4 (two form factors). K(q) = F2(q^2) P2 + F0(q^2) P0s, only the
     transverse spin-2 and spin-0 projectors.
  K5 (Einstein). As q -> 0, 4 F2/q^2 -> (1/6) x the regulator's
     stiffness integral, the heat kernel weight a1 = 1/6 per real
     component, and F0/F2 -> -2: the linearized Einstein-Hilbert
     operator. Between conserved sources its inverse is
     (T.T - T^2/2)/q^2: the massless pole, Newton's law with
     1/G = 16 pi N (stiffness), and gamma = (1 - r)/(2 + r) -> 1,
     r = F2/F0. Under the proper time cutoff 1/G = (N/6)/(2 pi t0)
     E2(m^2 t0), FACT 1's 1/G = Sigma/(2 pi t0).
  K6 (nothing else below the regulator). F2/q^2 > 0 > F0/q^2 on the
     whole Euclidean axis. Below the two particle threshold the only
     zero is q^2 = 0 once the regulator scale is ten times m. The
     departure of F0/F2 from -2 at |q| = m falls with the regulator
     scale: it drops by a factor 45 from L = 10 m to L = 100 m.
  K7 (no added state). Continued to timelike q^2 = -s, the kernel
     under the proper time cutoff has an absorptive part only above
     the two particle threshold, where |Im F2| = s^2 beta^5/(7680 pi),
     beta = sqrt(1 - 4 m^2/s), for every t0: the cut of the components
     alone. The cutoff adds no threshold of its own. The Pauli-Villars
     weight adds its own negative norm thresholds at 4 M_j^2.
  K8 (traceless form). On traceless h (det g = 1) the same kernel has
     the same five spin-2 eigenvalues, and the static exchange between
     conserved sources is Newton's law with the same G.
  K9 (locality fixes the ratio). With F2 = a q^2 and F0 = b q^2 the
     kernel q^2 (a P2 + b P0s) carries the term (q.h.q)^2/q^2 with
     coefficient (2a + b)/3, quartic in the direction of q. Every local
     two derivative term is at most quadratic in it. The q^2 term of a
     generally covariant kernel analytic at q = 0 is local, so
     F0/F2 -> -2 at every order of the loop expansion: higher orders
     move the stiffness, not the structure. The check fits the
     kernel over directions of q at fixed |q| and measures its
     degree four content, against a kernel shifted off the ratio.
     The kernel is induced under a covariant regulator, so it is
     covariant for any form of the vacuum term, and the check
     holds on traceless h as well: the traceless form the lane reads
     admits no other ratio.

RESULT. For the carrier in the empty field at one loop, under both
regulator weights, the kernel carries only the spin-2 and spin-0 form
factors and annihilates every gauge direction. At low momentum the
stiffness is the heat kernel weight a1 = 1/6 per real component and
F0/F2 = -2, so the static exchange is Newton's law with gamma = 1 in
the full and the traceless form. A q^2 term off F0/F2 = -2 is
nonlocal and is excluded.

Usage: python3 tests/probes/composite_pole.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np
from scipy.special import exp1, expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import carrier_sakharov as cs  # noqa: E402

# Feynman parameter nodes on [0, 1], cosine mapped: the integrands peak at
# the endpoints at large q^2 t and are logarithmic there for massless rows.
_U, _W = np.polynomial.legendre.leggauss(200)
_U, _W = 0.5 * (_U + 1), 0.5 * _W
XG = 0.5 * (1 - np.cos(math.pi * _U))
WG = _W * 0.5 * math.pi * np.sin(math.pi * _U)
C16 = 1.0 / (16 * math.pi**2)
PV4 = ((1.0, 0.0), (-3.0, 1.0), (3.0, 2.0), (-1.0, 3.0))   # (c_j, s_j / L^2)
FACTORIAL = (1.0, 1.0, 2.0)
N_COMPONENTS = 6                                         # three complex fields


def pauli_villars(m2, lam):
    """J(a, n; X) for rho = sum c exp(-t s): (-1)^(a+1) sum c M^2n D^a ln D / a!."""
    M2 = [(c, m2 + s * lam * lam) for c, s in PV4]

    def J(a, n, X):
        return (-1)**(a + 1) * sum(c * M**n * (M + X)**a * np.log(M + X)
                                   for c, M in M2) / FACTORIAL[a]
    rho = sum(c * M * M * math.log(M) for c, M in M2) / (64 * math.pi**2)
    stiffness = sum(c * M * math.log(M) for c, M in M2) / (192 * math.pi**2)
    return J, rho, stiffness


def proper_time(m2, t0):
    """J(a, n; X) for rho = theta(t - t0)."""
    def J(a, n, X):
        D = X + m2

        def E(k):
            return t0**(1 - k) * expn(k, t0 * D)
        if n == 0:
            return E(a + 1)
        if n == 1:
            return -(a + 1) * E(a + 2) - X * E(a + 1)
        return (a + 1) * (a + 2) * E(a + 3) + 2 * (a + 1) * X * E(a + 2) + X * X * E(a + 1)
    rho = -expn(3, t0 * m2) / (32 * math.pi**2 * t0**2)
    return J, rho, expn(2, t0 * m2) / (192 * math.pi**2 * t0)


def gamma2(h, q, J):
    """Induced quadratic action per real component; q may be timelike (complex)."""
    h = 0.5 * (h + h.T)
    trh = np.trace(h)
    A = h - 0.5 * trh * np.eye(4)
    Aq = A @ q
    q2, qAq, b0 = (float(np.real(v)) for v in (q @ q, q @ Aq, Aq @ Aq))
    trA, trA2 = np.trace(A), np.trace(A @ A)
    X = XG * (1 - XG) * q2
    c0 = -XG * (1 - XG) * qAq
    b2 = (1 - 2 * XG)**2 * b0
    integrand = (6 * (trA * trA + 2 * trA2) / 24 * J(2, 0, X)
                 + 2 * (b2 / 4 + c0 * trA / 2) * J(1, 0, X)
                 - 0.5 * trh * trA * J(1, 1, X) + c0 * c0 * J(0, 0, X)
                 - c0 * trh * J(0, 1, X) + 0.25 * trh * trh * J(0, 2, X))
    I1 = 2 * C16 * float(J(2, 0, 0.0))
    return ((I1 / 4) * (0.5 * np.trace(h @ h) - 0.25 * trh * trh)
            - 0.25 * C16 * float(WG @ integrand))


def spin2(h, theta):
    return float(np.real(np.trace(theta @ h @ theta @ h) - np.trace(theta @ h)**2 / 3))


def spin0(h, theta):
    return float(np.real(np.trace(theta @ h)**2 / 3))


def form_factors(q, J):
    """F2, F0 of K(q) = Gamma2(q) - Gamma2(0), from a TT tensor and from theta."""
    theta = np.real(np.eye(4) - np.outer(q, q) / (q @ q))
    w, v = np.linalg.eigh(theta)
    e1, e2 = v[:, -1], v[:, -2]                          # two directions transverse to q
    tt = (np.outer(e1, e2) + np.outer(e2, e1)) / math.sqrt(2)
    zero = np.zeros(4)
    F2 = (gamma2(tt, q, J) - gamma2(tt, zero, J)) / spin2(tt, theta)
    F0 = (gamma2(theta, q, J) - gamma2(theta, zero, J)) / spin0(theta, theta)
    return F2, F0, theta


def audit(J, rho, stiffness, rng):
    """K1-K5 for one regulator weight."""
    zero = np.zeros(4)
    vac = 0.0
    for _ in range(5):
        h = rng.normal(size=(4, 4))
        h = h + h.T
        want = rho * (-0.25 * np.trace(h @ h) + 0.125 * np.trace(h)**2)
        vac = max(vac, abs(gamma2(h, zero, J) / want - 1))
    gap = (2 * C16 * float(J(2, 0, 0.0)) / 4) / (-rho)
    rows = []
    for size in (0.01, 0.1, 1.0):
        q = rng.normal(size=4)
        q *= size / np.linalg.norm(q)
        F2, F0, theta = form_factors(q, J)

        def K(a):
            return gamma2(a, q, J) - gamma2(a, zero, J)
        leak = resid = 0.0
        scale = abs(F2) + abs(F0)
        for _ in range(6):
            xi = rng.normal(size=4)
            hg = np.outer(q, xi) + np.outer(xi, q)
            h = rng.normal(size=(4, 4))
            h = h + h.T
            leak = max(leak, abs(K(hg)), abs(0.25 * (K(h + hg) - K(h - hg))))
            resid = max(resid, abs(K(h) - F2 * spin2(h, theta) - F0 * spin0(h, theta)))
        r = F2 / F0
        rows.append(dict(q=size, stiff_ratio=4 * F2 / size**2 / stiffness, f0_over_f2=F0 / F2,
                         gamma=(1 - r) / (2 + r), leak=leak / scale, residual=resid / scale))
    return dict(vacuum=vac, gap=gap, rows=rows)


def newton(J, stiffness, q=1e-3):
    """K5: exchange between two static sources, against -4 pi G M1 M2 / q^2."""
    qv = np.array([q, 0.0, 0.0, 0.0])
    F2, F0, theta = form_factors(qv, J)
    u = np.array([0.0, 0.0, 0.0, 1.0])                   # the time direction, u.q = 0
    T = np.outer(u, u)
    exchange = -(1.0 / 8.0) * (spin2(T, theta) / F2 + spin0(T, theta) / F0)
    G = 1.0 / (16 * math.pi * stiffness)
    return exchange / (-4 * math.pi * G / q**2)


def absorptive(s, t0, m2=1.0, eps=1e-9):
    """K7: F2 at q^2 = -s - i eps under the proper time cutoff (TT channel)."""
    x, w = np.polynomial.legendre.leggauss(2000)
    x, w = 0.5 * (x + 1), 0.5 * w

    def E3(z):
        e = exp1(z)
        for k in (1, 2):
            e = (np.exp(-z) - z * e) / k
        return e
    D = m2 - x * (1 - x) * s - 1j * eps
    return -0.125 * C16 * np.sum(w * (E3(t0 * D) - E3(t0 * m2 + 0j))) / t0**2


def trace_free_exchange(J, stiffness, q=1e-2):
    """K8: static exchange on the full and on the traceless space, over Newton."""
    qv, zero = np.array([q, 0.0, 0.0, 0.0]), np.zeros(4)
    basis = []
    for i in range(4):
        for j in range(i, 4):
            E = np.zeros((4, 4))
            E[i, j] = E[j, i] = 1.0 if i == j else 1 / math.sqrt(2)
            basis.append(E)

    def K(a):
        return gamma2(a, qv, J) - gamma2(a, zero, J)
    n = len(basis)
    M = np.array([[K(basis[a]) if a == b else 0.25 * (K(basis[a] + basis[b])
                                                      - K(basis[a] - basis[b]))
                   for b in range(n)] for a in range(n)])

    def vec(T):
        return np.array([np.trace(B @ T) for B in basis])
    u = np.array([0.0, 0.0, 0.0, 1.0])
    T, tr = vec(np.outer(u, u)), vec(np.eye(4))
    P = np.eye(n) - np.outer(tr, tr) / (tr @ tr)
    newton_law = -4 * math.pi / (16 * math.pi * stiffness) / q**2
    full = -0.125 * T @ np.linalg.pinv(M, rcond=1e-9) @ T / newton_law
    free = -0.125 * (P @ T) @ np.linalg.pinv(P @ M @ P, rcond=1e-9) @ (P @ T) / newton_law
    spin2_full = np.sort(np.linalg.eigvalsh(M))[-5:] / q**2
    spin2_free = np.sort(np.linalg.eigvalsh(P @ M @ P))[-5:] / q**2
    return full, free, float(np.max(abs(spin2_full / spin2_free - 1)))


def direction_content(J, h, q=1e-2, shift=0.0, n=40, seed=3):
    """K9: degree four share of K(h; q) over directions of q at fixed |q|."""
    rng = np.random.default_rng(seed)
    zero = np.zeros(4)
    dirs = rng.normal(size=(n, 4))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    vals, rows = [], []
    for u in dirs:
        qv = q * u
        k = gamma2(h, qv, J) - gamma2(h, zero, J)
        k += shift * q * q * (u @ h @ u)**2         # a term quartic in the direction
        vals.append(k)
        rows.append([u[i] * u[j] for i in range(4) for j in range(i, 4)])
    vals, rows = np.array(vals), np.array(rows)
    coef, *_ = np.linalg.lstsq(rows, vals, rcond=None)
    return float(np.linalg.norm(vals - rows @ coef) / np.linalg.norm(vals))


def timelike_ratio(J, s):
    q = np.array([1j * math.sqrt(s), 0.0, 0.0, 0.0])
    F2, F0, _ = form_factors(q, J)
    return F2 / (-s), F0 / (-s)


def run(report=print):
    rng = np.random.default_rng(0)
    m2 = 1.0
    weights = {"Pauli-Villars, L = 3 m": pauli_villars(m2, 3.0),
               "proper time, t0 = 1/(10 m)^2": proper_time(m2, 1e-2)}
    out = {}
    for name, (J, rho, stiff) in weights.items():
        a = audit(J, rho, stiff, rng)
        a["newton"] = newton(J, stiff)
        assert a["vacuum"] < 1e-9 and abs(a["gap"] - 1) < 1e-9
        for row in a["rows"]:
            assert row["leak"] < 1e-6 and row["residual"] < 1e-6      # round off at small q
        assert abs(a["rows"][0]["stiff_ratio"] - 1) < 1e-4
        assert abs(a["rows"][0]["f0_over_f2"] + 2) < 1e-3
        assert abs(a["rows"][0]["gamma"] - 1) < 1e-3
        assert abs(a["newton"] - 1) < 1e-4 and stiff > 0
        out[name] = a

    # The proper time weight: positive metric from a positive regulator.
    Jpt, rho_pt, stiff_pt = weights["proper time, t0 = 1/(10 m)^2"]
    I1_pt = 2 * C16 * float(Jpt(2, 0, 0.0))
    assert I1_pt > 0 > rho_pt
    t0 = 1e-2
    G_inverse = 16 * math.pi * N_COMPONENTS * stiff_pt
    fact1 = (N_COMPONENTS / 6) / (2 * math.pi * t0) * expn(2, t0 * m2)
    assert abs(G_inverse / fact1 - 1) < 1e-12
    assert abs(stiff_pt / cs.a_li([m2], 10.0, "pt") - 1) < 1e-12

    # K6: the Euclidean axis, and the region below threshold.
    Jpv = weights["Pauli-Villars, L = 3 m"][0]
    e = np.array([1.0, 0.0, 0.0, 0.0])
    signs = [(form_factors(x * e, J)[0] > 0, form_factors(x * e, J)[1] < 0)
             for J in (Jpv, Jpt) for x in np.geomspace(1e-2, 1e3, 60)]
    assert all(a and b for a, b in signs)
    departure, below = {}, {}
    for lam in (3.0, 10.0, 30.0, 100.0):
        J = pauli_villars(m2, lam)[0]
        F2, F0, _ = form_factors(e, J)
        departure[lam] = F0 / F2 + 2
        ratios = [timelike_ratio(J, s) for s in np.linspace(0.05, 3.999, 80)]
        below[lam] = (all(f2 > 0 for f2, _ in ratios), all(f0 < 0 for _, f0 in ratios))
    assert all(below[lam] == (True, True) for lam in (10.0, 30.0, 100.0))
    assert below[3.0] == (True, False)                   # regulator at the carrier scale
    scaling = math.log(departure[10.0] / departure[100.0]) / math.log(10.0)
    assert 1.5 < scaling < 2.1

    # K7: the absorptive part is the two particle cut, for every t0.
    cut = {}
    for s in (5.0, 20.0):
        want = s * s * (1 - 4 * m2 / s)**2.5 / (7680 * math.pi)
        cut[s] = [abs(absorptive(s, t).imag) / want for t in (1e-2, 1e-3)]
        assert max(abs(c - 1) for c in cut[s]) < 1e-4
    below_cut = abs(absorptive(3.0, 1e-2).imag)
    assert below_cut < 1e-8

    # K8: the traceless form.
    full, free, spin2_dev = trace_free_exchange(Jpt, stiff_pt)
    assert abs(full - 1) < 1e-4 and abs(free - 1) < 1e-4 and spin2_dev < 1e-8

    # K9: the q^2 term is local only at F0/F2 = -2.
    h9 = rng.normal(size=(4, 4))
    h9 = h9 + h9.T
    local = direction_content(Jpt, h9)
    off_ratio = direction_content(Jpt, h9, shift=0.1 * stiff_pt)
    assert local < 1e-4 and off_ratio > 1e-2
    h9_tl = h9 - np.trace(h9) / 4 * np.eye(4)                # the traceless space
    local_tl = direction_content(Jpt, h9_tl)
    off_ratio_tl = direction_content(Jpt, h9_tl, shift=0.1 * stiff_pt)
    assert local_tl < 1e-4 and off_ratio_tl > 1e-2

    for name, a in out.items():
        report(f"[{name}] K1 vacuum term structure to {a['vacuum']:.1e}; K2 gap equation "
               f"<d phi d phi>/V0 = {a['gap']:.12f} delta.")
        for row in a["rows"]:
            report(f"   |q| = {row['q']:5.2f} m: K3 gauge leak {row['leak']:.1e}; K4 residual "
                   f"{row['residual']:.1e}; K5 4F2/q^2 over the a1 stiffness = "
                   f"{row['stiff_ratio']:.6f}, F0/F2 = {row['f0_over_f2']:+.6f}, gamma = "
                   f"{row['gamma']:.6f}.")
        report(f"   K5 static exchange over -4 pi G M1 M2/q^2: {a['newton']:.6f}.")
    report(f"Proper time: I1 = {I1_pt:.4f} > 0, V0 = {-rho_pt:.4f} > 0; 1/G = "
           f"{G_inverse:.6f} = (N/6)/(2 pi t0) E2(m^2 t0) = {fact1:.6f} for N = {N_COMPONENTS}.")
    report("K6 F2/q^2 > 0 > F0/q^2 on the Euclidean axis from 0.01 m to 1000 m, both weights.")
    for lam in (3.0, 10.0, 30.0, 100.0):
        report(f"   L = {lam:5.0f} m: F0/F2 + 2 at |q| = m is {departure[lam]:+.2e}. Below "
               f"threshold F2/q^2 > 0: {below[lam][0]}, F0/q^2 < 0: {below[lam][1]}.")
    report(f"   departure falls as (m/L)^{scaling:.2f}.")
    report(f"K7 proper time, |Im F2| over s^2 beta^5/(7680 pi): s = 5: {cut[5.0][0]:.6f}, "
           f"{cut[5.0][1]:.6f} (t0 = 1e-2, 1e-3); s = 20: {cut[20.0][0]:.6f}, {cut[20.0][1]:.6f}; "
           f"below threshold {below_cut:.1e}.")
    report(f"K8 static exchange over Newton: full form {full:.6f}, traceless form {free:.6f}; "
           f"spin-2 eigenvalues agree to {spin2_dev:.1e}.")
    report(f"K9 degree four share of the kernel over directions of q: {local:.1e}. With a "
           f"term off the Einstein ratio added, {off_ratio:.2f}: the q^2 term is local only "
           f"at F0/F2 = -2, at every order. On traceless h {local_tl:.1e} and {off_ratio_tl:.2f}.")
    report("RESULT: for the carrier in the empty field, at one loop under a covariant "
           "regulator, the induced action of the composite metric is transverse at "
           "every momentum, has two form factors, and at low momentum is the "
           "linearized Einstein-Hilbert operator with stiffness a1 = 1/6 per real "
           "component. Its inverse between conserved sources is the massless pole: "
           "Newton's law and gamma = 1. The flat metric is the fluctuation term at "
           "the point where the vacuum term cancels. The proper time cutoff is such "
           "a regulator, and with it 1/G is FACT 1's Sigma/(2 pi t0).")
    return dict(audits=out, I1_proper_time=I1_pt, rho_proper_time=rho_pt,
                G_inverse=G_inverse, fact1=fact1, departure=departure, below=below,
                scaling=scaling, cut=cut, below_cut=below_cut,
                exchange_full=full, exchange_trace_free=free,
                local=local, off_ratio=off_ratio)


if __name__ == "__main__":
    run()
