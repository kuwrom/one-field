"""Sakharov induced tensor speed of the relativistic Z3 condensate.

PREMISES (each one is a premise, none is derived here):
  P1  the relativistic level: three complex Klein-Gordon fields phi_a,
      L = sum_a |d phi_a|^2 - m^2 |phi|^2 - V(phi), V the U(1) x Z3
      quartic whose nonrelativistic limit is the committed NLS
      (V(chi) = V_NLS(sqrt(2m) chi), g1/g0 = nls_soliton.G1_OVER_G0),
      condensed at mu_t = m + mu_NR (Fagnocchi et al. arXiv:1001.1044);
  P2  the metric variable is the actionless metric of the carrier
      written without a background (tests/probes/composite_metric.py).
      Its induced dynamics is computed: phi couples covariantly to
      g_mu nu with lapse 1, shift 0, g_ij = -(delta + h)_ij;
  P3  canonical quantization, one loop, zero temperature, hbar = c = 1.

SPECTRUM. With phi = exp(-i mu_t t)(Phi + eta), mu_t^2 = m^2 + 2 m mu_N
(mu_N = mu_NR (1 + mu_NR/2m) is the NLS chemical potential), the
fluctuation operator is block diagonal in the aligned Z3 frame:
    D(w, p) = [[p^2 + 2m a_u - w^2 - 2 mu_t w, 2m b],
               [2m b, p^2 + 2m a_v - w^2 + 2 mu_t w]],
(a_u, a_v, b) = (mu_N, mu_N, mu_N), (2r, r/2, r/2), (r/2, 2r, r/2), the
entries of gravity_tensor's A, B. det D = 0 gives six positive branches:
the phonon, w^2 = c_s^2 p^2, c_s^2 = mu_N/(m + 3 mu_N) (-> 1/3 at
mu_N >> m), the relative doublet, and three antiparticle like branches
(phonon partner at 2 sqrt(mu_t^2 + m mu_N)). All have dw/dp -> 1 at
p >> m, and the lower three reduce to clock_frequencies at mu_N << m.

TT SELF-ENERGY. The vertex is (1/2) h_ij T^ij = h_ij d_i phi* d_j phi (TT),
the same on both Nambu components, so with Euclidean G = D_E^-1 and
symmetric routing (e_ij p_i q_j = e_ij p_i p_j for TT e),
    Pi(K) = seagull - (1/2) int d^4p/(2pi)^4 (e.p.p)^2 tr G(p-K/2) G(p+K/2).
The seagull is K independent. Integrating by parts in p4 and p_z,
    A = (1/4) int (e.p.p)^2 tr (d_4 G)^2,   B = (1/4) int (e.p.p)^2 tr (d_z G)^2,
c_T^2 = B/A. For O(4)-invariant G the two are equal: c_T = 1.

REGULATOR. A = A_LI + dA with A_LI the Lorentz invariant reference
(B = 0: masses M_i^2 = m^2 + Hessian of V) and dA = cond - ref, which
is UV finite (integrand ~ p^-3, checked). Covariant regulators act on
A_LI: Pauli-Villars Tr log D -> sum_j c_j Tr log(D + s_j), c = (1,-2,1),
s = (0, L^2, 2 L^2), A_LI = sum_ij c_j (M_i^2+s_j) ln(M_i^2+s_j)/(192 pi^2)
(the PV copies carry the condensate, so the full PV answer is used).
Proper time, A_LI = sum_i L^2 [e^-y - y E1(y)]/(192 pi^2), y = M_i^2/L^2.
Both give A_LI = B_LI. A spatial cutoff |p| < L gives B/A = 5/7, and an
O(4) cutoff on the corotating loop momentum gives c_T^2 - 1 = 0.40 mu_t^2/L^2
in the free vacuum: neither is covariant for the phi quanta.

RESULTS (numbers printed by run()):
  dA -> A_NLS and dB -> B_NLS of induced_horava.py (physical units) as
  mu_N/m -> 0: the condensate's own tensor response is the Galilean one.
  dB - dA = -kappa_1 m^2 sqrt(mu_N/m) (1 + O(mu_N/m)), kappa_1 = 3.121e-4,
  so
      delta = c_T - 1 = 16 pi G (dB - dA) = -0.0157 m^(3/2) mu_N^(1/2) / M_Pl^2,
  with 1/(32 pi G) the total Lorentz invariant stiffness. Power in m/L:
  2.000. Power in mu_N/m: 0.500. Delta = -kappa (m/L)^2 sqrt(mu_N/m) with
  kappa = 0.0711 (PV), 0.0986 (proper time). c_T < c.
  delta table (PV on the full condensed loop; proper time in brackets,
  which regulates only A_LI and so is meaningful for m << L):
      m/L     mu_NR/m = 0.1           mu_NR/m = 0.01
      1       -5.92e-2 (-1.93e-1)     -1.85e-2 (-6.23e-2)
      0.1     -2.34e-4 (-3.30e-4)     -7.33e-5 (-1.03e-4)
      0.01    -2.24e-6 (-3.10e-6)     -7.07e-7 (-9.80e-7)
      1e-3    -2.23e-8 (-3.10e-8)     -7.06e-9 (-9.79e-9)
  No separation, 1/(xi L) = d10 = 2, L = 1/l_Pl (m/L = 4.36, 14.1):
  delta = -0.59, -0.87 (mu_NR/m = 0.1, 0.01). With the stiffness fixed
  at 1/(32 pi G_N) by the rest of the ledger instead, delta = -0.094,
  -0.31 at xi = l_Pl/d10.
  GW170817 (|delta| < 1e-15): m/L < 2.1e-7 (xi L > 1.0e7) at
  mu_NR/m = 0.1 and m/L < 3.8e-7 (xi L > 1.9e7) at 0.01. For any m,
  xi > 2.7e6 l_Pl (the minimum is the ultra relativistic condensate,
  dB - dA -> -1.43e-4 / xi^2). A healing length of l_Pl/2 is too short
  by a factor 5e6 in length, 13 orders of magnitude in delta. A uniform
  clock filling space is off the Lorentz invariant point. The vacuum of
  gate C is the empty field. The bound is computed on the uniform clock
  of P1, which has the gas sign. On the committed sign the medium exists
  only inside lumps, and closed_event_scale.py bounds the interior pixel
  of a lump by the carrier's Compton length.
  Horava lambda (k = 0 one loop kernel, fixed background): lambda = 1
  in the Lorentz invariant limit. Lambda - 1 = -(dK_tr + 2 dA)/(3 A_LI)
  = -6.7e-6, -2.1e-6 at m/L = 0.01 (mu_NR/m = 0.1, 0.01), -> 0 as
  (m/L)^2 sqrt(mu_N/m). The condensate kernel alone gives
  dK_tr/dA -> 5/2 at fixed background and 5/8 (lambda = 1/8, the value
  of induced_horava.py) with its corotating vertex. On the vacuum like
  modes that vertex adds a log divergent trace term ~ m mu_N ln L. For
  this reason lambda is quoted at fixed background.

Usage: python3 tests/probes/carrier_sakharov.py
"""

import math
import os
import sys

import numpy as np
from scipy.special import exp1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from gravity_tensor import clock_frequencies  # noqa: E402
from induced_horava import a_integrand, b_integrand, channels, integrate  # noqa: E402
from nls_soliton import G1_OVER_G0  # noqa: E402
from root import d10  # noqa: E402

R_OVER_MU = G1_OVER_G0 / (3 - G1_OVER_G0)          # r / mu_N on the clock
MU0 = 3 - G1_OVER_G0                               # mu_N at g0 rho = 1
PV = ((1.0, 0.0), (-2.0, 1.0), (1.0, 2.0))         # (c_j, s_j / L^2)
PRE = 0.25 / (2 * math.pi)**4 * (8 * math.pi / 15)  # (1/4) x TT angular measure
GW_BOUND = 1e-15


def gl_log(lo, hi, n=16, ppd=2):
    """Gauss-Legendre nodes, log uniform panels (ppd per decade)."""
    k = max(1, int(math.ceil(math.log10(hi / lo) * ppd)))
    x, w = np.polynomial.legendre.leggauss(n)
    e = np.log(np.logspace(math.log10(lo), math.log10(hi), k + 1))
    u = (e[:-1, None] + (e[1:] - e[:-1])[:, None] * (x + 1) / 2).ravel()
    return np.exp(u), (np.exp(u.reshape(k, n)) * ((e[1:] - e[:-1]) / 2)[:, None] * w).ravel()


def blocks(m, muN):
    """(a_u, a_v, b, a_u a_v - b^2) of the three Nambu blocks, and mu_t."""
    r = R_OVER_MU * muN
    return ([(muN, muN, muN, 0.0), (2 * r, r / 2, r / 2, 0.75 * r * r),
             (r / 2, 2 * r, r / 2, 0.75 * r * r)], math.sqrt(m * m + 2 * m * muN))


def hessian_masses(m, muN):
    mut = blocks(m, muN)[1]
    return [mut**2 + 2 * m * a for a in (muN, 2 * R_OVER_MU * muN, R_OVER_MU * muN / 2)]


def branches(p, m, muN):
    """Six positive corotating frequencies from det D(w, p) = 0."""
    bl, mut = blocks(m, muN)
    X, out = p * p, []
    for au, av, b, dd in bl:
        A1, A2 = X + mut**2 + 2 * m * au, X + mut**2 + 2 * m * av
        c0 = X * X + 2 * m * X * (au + av) + 4 * m * m * dd
        roots = np.roots([1.0, 0.0, -(2 * mut**2 + A1 + A2), -2 * mut * (A2 - A1), c0])
        assert max(abs(roots.imag)) < 1e-6 * max(1.0, max(abs(roots.real)))
        out += [w for w in roots.real if w > 0]
    return np.sort(np.array(out))


def propagator(p4, X, m, mut, blk, s):
    """Entries g11, g22, g12 of G = (D_E + s)^-1, cancellation free determinant."""
    au, av, b, dd = blk
    Y, Q = X + s + p4**2, 2 * mut * p4
    eu, ev = X + s + 2 * m * au, X + s + 2 * m * av
    det = Y * Y + 2 * m * Y * (au + av) + 4 * m * m * dd + 1j * Q * (eu - ev) + Q * Q
    return (ev + p4**2 + 1j * Q) / det, (eu + p4**2 - 1j * Q) / det, -2 * m * b / det


def tt_traces(p4, X, m, mut, blk, s):
    """Re tr (d_4 G)^2 and Re tr (d_X G)^2 for one block."""
    g11, g22, g12 = propagator(p4, X, m, mut, blk, s)
    a, b = 2 * (p4 - 1j * mut), 2 * (p4 + 1j * mut)
    q11, q22 = g11 * a * g11 + g12 * b * g12, g12 * a * g12 + g22 * b * g22
    q12 = g11 * a * g12 + g12 * b * g22
    s11, s22, s12 = g11**2 + g12**2, g12**2 + g22**2, g12 * (g11 + g22)
    return (q11**2 + q22**2 + 2 * q12**2).real, (s11**2 + s22**2 + 2 * s12**2).real


def grid(lo, hi, n=16, ppd=2):
    p, wp = gl_log(lo, hi, n, ppd)
    q, wq = gl_log(lo * 1e-2, hi * 10, n, ppd)
    return p, wp, np.r_[q, -q][:, None], np.r_[wq, wq][:, None]


def tt_condensate(m, muN, s=0.0, hi=1e4):
    """(dA, dB): condensate minus Lorentz invariant reference, UV finite."""
    bl, mut = blocks(m, muN)
    p, wp, P4, W4 = grid(1e-8 * m, hi * math.sqrt(m * m + s))
    X, dA, dB = (p**2)[None, :], 0.0, 0.0
    for blk in bl:
        c4, cX = tt_traces(P4, X, m, mut, blk, s)
        r4, rX = tt_traces(P4, X, m, mut, (blk[0], blk[1], 0.0, blk[0] * blk[1]), s)
        dA += PRE * np.sum(wp * p**6 * ((c4 - r4) * W4).sum(0))
        dB += PRE * (4 / 7) * np.sum(wp * p**8 * ((cX - rX) * W4).sum(0))
    return dA, dB


def tt_vacuum(M2, mut, lam, regulator):
    """Numerical (A, B) of the free reference: 'pv' or 'spatial' |p| < lam."""
    pv = regulator == "pv"
    p, wp, P4, W4 = grid(1e-6 * lam, (1e4 if pv else 1.0) * lam, 20, 3)
    X, fA, fB = (p**2)[None, :], 0.0, 0.0
    for Mi in M2:
        blk = (0.5 * (Mi - mut**2),) * 2 + (0.0, 0.25 * (Mi - mut**2)**2)
        for c, s in (PV if pv else ((1.0, 0.0),)):
            t4, tX = tt_traces(P4, X, 1.0, mut, blk, s * lam**2)
            fA, fB = fA + c * t4, fB + c * tX
    return (PRE * np.sum(wp * p**6 * (fA * W4).sum(0)),
            PRE * (4 / 7) * np.sum(wp * p**8 * (fB * W4).sum(0)))


def tt_o4_corotating(M2, mut, lam, n=400):
    """Free field, O(4) cutoff on the corotating loop momentum: B/A - 1."""
    xb, wb = np.polynomial.legendre.leggauss(n)
    beta = ((xb + 1) * math.pi / 2)[None, :]
    rho, wr = gl_log(1e-6 * lam, lam, 24, 4)
    P4, P = rho[:, None] * np.cos(beta), rho[:, None] * np.sin(beta)
    t4, tX = tt_traces(P4, P**2, 1.0, mut, (0.5 * (M2 - mut**2),) * 2 + (0.0, 0.25 * (M2 - mut**2)**2), 0.0)
    J = rho[:, None] * wr[:, None] * (wb * math.pi / 2)[None, :]
    return (4 / 7) * np.sum(J * P**8 * tX) / np.sum(J * P**6 * t4) - 1


def a_li(M2, lam, regulator):
    if regulator == "pv":
        return sum(c * (Mi + s * lam**2) * math.log(Mi + s * lam**2)
                   for Mi in M2 for c, s in PV) / (192 * math.pi**2)
    return sum(lam**2 * (math.exp(-Mi / lam**2) - Mi / lam**2 * exp1(Mi / lam**2))
               for Mi in M2) / (192 * math.pi**2)


def mm(A, B):
    return (A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3],
            A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3])


def kinetic_k0(m, muN, blk, s, corotate, scale=None):
    """Integrated k = 0 w^2 kernels (TT, trace e = delta/sqrt3) of one block.

    Vertex t (D_E(p4) - W^2/4 + bg) - (e.p.p), t = tr e / 2. Bg is the
    corotating response of the condensate (induced_horava convention).
    """
    mut = blocks(m, muN)[1]
    scale = scale or math.sqrt(m * m + s)        # PV sums need one common grid
    p, wp, P4, W4 = grid(1e-8 * scale, 1e4 * scale)
    X = (p**2)[None, :]
    g11, g22, g12 = propagator(P4, X, m, mut, blk, s)
    Q = 2 * mut * P4
    au, av, b = blk[:3]
    D = (X + s + 2 * m * au + P4**2 - 1j * Q, 2 * m * b + 0 * P4, 2 * m * b + 0 * P4,
         X + s + 2 * m * av + P4**2 + 1j * Q)
    if corotate:
        kth = 1 / (1 / mut + 2 * mut / (mut**2 - m * m))
        kn = 2 * mut * kth / (mut**2 - m * m)
        D = (D[0] - kth * (-2j * (P4 - 1j * mut)) - kn * 2 * m * (au + muN),
             D[1] * (1 - kn), D[2] * (1 - kn),
             D[3] - kth * (2j * (P4 + 1j * mut)) - kn * 2 * m * (av + muN))
    G, zero = (g11, g12, g12, g22), 0 * P4 + 0j
    Dp = (2 * (P4 - 1j * mut), zero, zero, 2 * (P4 + 1j * mut))
    I = (1 + zero, zero, zero, 1 + zero)
    GDG = mm(G, mm(Dp, G))
    G1 = tuple(-x for x in GDG)
    G2 = tuple(2 * x - 2 * y for x, y in zip(mm(GDG, mm(Dp, G)), mm(G, G)))
    tr = lambda A: (A[0] + A[3]).real
    C = lambda M, N: (tr(mm(G2, mm(M, mm(G, N)))) / 8 + tr(mm(G, mm(M, mm(G2, N)))) / 8
                      - tr(mm(G1, mm(M, mm(G1, N)))) / 4)
    cii, cdd, cdi = C(I, I), C(D, D), C(D, I) + C(I, D)
    t, S, S2 = math.sqrt(3) / 2, p**2 / math.sqrt(3), p**4 / 3
    ktt = -0.5 * (2 / 15) * p**4 * cii
    ktr = -0.5 * (t * t * cdd - t * S * cdi + S2 * cii
                  - t / 2 * (t * tr(mm(G, mm(G, D))) - S * tr(mm(G, G))))
    f = 4 * math.pi / (2 * math.pi)**4
    return np.array([f * np.sum(wp * p**2 * (k * W4).sum(0)) for k in (ktt, ktr)])


def kinetic_condensate(m, muN, corotate=False):
    """(dK_TT, dK_tr): condensate minus reference (same vertex convention)."""
    out = np.zeros(2)
    for blk in blocks(m, muN)[0]:
        ref = (blk[0], blk[1], 0.0, blk[0] * blk[1])
        out += kinetic_k0(m, muN, blk, 0.0, corotate) - kinetic_k0(m, muN, ref, 0.0, corotate)
    return out


def delta_full(m_over_lam, x_nr):
    """(delta_PV full, delta_PT, delta leading) at m = 1."""
    muN, lam = x_nr * (1 + x_nr / 2), 1 / m_over_lam
    M2 = hessian_masses(1.0, muN)
    parts = [(c, tt_condensate(1.0, muN, s * lam**2)) for c, s in PV]
    dA = sum(c * v[0] for c, v in parts)
    dB = sum(c * v[1] for c, v in parts)
    a_pv, a_pt = a_li(M2, lam, "pv"), a_li(M2, lam, "pt")
    dA0, dB0 = parts[0][1]
    return (math.sqrt((a_pv + dB) / (a_pv + dA)) - 1,
            math.sqrt((a_pt + dB0) / (a_pt + dA0)) - 1, (dB0 - dA0) / (2 * a_pv))


def run(report=print):
    # 1. Spectrum.
    m, x = 1.0, 0.1
    mut = blocks(m, x)[1]
    w0 = branches(0.0 + 1e-9, m, x)
    cs2 = (branches(1e-4, m, x)[0] / 1e-4)**2
    assert len(w0) == 6 and abs(cs2 / (x / (m + 3 * x)) - 1) < 1e-6
    assert min(abs(w0 - 2 * math.sqrt(mut**2 + m * x))) < 1e-8
    vg = (branches(1e3 + 1e-3, m, x) - branches(1e3 - 1e-3, m, x)) / 2e-3
    assert max(abs(vg - 1)) < 1e-5
    xs = 1e-4
    nr_err = max(max(abs(branches(p, m, xs)[:3] / np.sort(clock_frequencies(
        p, xs, R_OVER_MU * xs)) - 1)) for p in (0.003, 0.01, 0.03))
    assert nr_err < 1e-3
    report(f"1. Spectrum (mu_N/m = {x}): six branches, gaps at p = 0: {np.round(w0, 6)}; "
           f"c_s^2 = mu_N/(m + 3mu_N) = {cs2:.6f}. Dw/dp at p = 1e3 m: |v_g - 1| <= "
           f"{max(abs(vg - 1)):.1e}. Lower three = clock_frequencies at mu_N/m = 1e-4 "
           f"(rel. err. {nr_err:.1e})")

    # 2. Lorentz invariant reference: covariant regulator gives c_T = c and lambda = 1.
    M2, lam = hessian_masses(m, x), 3.0
    A_pv, B_pv = tt_vacuum(M2, mut, lam, "pv")
    li_err = max(abs(B_pv / A_pv - 1), abs(A_pv / a_li(M2, lam, "pv") - 1))
    assert li_err < 2e-6
    kin = sum(sum(c * kinetic_k0(m, x, (0.5 * (Mi - mut**2),) * 2 + (0.0, 0.25 * (Mi - mut**2)**2),
                                 s * lam**2, False, lam) for c, s in PV) for Mi in M2)
    lam_li = (1 - kin[1] / kin[0]) / 3
    assert abs(kin[0] / A_pv - 1) < 2e-6 and abs(lam_li - 1) < 2e-6
    report(f"2. Free reference, Pauli-Villars, corotating frame (mu_t = {mut:.4f}): "
           f"A = {A_pv:.9e}, B/A - 1 = {B_pv / A_pv - 1:.1e}, analytic A_LI agrees "
           f"({li_err:.1e}); k = 0 kernel: trace/TT = {kin[1] / kin[0]:.8f}, "
           f"lambda = {lam_li:.8f}")
    A_sp, B_sp = tt_vacuum(M2, mut, 1e3, "spatial")
    o4 = tt_o4_corotating(1.0, 0.9, 100.0)
    assert abs(B_sp / A_sp - 5 / 7) < 1e-4 and o4 > 0
    report(f"   spatial cutoff |p| < L: B/A = {B_sp / A_sp:.6f} (5/7), delta = "
           f"{math.sqrt(B_sp / A_sp) - 1:+.4f} at every separation; O(4) cutoff on the "
           f"corotating momentum, free field M = 1, mu = 0.9, L = 100: "
           f"(B/A - 1) L^2/mu^2 = {o4 * 1e4 / 0.81:.4f}")

    # 3. Condensate part: UV finite, Galilean limit = induced_horava.
    nls = {n: (integrate(a_integrand, v), integrate(b_integrand, v))
           for n, v in channels(MU0, G1_OVER_G0).items()}
    A_nls = sum(v[0] for v in nls.values()) / math.sqrt(MU0)       # x m^2 sqrt(mu_N/m)
    B_nls = sum(v[1] for v in nls.values()) / MU0**1.5            # x m^2 (mu_N/m)^1.5
    scan = {x: tt_condensate(1.0, x) for x in (1.0, 0.1, 0.01, 1e-3, 1e-4)}
    for x, (dA, dB) in scan.items():
        report(f"3. mu_N/m = {x:6.0e}: dA = {dA:.6e}, dB = {dB:.6e}, dA/A_NLS = "
               f"{dA / (A_nls * x**0.5):.5f}, dB/B_NLS = {dB / (B_nls * x**1.5):.5f}, "
               f"(dB - dA)/sqrt(x) = {(dB - dA) / x**0.5:.6e}")
    assert abs(scan[1e-4][0] / (A_nls * 1e-2) - 1) < 1e-3
    assert abs(scan[1e-4][1] / (B_nls * 1e-6) - 1) < 1e-2
    f_pow = math.log((scan[1e-3][1] - scan[1e-3][0]) / (scan[1e-4][1] - scan[1e-4][0])) / math.log(10)
    kappa1 = -(scan[1e-4][1] - scan[1e-4][0]) / 1e-2
    tail = [tt_condensate(1.0, 0.1, hi=h)[0] for h in (1e1, 1e2)]
    uv_ratio = (tail[1] - tail[0]) / (scan[0.1][0] - tail[1])
    assert abs(f_pow - 0.5) < 0.01 and 50 < uv_ratio < 200
    report(f"   A_NLS, B_NLS live from induced_horava. Power in mu_N/m = {f_pow:.4f}; "
           f"kappa_1 = {kappa1:.5e}. Truncation at 1e1/1e2/inf m: tail ratio "
           f"{uv_ratio:.1f} (integrand ~ p^-3, finite)")

    # 4. delta table.
    table = {}
    for xnr in (0.1, 0.01):
        for ml in (1.0, 0.1, 0.01, 1e-3):
            table[(xnr, ml)] = delta_full(ml, xnr)
            d = table[(xnr, ml)]
            report(f"4. mu_NR/m = {xnr:5}, m/L = {ml:6.0e}: delta_PV = {d[0]:+.4e}, "
                   f"delta_PT = {d[1]:+.4e}, leading (dB dA)/(2A_LI) = {d[2]:+.4e}")
    pw = {xnr: math.log(table[(xnr, 0.01)][0] / table[(xnr, 1e-3)][0]) / math.log(10)
          for xnr in (0.1, 0.01)}
    assert all(abs(v - 2) < 1e-3 for v in pw.values())
    k_pv, k_pt = 16 * math.pi**2 * kappa1 / math.log(2), 32 * math.pi**2 * kappa1
    k_G = 16 * math.pi * kappa1
    report(f"   power in m/L: {pw}. Delta = -kappa (m/L)^2 sqrt(mu_N/m), kappa = "
           f"{k_pv:.4f} (PV), {k_pt:.4f} (PT). Regulator free: delta = 16 pi G (dB - dA) "
           f"= -{k_G:.4f} m^1.5 mu_N^0.5 / M_Pl^2")

    # 5. lambda (fixed background, the 1PI one loop kernel) and the NR cross check.
    kc = {x: kinetic_condensate(1.0, x * (1 + x / 2)) for x in (0.1, 0.01)}
    for xnr, ml in ((0.1, 0.01), (0.01, 0.01), (0.1, 1e-3), (0.01, 1e-3)):
        k = kc[xnr]
        a = a_li(hessian_masses(1.0, xnr * (1 + xnr / 2)), 1 / ml, "pv")
        lam_h = (1 - (-2 * a + k[1]) / (a + k[0])) / 3
        report(f"5. mu_NR/m = {xnr}, m/L = {ml:.0e}: dK_tr/dA = {k[1] / k[0]:.4f}, "
               f"lambda - 1 = {lam_h - 1:+.4e}")
    k_fix, k_rot = kinetic_condensate(1.0, 1e-4), kinetic_condensate(1.0, 1e-4, True)
    lam_rot = (1 - k_rot[1] / k_rot[0]) / 3
    assert abs(k_fix[1] / k_fix[0] - 2.5) < 0.01 and abs(lam_rot - 0.125) < 0.005
    report(f"   condensate kernel alone at mu_N/m = 1e-4: trace/TT = {k_fix[1] / k_fix[0]:.4f} "
           f"(fixed background), {k_rot[1] / k_rot[0]:.4f} corotating -> lambda = "
           f"{lam_rot:.4f} (induced_horava: 1/8)")

    # 6. No separation and the GW170817 bound.
    nosep = {}
    for xnr in (0.1, 0.01):
        muN = xnr * (1 + xnr / 2)
        ml = d10 / math.sqrt(2 * muN)             # 1/(xi L) = d10, 1/xi = sqrt(2 m mu_N)
        dA, dB = scan.get(muN) or tt_condensate(1.0, muN)
        nosep[xnr] = (delta_full(ml, xnr)[0], 16 * math.pi * (dB - dA) * ml**2)
        report(f"6. no separation, mu_NR/m = {xnr}: m/L = {ml:.3f}, delta_PV = "
               f"{nosep[xnr][0]:+.4f}. With G = G_N at xi = l_Pl/d10: delta = {nosep[xnr][1]:+.4f}")
    need = {}
    for xnr in (0.1, 0.01):
        d3 = abs(table[(xnr, 1e-3)][0])
        ml = 1e-3 * math.sqrt(GW_BOUND / d3)
        need[xnr] = (ml, 1 / (ml * math.sqrt(2 * xnr * (1 + xnr / 2))))
        report(f"   GW170817 at mu_NR/m = {xnr}: m/L < {ml:.2e}, xi L > {need[xnr][1]:.2e}")
    dA, dB = tt_condensate(1.0, 1e3, hi=1e6)
    g_ur = -(dB - dA) / 1e3 / 2                  # (dA - dB) xi^2 at mu_N >> m
    xi_min = math.sqrt(16 * math.pi * g_ur / GW_BOUND)
    report(f"   any m: dB - dA -> -{g_ur:.4e}/xi^2 (ultra relativistic minimum), so "
           f"xi > {xi_min:.2e} l_Pl. A healing length l_Pl/{d10} is short by a factor {xi_min * d10:.1e}")
    report("RESULT: with a covariant regulator and no condensate, c_T = c exactly. The "
           "condensate's Galilean inertia makes c_T < c by 16 pi G (dA - dB) ~ "
           "(m/L)^2 sqrt(mu_N/m). At a Planckian healing length delta = O(1). The "
           "vacuum of gate C is the empty field.")
    return dict(kappa1=kappa1, kappa_pv=k_pv, kappa_pt=k_pt, kappa_G=k_G, mu_power=f_pow,
                m_power=pw, table=table, nosep=nosep, gw_need=need, xi_min=xi_min,
                lambda_li=lam_li, lambda_corotating_nr=lam_rot, spatial_ratio=B_sp / A_sp,
                o4_artifact=o4, scan=scan, A_nls=A_nls, B_nls=B_nls)


if __name__ == "__main__":
    run()
