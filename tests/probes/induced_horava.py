"""One loop induced spatial metric action of the Z3-NLS clock vacuum.

HYPOTHESIS UNDER TEST. A uniform clock medium carries no tensor at first
order in the field (interference/gravity.py, gate T). The candidate takes the spatial
metric h_ij of the gauged Galilean symmetry (Newton-Cartan background,
Son arXiv:1306.0638 eq. (1); Jensen arXiv:1408.6855) as the variable
and lets its action be induced at one loop by the substrate's own
Bogoliubov quasiparticles. The result is a Horava type action
(Horava arXiv:0901.3775 eq. (2.11): K_ij K^ij - lambda K^2 plus a
potential). This probe computes its quadratic coefficients.

PREMISES (each one is a premise, none is derived here):
  P1  the gauging itself: psi couples to h_ij(t, x) as in Son eq. (1),
      measure sqrt(h), kinetic (1/2) h^ij d_i psi* d_j psi, contact
      interactions times sqrt(h). Lapse 1, shift 0, A_mu = 0;
  P2  canonical bosonic quantization of the quadratic fluctuations;
  P3  one loop (Gaussian fluctuations) about the clock vacuum
      phi_a = sqrt(rho) w^a, units hbar = m = 1, g0 rho = 1, g1 = g0/sqrt(2)
      (nls_soliton.G1_OVER_G0, the standing identification).

COUPLING. With the half density chi = h^(1/4) psi the time derivative
term (i/2)(chi* dchi/dt - c.c.) carries no metric, [chi, chi*] is
metric free, the number density is |chi|^2 and
    L = (i/2)(chi* chi' - c.c.) + mu |chi|^2
        - (1/2) sqrt(h) h^ij d_i(h^-1/4 chi*) d_j(h^-1/4 chi) - h^-1/2 V(chi).
Expanding h_ij = delta_ij + eps_ij gives S1 = (1/2) int eps_ij O_ij with
    O_ij = d_i eta* . d_j eta + delta_ij U(eta)       (pair part)
    U = eta* A eta + (1/2)(eta* B eta* + c.c.)       (A, B of gravity_tensor)
plus, for the trace only, the tree vertex (1/2) tr(eps) (mu + k^2/4) drho.
U uses the Hessian of V - mu N: a uniform h(t) rotates the condensate
phase (N is conserved), and the corotating frame removes mu |eta|^2.
The TT seagull is -(1/2)(eps^2)_ij <d_i eta* d_j eta> - (1/4) tr eps^2 <V>:
momentum and frequency independent, so it enters only Pi(0, 0).

INDUCED ACTION. Gamma = (1/2) int eps_ij Pi_ijkl eps_kl with
    Pi = seagull + (1/4) chi_OO,
    chi_OO(w, k) = int d^3p/(2pi)^3 (1/2) sum_nm |F_nm|^2 2E/(E^2 - w^2),
    F_nm = U_n(q)* M V_m(p)* + U_m(p)* M V_n(q)* + U_n* N U_m* + V_n* N* V_m*,
M = e_ij p_i q_j + tr(e) A, N = tr(e) B, E = E_n(q) + E_m(p). Only the
neutral pairs (0,0) and (+,-) survive; (+,+), (-,-), (0,+/-) vanish.
At k = 0 every channel gives the same universal vertex
    F = -(b/E_p) p_i p_j (e_ij - (1/2) delta_ij tr e),
so the w^2 kernel is proportional to tr(eps')^2 - (1/8)(tr eps')^2:
lambda = 1/8 exactly, for every channel, coupling and cutoff. In d
dimensions the same algebra gives lambda = (1 + d/2 - d^2/4)/2, which
is 1/2 = 1/d (the anisotropic Weyl value) in d = 2 where the contact
coupling is scale invariant.

TT SECTOR. Pi_TT = M2 + A w^2 - B k^2 + O(w^4, k^4 log k). A is a
positive spectral sum (healthy kinetic term); c_T^2 = B / A.

RESULT (clock vacuum, numbers printed by run()):
  A = 4.7266e-4 > 0,  B = 5.2026e-3 > 0,
  c_T = 3.3177 = 2.191 c_ph = 3.529 c_rel,  lambda = 1/8.
  Phonon channel alone: B/A = 5 mu to 1e-9 (5 c_ph^2).
  UV: A's integrand falls as p^-4, B's as p^-2 (B converges like 1/Lambda,
  so the free particle band p >> 1/xi supplies most of B). The
  Bogoliubov p^4/4 term is the regulator and no cutoff is needed.
  A hard cutoff on the relative momentum instead gives c_T/c_ph = 1.236,
  1.495, 1.855, 2.105 at 1, 2, 5, 20 / xi. The p << 1/xi shell ratio is
  8/7 mu, so no cutoff places c_T on the phonon cone.
  M2 (the TT pressure term) is linearly divergent and is not computed.
  Scalar sector: the trace couples linearly to the condensate density,
  so the phonon is a tree level scalar pole. With lambda = 1/8 outside
  [1/3, 1] the Horava scalar graviton is also present, nonghost.

VERDICT. The induced tensor cone is none of the medium's cones:
c_T / c_ph - 1 = 1.19 and c_T / c_rel - 1 = 2.53, against the GW170817
bound 1e-15 on |c_T/c_light - 1|, whichever cone would carry light.
lambda = 1/8, not the GR value 1. The tensor of gate T is the gradient
bilinear of the carrier in the empty field, induced at c_T = c with
lambda = 1 (tests/probes/composite_metric.py, carrier_sakharov.py).

Usage: python3 tests/probes/induced_horava.py
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from gravity_tensor import (bdg_matrix, clock_background,  # noqa: E402
                            tensor_spectral_density)
from nls_soliton import G1_OVER_G0  # noqa: E402

TT_ANGLE = 2.0 / 15.0          # sphere average of (e_ij n_i n_j)^2, e TT unit
TT_COS2 = 1.0 / 7.0            # <(e n n)^2 cos^2> / <(e n n)^2>, k along e's normal


def channels(mu, r):
    """(a, b, pair symmetry factor) of the two neutral pair channels."""
    return dict(phonon=(mu, mu, 0.5), relative=(5 * r / 4, r / 2, 1.0))


def coherence(X, a, b):
    """E, u, v of a Bogoliubov channel at p^2 = X (relative: splitting removed)."""
    eps = X / 2
    E = math.sqrt(eps * (eps + 2 * a) + (a - b) * (a + b))
    s = eps + a + E
    return E, math.sqrt(s / (2 * E)), -b / math.sqrt(2 * E * s)


def pair_kernel(X, Y, a, b):
    """(u_P v_Q + v_P u_Q)^2 2/(E_P + E_Q) at P^2 = X, Q^2 = Y, w = 0."""
    Ex, ux, vx = coherence(X, a, b)
    Ey, uy, vy = coherence(Y, a, b)
    return (ux * vy + vx * uy)**2 * 2 / (Ex + Ey)


def a_integrand(p, a, b, s):
    """d A / dp per TT polarization: (1/4) s (e p p)^2 (b/E)^2 2/(2E)^3."""
    E = coherence(p * p, a, b)[0]
    return 0.25 * s / (2 * math.pi**2) * p * p * TT_ANGLE * p**4 * (b / E)**2 / (4 * E**3)


def b_integrand(p, a, b, s):
    """d B / dp per TT polarization from the exact k^2 Taylor coefficient.

    Symmetric routing P, Q = p +/- k/2 keeps e_ij q_i p_j = e_ij p_i p_j
    for TT e. With x = p^2: P^2 = x + pk cos + k^2/4, Q^2 = x - pk cos + k^2/4.
    """
    x = p * p
    h = 1e-3 * x
    along = lambda t: pair_kernel(x + t, x + t, a, b)
    across = lambda t: pair_kernel(x + t, x - t, a, b)
    d1 = (-along(2*h) + 8*along(h) - 8*along(-h) + along(-2*h)) / (12 * h)
    d2 = (-across(2*h) + 16*across(h) - 30*across(0) + 16*across(-h)
          - across(-2*h)) / (12 * h * h)
    dchi = s / (2 * math.pi**2) * x * TT_ANGLE * x * x * (d1 / 4 + TT_COS2 * x * d2 / 2)
    return -0.25 * dchi


def integrate(f, args, upper=np.inf):
    edges = [1e-9, 0.3, 1.0, 3.0, 10.0, 30.0, np.inf]
    pieces = [(lo, min(hi, upper)) for lo, hi in zip(edges, edges[1:]) if lo < upper]
    return sum(quad(f, lo, hi, args=args, limit=400, epsabs=0, epsrel=1e-10)[0]
               for lo, hi in pieces)


def chi_direct(k, a, b, s, cutoff=60.0):
    """TT static susceptibility at momentum k by 2D quadrature (check only)."""
    def radial(p):
        def ang(c):
            P2 = p * p + p * k * c + k * k / 4
            Q2 = p * p - p * k * c + k * k / 4
            return p**6 * (1 - c * c)**2 / 4 * pair_kernel(P2, Q2, a, b)
        return quad(ang, -1, 1, epsabs=0, epsrel=1e-11)[0]
    edges = (1e-9, 1.0, 5.0, 20.0, cutoff)
    return s / (4 * math.pi**2) * sum(quad(radial, lo, hi, epsabs=0, epsrel=1e-10,
                                           limit=200)[0] for lo, hi in zip(edges, edges[1:]))


def nambu_modes(A, B, p):
    """Positive norm eigenvectors (U, V), |U|^2 - |V|^2 = 1, of the 6x6 operator."""
    w, X = np.linalg.eig(bdg_matrix(A, B, p * np.array([0.0, 0.0, 1.0])))
    out = []
    for i in range(6):
        x = X[:, i]
        norm = np.vdot(x[:3], x[:3]).real - np.vdot(x[3:], x[3:]).real
        if norm > 0:
            out.append((w[i].real, x[:3] / math.sqrt(norm), x[3:] / math.sqrt(norm)))
    assert len(out) == 3
    return sorted(out, key=lambda t: t[0])


def pair_amplitudes(A, B, p, e_pp, tr_e):
    """k = 0 pair amplitudes F_nm of O = e_ij O_ij from raw eigenvectors."""
    modes = nambu_modes(A, B, p)
    M, N = e_pp * np.eye(3) + tr_e * A, tr_e * B
    rows = []
    for En, Un, Vn in modes:
        for Em, Um, Vm in modes:
            F = (Un.conj() @ M @ Vm.conj() + Um.conj() @ M @ Vn.conj()
                 + Un.conj() @ N @ Um.conj() + Vn.conj() @ N.conj() @ Vm.conj())
            rows.append((En + Em, abs(F)**2))
    return rows, modes


def second_order_check(A, B, p, n, e, step=1e-4):
    """Exact E0''(a) of the deformed quadratic Hamiltonian vs <H''> - pair sum."""
    n = np.asarray(n, float) / np.linalg.norm(n)
    pv = p * n

    def kernel(x):
        hinv = np.linalg.inv(np.eye(3) + x * e)
        scale = np.linalg.det(np.eye(3) + x * e)**-0.5
        return pv @ hinv @ pv / 2 * np.eye(3) + scale * A, scale * B

    def e0(x):
        K, Bx = kernel(x)
        L = np.block([[K, Bx], [-Bx.conj(), -K.conj()]])
        w = np.sort(np.linalg.eigvals(L).real)[3:]
        return 0.5 * (w.sum() - np.trace(K).real)

    numeric = (e0(step) - 2 * e0(0) + e0(-step)) / step**2
    s2 = np.trace(e)**2 / 4 + np.trace(e @ e) / 2
    K2, B2 = (pv @ e @ e @ pv) * np.eye(3) + s2 * A, s2 * B
    rows, modes = pair_amplitudes(A, B, p, pv @ e @ pv, np.trace(e))
    seagull = sum((Vn @ K2 @ Vn.conj()).real + (Vn @ B2 @ Un.conj()).real
                  for _, Un, Vn in modes)
    return numeric, seagull - 0.25 * sum(w / E for E, w in rows)


def kinetic_coefficient(A, B, e, nodes=6, azimuths=8):
    """w^2 coefficient (1/4) int (1/2) sum |F|^2 2/E^3 for a uniform deformation e.

    Product rule (Gauss-Legendre in cos x uniform azimuth) is exact for the
    quartic angular dependence of |F|^2.
    """
    c, wl = np.polynomial.legendre.leggauss(nodes)
    phis = 2 * np.pi * np.arange(azimuths) / azimuths
    dirs, wts = [], []
    for x, w in zip(c, wl):
        for f in phis:
            st = math.sqrt(1 - x * x)
            dirs.append(np.array([st * math.cos(f), st * math.sin(f), x]))
            wts.append(w / azimuths)

    def radial(p):
        total = 0.0
        for wt, d in zip(wts, dirs):
            rows, _ = pair_amplitudes(A, B, p, p * p * (d @ e @ d), np.trace(e))
            total += wt / 2 * sum(0.5 * w * 2 / E**3 for E, w in rows)
        return 0.25 * p * p / (2 * math.pi**2) * total
    return sum(quad(radial, lo, hi, limit=200, epsabs=0, epsrel=1e-8)[0]
               for lo, hi in ((1e-6, 1.0), (1.0, 4.0), (4.0, 15.0), (15.0, np.inf)))


def run(report=print):
    phi, mu, r, A, B = clock_background()
    c_ph, c_rel = math.sqrt(mu), math.sqrt(5 * r / 4)
    xi_inv = math.sqrt(2 * mu)                     # 1/xi = sqrt(2 m mu)/hbar
    ch = channels(mu, r)
    report(f"Clock vacuum: mu = {mu:.9f}, r = {r:.9f}, c_ph = {c_ph:.9f}, "
           f"c_rel = {c_rel:.9f}, 1/xi = {xi_inv:.6f} (g1/g0 = {G1_OVER_G0:.6f})")

    # 1. Raw eigenvectors reproduce the channel coherence factors. Only
    #    neutral pairs (0,0) and (+,-) carry TT weight.
    worst = 0.0
    for p in (0.3, 1.0, 2.5):
        rows, _ = pair_amplitudes(A, B, p, 1.0, 0.0)
        live = sorted((round(E, 9), w) for E, w in rows if w > 1e-12)
        assert len(live) == 3 and len(rows) == 9
        for name, (a, b, _) in ch.items():
            E, u, v = coherence(p * p, a, b)
            hits = [w for Et, w in live if abs(Et - 2 * E) < 1e-7]
            assert len(hits) == (1 if name == "phonon" else 2)
            worst = max(worst, max(abs(w - (2 * u * v)**2) for w in hits))
    assert worst < 1e-9
    report(f"1. 6x6 eigenvector pair amplitudes = (2uv)^2 per channel "
           f"(error {worst:.1e}); (+,+), (-,-), (0,+/-) vanish.")

    # 2. Second order perturbation theory equals the exact curvature of the
    #    deformed zero point energy, TT and trace vertices alike.
    e_tt = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]]) / math.sqrt(2)
    e_tr = np.eye(3) / math.sqrt(3)
    e_mix = np.array([[1.0, 0.3, 0], [0.3, 0, -0.2], [0, -0.2, 0.5]])
    hf = 0.0
    for p, n, e in ((0.4, (0, 0, 1), e_tr), (1.3, (1, 2, 0.5), e_tt),
                    (2.2, (0.3, -1, 1), e_mix), (0.9, (1, 0, 0), e_mix)):
        num, pert = second_order_check(A, B, p, n, e)
        hf = max(hf, abs(num - pert) / max(abs(num), 1e-12))
    assert hf < 1e-5
    report(f"2. Vertex check: exact E0'' vs <H''> - pair sum, max rel. error {hf:.1e}.")

    # 3. TT coefficients A (w^2) and B (k^2), channel by channel.
    coef = {}
    for name, (a, b, s) in ch.items():
        coef[name] = (integrate(a_integrand, (a, b, s)), integrate(b_integrand, (a, b, s)))
    A_tt = sum(v[0] for v in coef.values())
    B_tt = sum(v[1] for v in coef.values())
    assert A_tt > 0 and B_tt > 0
    cT = math.sqrt(B_tt / A_tt)
    ratio_ph = coef["phonon"][1] / coef["phonon"][0] / mu
    assert abs(ratio_ph - 5) < 1e-7
    for name, (Ac, Bc) in coef.items():
        report(f"3. {name:8s}: A = {Ac:.6e}, B = {Bc:.6e}, B/A = {Bc / Ac:.9f}")
    report(f"   total  : A = {A_tt:.6e} > 0, B = {B_tt:.6e} > 0. Phonon channel "
           f"B/A = {ratio_ph:.9f} mu")

    # 3a. A against the existing k = 0 TT spectral density (two polarizations).
    A_spec = 0.0
    for chan in ("phase", "relative"):
        lo = 1e-6 if chan == "phase" else 2 * coherence(0, *ch["relative"][:2])[0]
        f = lambda w: 0.25 * 0.5 * tensor_spectral_density(0, w, mu, r, chan) * 2 / w**3
        A_spec += sum(quad(f, x0, x1, limit=200, epsabs=0, epsrel=1e-9)[0]
                      for x0, x1 in ((lo, lo + 2), (lo + 2, lo + 20), (lo + 20, np.inf)))
    assert abs(A_spec / A_tt - 1) < 1e-6
    report(f"   A from gravity_tensor spectral density: {A_spec:.6e} (agrees)")

    # 3b. B against a direct 2D quadrature of chi(0, k) at a common cutoff.
    B60 = sum(integrate(b_integrand, v, upper=60.0) for v in ch.values())
    k1, k2 = 0.1, 0.05
    direct = []
    for k in (k1, k2):
        diff = sum(chi_direct(k, *v) - chi_direct(0, *v) for v in ch.values())
        direct.append(-0.25 * diff / k**2)
    B_rich = (4 * direct[1] - direct[0]) / 3
    assert abs(B_rich / B60 - 1) < 2e-4
    report(f"   B at |p| < 60: Taylor {B60:.7e}, direct 2D (Richardson) {B_rich:.7e}")

    # 4. UV structure and the healing scale share.
    tails = {}
    for name, (a, b, s) in ch.items():
        pa = [a_integrand(p, a, b, s) for p in (100.0, 300.0)]
        pb = [b_integrand(p, a, b, s) for p in (100.0, 300.0)]
        tails[name] = (math.log(pa[1] / pa[0]) / math.log(3),
                       math.log(pb[1] / pb[0]) / math.log(3))
        assert abs(tails[name][0] + 4) < 0.02 and abs(tails[name][1] + 2) < 0.02
    report(f"4. Large-p integrand powers (A, B): {tails}. Both converge, "
           "B only as 1/Lambda.")
    cut = {}
    for n_xi in (1, 2, 5, 20):
        lam = n_xi * xi_inv
        Ac = sum(integrate(a_integrand, v, upper=lam) for v in ch.values())
        Bc = sum(integrate(b_integrand, v, upper=lam) for v in ch.values())
        cut[n_xi] = (Ac / A_tt, Bc / B_tt, math.sqrt(abs(Bc / Ac)) / c_ph)
        report(f"   hard cutoff {n_xi:2d}/xi: A share {cut[n_xi][0]:.4f}, "
               f"B share {cut[n_xi][1]:+.4f}, c_T/c_ph {cut[n_xi][2]:.4f}")

    # 5. Scalar sector: lambda from rotational structure at k = 0.
    kin = {name: kinetic_coefficient(A, B, e)
           for name, e in (("tt", e_tt), ("trace", e_tr),
                           ("diag100", np.diag([1.0, 0, 0])))}
    assert abs(kin["tt"] / A_tt - 1) < 1e-6
    lam_h = (1 - kin["trace"] / kin["tt"]) / 3
    beta = kin["tt"] / 2
    alpha = (kin["trace"] - kin["tt"]) / 3
    assert abs(lam_h - 0.125) < 1e-7
    assert abs(kin["diag100"] - (alpha + 2 * beta)) < 1e-7 * kin["tt"]
    scalar_kin = (1 - lam_h) / (1 - 3 * lam_h)
    report(f"5. k = 0 w^2 kernel: TT {kin['tt']:.6e}, trace {kin['trace']:.6e} "
           f"= {kin['trace'] / kin['tt']:.9f} TT. Lambda = {lam_h:.9f} (1/8)")
    report(f"   isotropic form alpha(tr)^2 + 2 beta tr(.^2) holds on diag(1,0,0); "
           f"Horava scalar kinetic factor (1-l)/(1-3l) = {scalar_kin:.6f} > 0")

    # 5a. Tree level scalar pole: drho couples only to the phonon.
    res_err = 0.0
    for p in (0.2, 1.0, 3.0):
        for i, (E, U, V) in enumerate(nambu_modes(A, B, p)):
            amp = abs(phi.conj() @ V.conj() + phi @ U.conj())**2
            E0, _, _ = coherence(p * p, mu, mu)
            want = 3 * (p * p / 2) / E0 if abs(E - E0) < 1e-9 else 0.0
            res_err = max(res_err, abs(amp - want))
    assert res_err < 1e-9
    report("   trace vertex (1/2) tr(eps)(mu + k^2/4) drho: one phonon residue "
           f"= n_tot eps_p/E_p (error {res_err:.1e}). Scalar sector carries the "
           "phonon pole at tree level, static kernel (n_tot/8)(mu + k^2/4)(tr eps)^2")

    report(f"RESULT: c_T = {cT:.6f} = {cT / c_ph:.6f} c_ph = {cT / c_rel:.6f} c_rel; "
           f"A > 0, B > 0, lambda = {lam_h:.6f}.")
    report(f"  |c_T/c - 1| = {cT / c_ph - 1:.3f} on the phonon cone, "
           f"{cT / c_rel - 1:.3f} on the relative cone; GW170817 needs < 1e-15. "
           "lambda = 1/8 is not the GR value 1.")
    return dict(mu=mu, r=r, c_ph=c_ph, c_rel=c_rel, A=A_tt, B=B_tt, c_T=cT,
                cT_over_cph=cT / c_ph, cT_over_crel=cT / c_rel, channels=coef,
                phonon_ratio=ratio_ph, lam=lam_h, kinetic=kin, uv_powers=tails,
                cutoff_shares=cut, vertex_check=hf, A_spectral=A_spec,
                B_direct=B_rich, scalar_kinetic_factor=scalar_kin)


if __name__ == "__main__":
    run()
