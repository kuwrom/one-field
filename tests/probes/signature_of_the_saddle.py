"""The signature of the saddle: which constant metric carries a clock.

THE QUESTION. The carrier written without a background metric is the
theory of the determinant density sqrt|det h| / U (composite_metric.py
C1). That density is GL(4) covariant, so a constant saddle
h_vac = diag(sigma_0, ..., sigma_3) of any sign pattern solves the same
gap equation as any other: the form does not know a signature. The
signature is read off the saddle. This probe computes which saddle
carries a carrier that persists, and what the positive norm vacuum of
that carrier puts into the bilinear.

PREMISES. P1 the determinant form with U > 0 (composite_metric.py). P2
the fluctuation bilinear s = sum_a (d conj(phi_a) d phi_a + c.c.) of the
carrier about the saddle, to first order in s. P3 for S4 the free
complex field of mass m in its Fock vacuum, with a regulator that is
positive mode by mode but picks a frame, against the covariant
proper time weight of composite_pole.py.

STATEMENTS (the numbers below check each one).

  S1 (the carrier propagates on the inverse saddle). To first order in
     s, sqrt|det(h_vac + s)| / U = sqrt|det h_vac| / U x [1 + (1/2)
     tr(h_vac^-1 s)] - (sqrt|det h_vac| / U^2) dU: the kinetic form is
     (sqrt|det h_vac| / U) sum_mu (1/sigma_mu) |d_mu phi|^2 and the
     potential enters with a minus sign. Checked by finite differences
     on all 16 sign patterns.
  S2 (energy positivity picks the signature). With x^tau as the
     evolution parameter the canonical Hamiltonian density of S1 is
     bounded below iff sigma_tau > 0 > sigma_mu for every mu != tau
     (the inverse saddle has the signs of the saddle). Of the 64
     (pattern, tau) cases exactly 4 pass: the Lorentzian patterns with
     tau on the odd sign. In relativistic_carrier.py this is K > 0 and
     C > 0, the one place a sign is fixed, and it is fixed by the
     carrier having a ground state. U -> -U would flip every sign at
     once, so the overall sign is the potential's.
  S3 (well posedness, Hadamard). On the Euclidean and ultrahyperbolic
     saddles, and on a Lorentzian saddle evolved along one of its
     majority directions, the linear modes exp(i k.x) have growth rates
     unbounded in |k|: max Im omega over |k| <= K grows like K (at
     K = 2, 10, 100 it is about 1.7, 9.9, 99.9 in units of m), and no
     initial data fix a past.  The Lorentzian saddle of the opposite
     overall sign, evolved along its odd direction, is hyperbolic and
     well posed. Its growth is sqrt(m^2 - k^2), bounded by m, in a
     tachyonic band |k| < m. S2 rejects it because the carrier has no
     ground state there.  So the type of the equation fixes the Lorentzian
     signature (only a hyperbolic equation carries a record of its
     past), and the ground state fixes its orientation.  On the
     passing saddle Im omega = 0 for every k. A clock, exp(-i Omega t) f(r) with Omega < m, needs the
     hyperbolic dispersion, so the lump of lump_in_empty_field.py lives
     on that saddle alone.
  S4 (what the positive norm vacuum puts into the bilinear). With a
     frame fixed positive regulator f(|k|/Lambda) the Fock vacuum gives
     <d_mu conj(phi) d_nu phi> = diag(A, B, B, B), both positive. Split
     into the trace part and the traceless frame part,
         h = C eta + D (u u - eta/4),  eta = diag(+,-,-,-), u = (1,0,0,0):
     C = (A - 3B)/2 = m^2 <|phi|^2> / 2 by the virial identity
     omega^2 - k^2 = m^2 mode by mode. C is positive, two powers of
     Lambda below D (it carries the m^2 Lambda^2 of <|phi|^2>, not
     Lambda^4), and has the orientation that passes S2. D grows as
     Lambda^4 and belongs to the frame. D gives a
     cutoff vacuum energy the equation of state p = rho/3. Under the
     covariant proper time weight the frame part is absent by
     construction. The Euclidean value (I1/4) delta of composite_pole.py
     K2 splits as contact minus field:
     I1 = [exp(-m^2 t0)/t0^2 - m^2 t0 E2(m^2 t0)/t0] / (16 pi^2). The
     field part -m^2 <phi^2>_E continues to +m^2 <|phi|^2> eta, the
     same tensor as the frame computation's C eta. So the saddle's
     signature and orientation are the carrier's own. The vacuum term
     tuning of K2 absorbs the contact part.
  S5 (one cone). On the passing saddle h_vac = diag(a_0, -a_1, -a_1,
     -a_1) the characteristic cone of S1's kinetic form is c^2 = C/K =
     a_0/a_1, the cone of g = h/U, independent of U and of the overall
     scale: the composite metric's cone is the carrier's cone of gate C
     (one_cone.py, relativistic_carrier.py).

RESULT. Of the 64 (pattern, tau) cases, exactly 4 have a Hamiltonian
bounded below: the (1,3) saddles with tau on the odd sign, where Im
omega = 0 at every sampled k. Every other case grows without bound in
|k| or has no ground state. On the passing saddle the positive norm
vacuum fills the trace part of the bilinear with C eta, C > 0, and the
kinetic cone is the cone of g = h/U (c^2 = 4 at both checked points).

Usage: python3 tests/probes/signature_of_the_saddle.py     (a few seconds)
"""

import itertools
import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import composite_pole as cp  # noqa: E402

PATTERNS = list(itertools.product((1.0, -1.0), repeat=4))


def gram(grad):
    """s_mu nu = sum_a (d_mu conj(phi_a) d_nu phi_a + c.c.) for grad[a, mu]."""
    grad = np.asarray(grad, dtype=complex)
    return 2 * (grad.conj().T @ grad).real


def density(h, U):
    return math.sqrt(abs(np.linalg.det(h))) / U


def first_order(rng, U=1.3, eps=1e-6):
    """S1: finite difference first order coefficients on every pattern."""
    worst_kin = worst_pot = 0.0
    for pat in PATTERNS:
        hv = np.diag(pat)
        grad = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
        s = gram(grad)
        d0 = density(hv, U)
        fd = (density(hv + eps * s, U) - density(hv - eps * s, U)) / (2 * eps)
        kin = d0 * sum(pat[m] * np.sum(abs(grad[:, m])**2) for m in range(4))
        worst_kin = max(worst_kin, abs(fd / kin - 1))
        dU = 0.37
        fdU = (density(hv, U + eps * dU) - density(hv, U - eps * dU)) / (2 * eps)
        worst_pot = max(worst_pot, abs(fdU / (-d0 * dU / U) - 1))
    return worst_kin, worst_pot


def hamiltonian_count(U=1.3, M2=1.0, samples=200, seed=1):
    """S2: (pattern, tau) cases whose Hamiltonian density is bounded below."""
    rng = np.random.default_rng(seed)
    passing, unbounded = [], 0
    for pat in PATTERNS:
        for tau in range(4):
            # L = k_tau |d_tau phi|^2 + sum_i k_i |d_i phi|^2 - m |phi|^2, k_mu = sigma_mu
            # (the common factor sqrt|det h_vac| / U is positive), m = M2 / U > 0
            k = np.array(pat)
            kt, ks = k[tau], np.delete(k, tau)
            m = M2 / U
            # H = |pi|^2 / k_tau - sum_i k_i |d_i phi|^2 + m |phi|^2
            bounded = kt > 0 and np.all(ks < 0)
            if bounded:
                passing.append((pat, tau))
            # witness: a mode with H < 0 that scales as its amplitude squared
            lowest = np.inf
            for _ in range(samples):
                pi, grads, amp = rng.normal(size=2), rng.normal(size=(3, 2)), rng.normal(size=2)
                H = (pi @ pi) / kt - sum(ks[i] * (grads[i] @ grads[i]) for i in range(3)) \
                    + m * (amp @ amp)
                lowest = min(lowest, H)
            if (lowest < 0) != (not bounded):
                raise RuntimeError("positivity witness disagrees with the sign rule")
            unbounded += lowest < 0
    return passing, unbounded


def growth(pat, tau, K, M2=1.0, n=60):
    """S3: max Im omega over |k| <= K on the saddle pat, evolution along tau."""
    k = np.array(pat)
    kt, ks = k[tau], np.delete(k, tau)
    grid = np.linspace(-K, K, n)
    worst = 0.0
    for kx in grid:
        for ky in grid:
            for kz in grid:
                q = np.array([kx, ky, kz])
                if q @ q > K * K:
                    continue
                # sigma_tau omega^2 + sum_i sigma_i k_i^2 = M2  (from sum_mu sigma_mu k_mu^2 = M2)
                w2 = (M2 - ks @ (q * q)) / kt
                worst = max(worst, math.sqrt(-w2) if w2 < 0 else 0.0)
    return worst


def vacuum_bilinear(m, Lam, kind):
    """S4: A = <|d_t phi|^2>, B = <|d_x phi|^2>, <|phi|^2> in the Fock vacuum."""
    if kind == "sharp":
        f = lambda k: 1.0 if k <= Lam else 0.0
        top = Lam
    else:
        f = lambda k: math.exp(-(k / Lam)**2)
        top = 12 * Lam
    meas = 1.0 / (2 * math.pi**2)                       # d^3k/(2 pi)^3 = k^2 dk meas
    w = lambda k: math.sqrt(k * k + m * m)
    A = quad(lambda k: meas * k * k * w(k) * f(k) / 2, 0, top, limit=400)[0]
    B = quad(lambda k: meas * k**4 * f(k) / (6 * w(k)), 0, top, limit=400)[0]
    phi2 = quad(lambda k: meas * k * k * f(k) / (2 * w(k)), 0, top, limit=400)[0]
    # the trace part from one integrand: (omega^2 - k^2) / (2 omega) = m^2 / (2 omega)
    AB = quad(lambda k: meas * k * k * f(k) * (w(k) / 2 - k * k / (2 * w(k))), 0, top,
              limit=400)[0]
    h = 2 * np.diag([A, B, B, B])                       # the c.c. factor
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    C = AB / 2                                          # trace part, tr(eta h)/4
    D = 4.0 / 3.0 * (h[0, 0] - C)                       # traceless frame part
    frame = C * eta + D * (np.diag([1.0, 0, 0, 0]) - eta / 4)
    return dict(A=A, B=B, phi2=phi2, C=C, D=D, virial=(AB - m * m * phi2) / AB,
                split=(A - 3 * B - AB) / A, recon=np.max(abs(frame - h)) / A, frame=frame)


def euclidean_split(m2, t0):
    """S4: I1 of composite_pole K2 as contact minus field part."""
    J, _, _ = cp.proper_time(m2, t0)
    I1 = 2 * cp.C16 * float(J(2, 0, 0.0))
    contact = math.exp(-m2 * t0) / (16 * math.pi**2 * t0**2)
    phi2_E = expn(2, m2 * t0) / (16 * math.pi**2 * t0)
    return I1, contact, m2 * phi2_E


def cone(a0, a1, U, rng, eps=1e-6):
    """S5: c^2 = C/K of the first order kinetic form on diag(a0, -a1, -a1, -a1)."""
    hv = np.diag([a0, -a1, -a1, -a1])
    d0 = density(hv, U)

    def coeff(mu):
        grad = np.zeros((3, 4), dtype=complex)
        grad[0, mu] = 1.0
        s = gram(grad)
        return (density(hv + eps * s, U) - density(hv - eps * s, U)) / (2 * eps)
    K, C = coeff(0), -coeff(1)
    g = hv / U
    null = np.array([1.0, math.sqrt(C / K), 0.0, 0.0])    # (t, x) on the cone of the kinetic form
    return C / K, abs(null @ g @ null) / abs(g[0, 0]), d0


def run(report=print):
    rng = np.random.default_rng(0)
    kin, pot = first_order(rng)
    assert kin < 1e-6 and pot < 1e-6
    report(f"S1: first order kinetic coefficient sqrt|det h_vac|/U x sigma_mu on all 16 "
           f"patterns, worst relative deviation {kin:.1e}. Potential term -sqrt|det|/U^2 dU, "
           f"{pot:.1e}.")

    passing, unbounded = hamiltonian_count()
    assert len(passing) == 4 and unbounded == 60
    assert all(pat[tau] > 0 and sum(pat) == -2 for pat, tau in passing)
    report(f"S2: of 64 (pattern, tau) cases {len(passing)} have a Hamiltonian bounded below: "
           + ", ".join(f"{tuple(int(x) for x in pat)} along x^{tau}" for pat, tau in passing)
           + f". The other {unbounded} are unbounded below (random modes with H < 0).")

    cases = {"(+,-,-,-) along the +": ((1, -1, -1, -1), 0),
             "(+,-,-,-) along a -": ((1, -1, -1, -1), 1),
             "(+,+,+,+)": ((1, 1, 1, 1), 0),
             "(+,+,-,-) along a +": ((1, 1, -1, -1), 0),
             "(-,+,+,+) along the -": ((-1, 1, 1, 1), 0)}
    rates = {}
    for name, (pat, tau) in cases.items():
        rates[name] = [growth(pat, tau, K) for K in (2.0, 10.0, 100.0)]
    assert max(rates["(+,-,-,-) along the +"]) == 0.0
    for name in list(cases)[1:4]:
        r = rates[name]
        assert r[2] > 9 * r[1] > 4 * r[0] > 0
    # growth = sqrt(m^2 - k^2) <= m, attained at k = 0. The |k| <= 100 grid steps over the band
    r_t = rates["(-,+,+,+) along the -"]
    assert all(x <= 1.0 for x in r_t) and r_t[0] > 0.99
    report("S3: max Im omega over |k| <= 2, 10, 100: "
           + "; ".join(f"{n}: " + ", ".join(f"{x:.2f}" for x in r) for n, r in rates.items())
           + ". Growth unbounded in |k| on the nonhyperbolic readings (Hadamard), bounded by m "
           "on the opposite sign Lorentzian saddle (well posed, no ground state, S2), zero on "
           "the passing one.")

    m = 1.0
    vac = {}
    for kind in ("sharp", "smooth"):
        for Lam in (10.0, 100.0):
            v = vacuum_bilinear(m, Lam, kind)
            assert abs(v["virial"]) < 1e-8 and abs(v["split"]) < 1e-7 and v["recon"] < 1e-7
            assert v["C"] > 0 and v["D"] > 0
            vac[(kind, Lam)] = v
    sharp = (vac[("sharp", 100.0)]["D"] / vac[("sharp", 10.0)]["D"],
             vac[("sharp", 100.0)]["C"] / vac[("sharp", 10.0)]["C"])
    assert 0.9e4 < sharp[0] < 1.1e4 and 90 < sharp[1] < 110
    for (kind, Lam), v in vac.items():
        report(f"S4 ({kind} cutoff, Lambda = {Lam:g} m): A = {v['A']:.4g}, B = {v['B']:.4g}; "
               f"C = {v['C']:.4f} = m^2 <|phi|^2>/2 (virial residual {v['virial']:.1e}), "
               f"D = {v['D']:.4g}.")
    report(f"S4: from Lambda = 10 m to 100 m the frame part D grows x{sharp[0]:.0f} "
           f"(Lambda^4), the trace part C x{sharp[1]:.1f} (the m^2 Lambda^2 of <|phi|^2>). "
           "C > 0 on eta = diag(+,-,-,-): the orientation that passes S2.")
    I1, contact, field = euclidean_split(1.0, 1e-2)
    assert abs(I1 - (contact - field)) < 1e-12 * abs(I1)
    report(f"S4 (proper time, t0 = 1e-2): I1 = {I1:.4f} = contact {contact:.4f} - field "
           f"{field:.4f}. The field part m^2 <phi^2>_E continues to the C eta above, the "
           "contact part is the regulator's and is what K2's vacuum term tuning absorbs; "
           "no frame part exists under a covariant weight.")

    c2, on_cone, _ = cone(2.0, 0.5, 1.7, rng)
    c2b, on_cone_b, _ = cone(6.0, 1.5, 0.4, rng)
    assert abs(c2 - 4.0) < 1e-6 and abs(c2b - 4.0) < 1e-6
    assert on_cone < 1e-6 and on_cone_b < 1e-6
    report(f"S5: on diag(a0, -a1, -a1, -a1) the kinetic form's cone is c^2 = C/K = {c2:.6f} "
           f"= a0/a1 at (a0, a1, U) = (2, 0.5, 1.7) and {c2b:.6f} at (6, 1.5, 0.4). That "
           f"vector is null for g = h/U (|g(v,v)|/|g_00| = {on_cone:.1e}). One cone.")
    report("RESULT: the determinant form has saddles of every signature. The carrier "
           "has a ground state and a well posed evolution on the (1,3) saddles alone, "
           "oriented with the potential's sign, and its own positive norm vacuum fills "
           "the bilinear with that tensor: the signature is the carrier's, not a posit.")
    return dict(first_order=(kin, pot), passing=passing, unbounded=unbounded, rates=rates,
                vacuum=vac, euclidean=(I1, contact, field), cone=(c2, c2b, on_cone, on_cone_b))


if __name__ == "__main__":
    run()
