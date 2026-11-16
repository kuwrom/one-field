"""The carrier level as a relativistic field theory: one cone in the UV.

THE LEVEL. Three complex scalars phi_a (a in Z3), hbar = 1, c explicit,
    L = sum_a |d_t phi_a|^2 - c^2 |grad phi_a|^2 - m^2 c^4 |phi_a|^2 - V,
    V = (G0/2)(sum_a |phi_a|^2)^2 + G1 sum_k |phi_k|^2 Re(phi*_{k+1} phi_{k+2}),
symmetric under global U(1), Z3 and the S3 reflection (alpha = 0).

PREMISES. P1 locality (the framework's axiom: a density in the fields and
finitely many derivatives at one point). P2 the two derivative truncation.
P3 isotropy. P4 second order in time (a carrier has both frequency signs;
envelope_universality.py uses such a carrier). P5 one kinetic
normalization for the three components (the NLS's own diagonal kinetic
term; Z3 x S3 allow circulant kinetic matrices with c_0 != c_1 = c_2).
P6 the two invariant V (the NLS's own posit, carried over unchanged).
Given P1-P5, Lorentz invariance is not a further posit. With U(1), a
first order term i b (phi* d_t phi - c.c.) is a chemical potential that a
time dependent phase removes. And K|phi'|^2 - C|grad phi|^2 - M|phi|^2 - V
is boost invariant with c^2 = C/K. P1 supplies locality only.

1. NR LIMIT (exact). phi_a = exp(-i m c^2 t) psi_a / (c sqrt(2m)) gives
    L = (i/2)(psi* psi' - c.c.) - |grad psi|^2/(2m) - w(psi) + |psi'|^2/(2 m c^2),
    w = (g0/2) rho^2 + g1 sum_k rho_k Re(psi*_{k+1} psi_{k+2}),
    g_i = G_i / (4 m^2 c^4)   (g_i = G_i/(4 m^2) at c = 1),  g1/g0 = G1/G0,
    i psi' = -grad^2 psi/(2m) + dw/dpsi* + psi''/(2 m c^2).
V is U(1) invariant, so V(exp(-i theta) phi) = V(phi) identically: no
harmonic, no averaging, no generated invariant. The one correction is
psi''/(2mc^2). A REAL carrier u_a = (psi_a e^{-i theta} + c.c.)/sqrt(2m)
needs the rotating wave average, which sends a real quartic T_abcd u^4 to
the fully symmetric tensor 6 T_abcd psi_a psi_b psi*_c psi*_d/(4m^2):
    <(sum u^2)^2>          = [4 rho^2 + 2 |sum psi_a^2|^2] / (4m^2),
    <sum_k u_k^2 u_k+1 u_k+2> = [4 W1 + 2 sum_k Re(psi_k^2 psi*_k+1 psi*_k+2)] / (4m^2).
The U(1) x S3 x C quartics form a 6-dimensional space. The real carrier
image is its 4-dimensional fully symmetric part, and the NLS plane
{rho^2, W1} is complementary to it: no real quartic carrier yields the
Z3-NLS. Complex carrier fields are forced by the NLS, not chosen.

2. CONDENSATE. psi = Psi exp(-i mu_NR t), Psi_a = sqrt(rho) w^a (the clock),
mu_tot = m c^2 + mu_NR, mu_tot^2 = m^2 c^4 + 2 m c^2 mu_w, mu_w = (3g0 - g1) rho.
In the aligned frame (gravity_tensor.py: A = diag(mu_w, 2r, r/2), B_00 =
mu_w, B_12 = r/2) each channel (a, d, b) obeys, with s = c^2 p^2,
    [(omega + mu_tot)^2 - s - M_a^2] [(omega - mu_tot)^2 - s - M_d^2] = beta^2,
    M_x^2 = mu_tot^2 + 2 m c^2 x = m^2 c^4 + 2 m c^2 (mu_w + x),  beta = 2 m c^2 b:
a particle Klein-Gordon shell and an antiparticle shell, split by 2 mu_tot,
coupled by the anomalous term. Channels (mu_w, mu_w, mu_w), (2r, r/2, r/2),
(r/2, 2r, r/2): 6 positive branches, 3 particle like (E = omega + mu_tot),
3 antiparticle like near omega = 2 m c^2 (E = omega - mu_tot).
 (i) phonon: (s - w^2)(s - w^2 + 2(mu_tot^2 - m^2c^4)) = 4 mu_tot^2 w^2,
     c_s^2 = c^2 (mu_tot^2 - m^2c^4)/(3 mu_tot^2 - m^2c^4) = c^2 kappa/(1 + 3 kappa),
     kappa = mu_w/(m c^2): c_s = c_ph/sqrt(1 + 3 kappa) -> sqrt(mu_w/m) (NR),
     -> c/sqrt(3) (ultra relativistic); Fagnocchi et al. arXiv:1001.1044
     eq. (46), c_s^2 = c^2 b/(1+b), with b = (mu_tot^2 - m^2c^4)/(2 mu_tot^2).
 (ii) as c -> infinity the channel equation becomes (eps + a - w)(eps + d + w)
     = b^2, which is one_cone.py's Omega^2 = E0^2 + (5r/4) p^2 + p^4/4,
     Omega = omega -/+ 3r/4. The residual is O(kappa) (checked below).
 (iii) at p >> mc every branch reads E^2 = c^2 p^2 + M_x^2 -/+ beta^2/(4 mu_tot c p)
     + ... (particle -, antiparticle +): one cone, with medium masses M_x
     (Hartree plus exchange), crossover |v/c - 1| < delta at p ~ M_x/(c sqrt(2 delta)).

3. GATE C AT THIS LEVEL. The principal symbol is (omega^2 - c^2 p^2) I_3 in
every background (V has no derivatives), so fronts of any disturbance,
knot or quasiparticle, move at exactly c: one causal cone. A localized,
stationary knot's excess four momentum (von Laue) transforms under carrier
boosts together with its condensate. A knot comoving with its condensate
has E = gamma M c^2 exactly. A knot moving THROUGH the condensate is a
different solution, and the condensate sets its dispersion. c_s and c_rel are internal
speeds of the medium, not the light cone. In the empty field the
carrier's cone is the only cone (interference/gravity.py, gate C). In a
condensate the only gapless mode is the phonon, capped at c/sqrt(3) < c.
If light were the phonon, every carrier quantum at p >> mc would outrun it.

Units in the code: hbar = m = 1, g0 rho = 1, g1 = g0/sqrt(2) (NLS units of
gravity_tensor.clock_background). Nu = mu_NR/(m c^2) sets c^2 = mu_w/kappa,
kappa = nu + nu^2/2. Nothing is fitted. Nu is a scanned input.

Usage: python3 tests/probes/relativistic_carrier.py     (a few seconds)
"""

import itertools
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from gravity_tensor import clock_background, clock_frequencies  # noqa: E402
from nls_soliton import interaction_matrix  # noqa: E402

NUS = (0.01, 0.1)
DELTAS = (1e-2, 1e-4)


def potential(f, G0, G1):
    """The two invariant quartic V (equally w with g in place of G)."""
    rho = np.abs(f)**2
    return G0 / 2 * rho.sum()**2 + G1 * sum(
        rho[k] * np.real(np.conj(f[(k + 1) % 3]) * f[(k + 2) % 3]) for k in range(3))


def nr_map_check(rng):
    """Pointwise identity L_rel(phi) = L_NLS(psi) + |psi'|^2/(2 m c^2)."""
    worst = 0.0
    for _ in range(20):
        m, c, G0, G1, t = rng.uniform(0.3, 3, 5)
        psi, dpsi = (rng.normal(size=3) + 1j * rng.normal(size=3) for _ in "ab")
        grad = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        z = np.exp(-1j * m * c * c * t) / (c * math.sqrt(2 * m))
        phi, dphi, gphi = z * psi, z * (dpsi - 1j * m * c * c * psi), z * grad
        l_rel = (np.sum(abs(dphi)**2) - c * c * np.sum(abs(gphi)**2)
                 - m * m * c**4 * np.sum(abs(phi)**2) - potential(phi, G0, G1))
        g0, g1 = G0 / (4 * m * m * c**4), G1 / (4 * m * m * c**4)
        l_nls = (-np.imag(np.vdot(psi, dpsi)) - np.sum(abs(grad)**2) / (2 * m)
                 - potential(psi, g0, g1))
        worst = max(worst, abs(l_rel - l_nls - np.sum(abs(dpsi)**2) / (2 * m * c * c)))
    return worst


def wirtinger_check(rng, g0=1.0, g1=1 / math.sqrt(2), h=1e-6):
    """dw/dpsi* equals the repository's interaction_matrix(psi) psi."""
    psi = rng.normal(size=3) + 1j * rng.normal(size=3)
    grad = np.array([(potential(psi + h * e, g0, g1) - potential(psi - h * e, g0, g1)
                      + 1j * (potential(psi + 1j * h * e, g0, g1)
                              - potential(psi - 1j * h * e, g0, g1))) / (4 * h)
                     for e in np.eye(3)])
    return np.max(abs(grad - interaction_matrix(psi, g0, g1) @ psi))


def rotating_wave(rng, m=1.3, n=32):
    """Complex carrier: V(e^{-i th} phi) is constant. Real carrier: the
    average of its quartics carries the extra invariants."""
    psi = rng.normal(size=3) + 1j * rng.normal(size=3)
    th = 2 * np.pi * np.arange(n) / n
    vc = [potential(np.exp(-1j * x) * psi, 1.0, 0.7) for x in th]
    u = [2 * np.real(psi * np.exp(-1j * x)) / math.sqrt(2 * m) for x in th]
    r0 = np.mean([np.sum(v * v)**2 for v in u])
    r1 = np.mean([sum(v[k]**2 * v[(k + 1) % 3] * v[(k + 2) % 3] for k in range(3))
                  for v in u])
    rho = np.abs(psi)**2
    w1 = sum(rho[k] * np.real(np.conj(psi[(k + 1) % 3]) * psi[(k + 2) % 3]) for k in range(3))
    w2 = sum(np.real(psi[k]**2 * np.conj(psi[(k + 1) % 3] * psi[(k + 2) % 3])) for k in range(3))
    e0 = abs(r0 - (4 * rho.sum()**2 + 2 * abs(np.sum(psi**2))**2) / (4 * m * m))
    e1 = abs(r1 - (4 * w1 + 2 * w2) / (4 * m * m))
    return float(np.ptp(vc)), float(max(e0, e1))


def invariant_spaces():
    """Dimensions of the U(1) x S3 x C quartic space, its fully symmetric
    (real carrier) part, and the rank of that part joined with the NLS plane."""
    def pair_sym(t):
        t = (t + t.transpose(1, 0, 2, 3)) / 2
        t = (t + t.transpose(0, 1, 3, 2)) / 2
        return (t + t.transpose(2, 3, 0, 1)) / 2

    def full_sym(t):
        return sum(t.transpose(q) for q in itertools.permutations(range(4))) / 24

    def s3_average(t):
        out = 0
        for q in itertools.permutations(range(3)):
            P = np.eye(3)[list(q)]
            out = out + np.einsum('ai,bj,ck,dl,ijkl->abcd', P, P, P, P, t)
        return out / 6

    raw = []
    for idx in itertools.product(range(3), repeat=4):
        t = np.zeros((3,) * 4)
        t[idx] = 1
        raw.append(s3_average(pair_sym(t)).ravel())
    _, s, vt = np.linalg.svd(np.array(raw))
    inv = vt[:int((s > 1e-10).sum())]
    sym = np.array([full_sym(v.reshape((3,) * 4)).ravel() for v in inv])
    _, s2, vt2 = np.linalg.svd(sym)
    sym = vt2[:int((s2 > 1e-10).sum())]
    d = np.eye(3)
    c_rho = pair_sym(np.einsum('ac,bd->abcd', d, d))
    c_w1 = np.zeros((3,) * 4)
    for k in range(3):
        c_w1[k, (k + 2) % 3, k, (k + 1) % 3] += 0.5
        c_w1[k, (k + 1) % 3, k, (k + 2) % 3] += 0.5
    c_w1 = pair_sym(c_w1)
    joined = np.linalg.matrix_rank(np.vstack([sym, c_rho.ravel(), c_w1.ravel()]), 1e-10)
    off_sym = np.array([(c - full_sym(c)).ravel() for c in (c_rho, c_w1)]).T
    return len(inv), len(sym), joined, float(np.linalg.svd(off_sym, compute_uv=False).min())


def parameters(nu, mu_w):
    kappa = nu + nu * nu / 2
    c2 = mu_w / kappa
    return kappa, c2, c2 * (1 + nu)          # kappa, c^2, mu_tot (m = 1)


def carrier_operator(A, B, p, c2, mu_tot):
    """Linearized quadratic problem [K - omega Lam - omega^2/(2c^2)] v = 0 as
    a 12 x 12 first order operator on (v, omega v). Raw component basis."""
    eye = np.eye(3)
    K = np.block([[A + p * p / 2 * eye, B], [B.conj(), A.conj() + p * p / 2 * eye]])
    Lam = mu_tot / c2 * np.diag([1.0, 1, 1, -1, -1, -1])
    Z = np.zeros((6, 6))
    return np.block([[Z, np.eye(6)], [2 * c2 * K, -2 * c2 * Lam]])


def channels(mu_w, r):
    return dict(phase=(mu_w, mu_w, mu_w), upper=(2 * r, r / 2, r / 2),
                lower=(r / 2, 2 * r, r / 2))


def channel_parts(w, p, c2, mu_tot, a, d, b):
    s = c2 * p * p
    P = (w + mu_tot)**2 - s - (mu_tot**2 + 2 * c2 * a)
    Q = (w - mu_tot)**2 - s - (mu_tot**2 + 2 * c2 * d)
    return P, Q, (2 * c2 * b)**2


def channel_roots(p, c2, mu_tot, a, d, b):
    """Positive roots (particle, antiparticle) of one channel, Newton polished.
    The phase channel (a = d = b) is even in omega: closed form, rationalized."""
    s = c2 * p * p
    if a == d == b:
        D = 3 * mu_tot**2 - c2 * c2
        root = math.sqrt(4 * mu_tot**2 * s + D * D)
        low = s * (s + 2 * mu_tot**2 - 2 * c2 * c2) / (s + D + root)
        return np.sqrt([low, s + D + root])
    P = np.array([1.0, 2 * mu_tot, -s - 2 * c2 * a])
    Q = np.array([1.0, -2 * mu_tot, -s - 2 * c2 * d])
    poly = np.polymul(P, Q)
    poly[-1] -= (2 * c2 * b)**2
    roots = np.sort(np.roots(poly).real)[2:]
    for _ in range(4):
        Pv, Qv, _b2 = channel_parts(roots, p, c2, mu_tot, a, d, b)
        roots = roots - (Pv * Qv - _b2) / (2 * (roots + mu_tot) * Qv + 2 * (roots - mu_tot) * Pv)
    return roots


# branch -> (channel, antiparticle?). An antiparticle carries the HF entry d
# of its channel, so antilower (HF r/2) is the conjugate of lower.
BRANCHES = {"phonon": ("phase", 0), "lower": ("lower", 0), "upper": ("upper", 0),
            "anti-phonon": ("phase", 1), "anti-lower": ("upper", 1),
            "anti-upper": ("lower", 1)}


def medium_mass2(name, c2, mu_tot, chans):
    ch, anti = BRANCHES[name]
    a, d, _ = chans[ch]
    return mu_tot**2 + 2 * c2 * (d if anti else a)


def branch(name, p, c2, mu_tot, chans):
    """Lab energy E and group velocity of one of the six branches."""
    ch, anti = BRANCHES[name]
    a, d, b = chans[ch]
    w = channel_roots(p, c2, mu_tot, a, d, b)[anti]
    P, Q, _ = channel_parts(w, p, c2, mu_tot, a, d, b)
    vg = (P + Q) * 2 * c2 * p / (2 * (w + mu_tot) * Q + 2 * (w - mu_tot) * P)
    return (w - mu_tot if anti else w + mu_tot), vg


def crossover(name, c2, mu_tot, chans, delta, kind):
    c = math.sqrt(c2)

    def dev(p):
        E, vg = branch(name, p, c2, mu_tot, chans)
        return abs(E / (c * p) - 1 if kind == "phase" else 1 - vg / c) - delta

    grid = c * np.logspace(-2, 3, 300)
    vals = np.array([dev(p) + delta for p in grid])
    i = int(np.nonzero(vals >= delta)[0].max())
    return brentq(dev, grid[i], grid[i + 1], xtol=1e-12, rtol=1e-12) / c


def run(report=print):
    rng = np.random.default_rng(7)
    _, mu_w, r, A, B = clock_background()
    chans = channels(mu_w, r)

    # 1. NR map and the complex carrier argument.
    e_map, e_wirt = nr_map_check(rng), wirtinger_check(rng)
    spread_complex, e_real = rotating_wave(rng)
    n_inv, n_sym, joined, off = invariant_spaces()
    assert e_map < 1e-11 and e_wirt < 1e-8 and spread_complex < 1e-12 and e_real < 1e-12
    assert (n_inv, n_sym, joined) == (6, 4, 6) and off > 0.3
    report("1. phi = e^{-imc^2 t} psi/(c sqrt(2m)): L_rel = L_NLS + |psi'|^2/(2mc^2) "
           f"pointwise (error {e_map:.1e}), g_i = G_i/(4 m^2 c^4). Dw/dpsi* = "
           f"interaction_matrix psi ({e_wirt:.1e})")
    report(f"   complex carrier: V(e^{{-i th}} phi) spread {spread_complex:.1e} (no "
           f"harmonics). Real carrier averages match the extra invariant formulas ({e_real:.1e})")
    report(f"   U(1) x S3 x C quartics: {n_inv}. Real carrier image {n_sym} (fully "
           f"symmetric). With the NLS plane rank {joined}: complementary, least "
           f"off symmetric singular value {off:.4f}")

    # 2. The six branch spectrum. Operator versus channel equations.
    out = dict(nr_map_error=e_map, invariants=n_inv, real_image=n_sym, per_nu={})
    for nu in NUS:
        kappa, c2, mu_tot = parameters(nu, mu_w)
        c = math.sqrt(c2)
        err = 0.0
        for p in (0.05, 1.0, 5.0, 3 * c):
            eig = np.linalg.eigvals(carrier_operator(A, B, p, c2, mu_tot))
            assert max(abs(eig.imag)) < 1e-6 * c2
            pos = np.sort(eig.real)[6:]
            chan = np.sort(np.concatenate([channel_roots(p, c2, mu_tot, *v)
                                           for v in chans.values()]))
            err = max(err, float(np.max(abs(pos - chan) / pos)))
        assert err < 1e-9
        gaps = {n: branch(n, 1e-9 * c, c2, mu_tot, chans)[0] for n in BRANCHES}
        # (i) phonon speed: channel slope, closed form, Fagnocchi eq. (46).
        p0 = 1e-4
        eig0 = np.linalg.eigvals(carrier_operator(A, B, p0, c2, mu_tot)).real
        slope = np.min(eig0[eig0 > 0]) / p0
        cs = c * math.sqrt(kappa / (1 + 3 * kappa))
        cs_fag = c * math.sqrt((mu_tot**2 - c2**2) / (3 * mu_tot**2 - c2**2))
        assert abs(slope / cs - 1) < 1e-6 and abs(cs_fag / cs - 1) < 1e-12
        # (iii) every branch -> cp, with medium masses M_x.
        cross = {n: {(k, dl): crossover(n, c2, mu_tot, chans, dl, k)
                     for k in ("phase", "group") for dl in DELTAS} for n in BRANCHES}
        big = 1e3 * c
        asym, vmax = {}, 0.0
        for n in BRANCHES:
            M2 = medium_mass2(n, c2, mu_tot, chans)
            E, _ = branch(n, big, c2, mu_tot, chans)
            asym[n] = 2 * (c * big)**2 * (E / (c * big) - 1) / M2
            vmax = max(vmax, max(branch(n, q, c2, mu_tot, chans)[1] / c
                                 for q in c * np.logspace(-3, 3, 400)))
        assert all(abs(v - 1) < 1e-3 for v in asym.values())
        assert vmax < 1
        report(f"2. nu = mu_NR/(mc^2) = {nu}: kappa = {kappa:.6f}, c = {c:.6f}, "
               f"mu_tot = {mu_tot:.6f}, 1/xi = {math.sqrt(2 * mu_w):.4f} = "
               f"{math.sqrt(2 * kappa):.4f} mc. Operator vs channel roots {err:.1e}")
        report("   p -> 0 lab energies E (particle: omega + mu_tot, antiparticle: "
               "omega - mu_tot), units of mc^2: "
               + ", ".join(f"{n} {gaps[n] / c2:.6f}" for n in BRANCHES))
        report(f"   (i) c_s = c sqrt(kappa/(1+3kappa)) = {cs:.9f} = {cs / c:.6f} c "
               f"(operator slope {slope:.9f}); NLS c_ph = {math.sqrt(mu_w):.6f}, "
               f"c_rel = {math.sqrt(1.25 * r):.6f} = {math.sqrt(1.25 * r) / c:.6f} c")
        report(f"   (iii) 2 c^2p^2 (E/cp - 1)/M_x^2 at p = 1000 mc: "
               + ", ".join(f"{n} {v:.6f}" for n, v in asym.items())
               + f". Max v_g/c over 1e-3..1e3 mc = {vmax:.9f}")
        estimate = {}
        for n in BRANCHES:
            row = cross[n]
            estimate[n] = math.sqrt(medium_mass2(n, c2, mu_tot, chans)) / c2
            assert abs(row[("phase", 1e-4)] * math.sqrt(2e-4) / estimate[n] - 1) < 1e-2
            report(f"     {n:11s} |v_ph/c - 1| < 1e-2, 1e-4 beyond p = "
                   f"{row[('phase', 1e-2)]:8.4f}, {row[('phase', 1e-4)]:9.4f} mc;"
                   f"  |v_g/c - 1|: {row[('group', 1e-2)]:8.4f}, "
                   f"{row[('group', 1e-4)]:9.4f} mc;  M_x/mc^2 = {estimate[n]:.6f}")
        out["per_nu"][nu] = dict(kappa=kappa, c=c, mu_tot=mu_tot, c_s=cs, c_s_over_c=cs / c,
                                 gaps=gaps, crossover=cross, asymptotic_ratio=asym,
                                 medium_mass=estimate,
                                 max_group_over_c=vmax, operator_error=err)

    # (ii) NR limit: particle branches -> clock_frequencies, residual O(kappa).
    resid = {}
    for nu in (1e-3, 1e-4):
        kappa, c2, mu_tot = parameters(nu, mu_w)
        worst = kg = 0.0
        for p in (0.03, 0.4, 1.1, 3.0):
            w = np.array([branch(n, p, c2, mu_tot, chans)[0] - mu_tot
                          for n in ("phonon", "lower", "upper")])
            worst = max(worst, float(np.max(abs(w / clock_frequencies(p, mu_w, r) - 1))))
            target = (21 * r * r / 16) + 1.25 * r * p * p + p**4 / 4
            for om in (w[1] + 0.75 * r, w[2] - 0.75 * r):
                kg = max(kg, abs(om * om / target - 1))
        resid[nu] = (worst, kg, kappa)
    scale = resid[1e-3][0] / resid[1e-4][0]
    assert resid[1e-4][0] < 3e-4 and resid[1e-4][1] < 1e-3 and 9 < scale < 11
    report(f"   (ii) NR limit vs one_cone/gravity_tensor branches: max rel. residual "
           f"{resid[1e-3][0]:.2e} (kappa {resid[1e-3][2]:.1e}), {resid[1e-4][0]:.2e} "
           f"(kappa {resid[1e-4][2]:.1e}), ratio {scale:.3f}: O(kappa); Omega^2 = "
           f"E0^2 + (5r/4)p^2 + p^4/4 residual {resid[1e-4][1]:.2e}")
    k_ur = parameters(10.0, mu_w)[0]
    report(f"   ultra relativistic check nu = 10: c_s^2/c^2 = {k_ur / (1 + 3 * k_ur):.6f} "
           "-> 1/3. The phonon never reaches the carrier cone")
    report("RESULT: given P1-P5 the carrier level is Lorentz invariant. Its NR limit is "
           "the Z3-NLS exactly up to psi''/(2mc^2). Complex fields are forced. All "
           "six branches end on E^2 = c^2p^2 + M^2 (one causal cone). c_s and c_rel are "
           "internal speeds of a condensate. In the empty field the carrier's cone is "
           "the only cone.")
    out.update(nr_residual=resid, nr_scaling=scale)
    return out


if __name__ == "__main__":
    run()
