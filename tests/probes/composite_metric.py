"""The carrier written without a background metric: the composite metric.

THE QUESTION. The carrier level (relativistic_carrier.py) contracts its
two derivatives with a flat eta. Written with no metric at all, the same
fields carry a symmetric tensor of their own, the gradient bilinear
    h_mu nu = sum_a (d_mu conj(phi_a) d_nu phi_a + c.c.),
and the action is a density built from h and the potential U alone.

PREMISES. P1 the carrier fields of relativistic_carrier.py, three complex
scalars, six real components, U = m^2 |phi|^2 + V. P2 no background
metric: the action is built from the fields and their first derivatives.
P3 for C4 only: one loop, free fields, a covariant Pauli-Villars
regulator standing for the physical one.

STATEMENTS (the numbers below check each one).

  C1 (elimination). For L = sqrt|g| (g^{mu nu} h_mu nu / 2 - U) in D = 4,
     stationarity in g gives g = h / U and L = sqrt|det h| / U, and the
     energy momentum tensor T_mu nu = h_mu nu - g_mu nu (g^{ab} h_ab / 2
     - U) vanishes identically (Carone, Erlich, Vaman arXiv:1610.08521,
     Eqs. 2.1, 2.5, 2.6; without clock and ruler fields Carone
     arXiv:1911.13166 Eq. 4.1 and Batz, Erlich, Mrini arXiv:2011.06541
     Eqs. 1.1, 2.1, 2.2). A metric with no action of its own and no
     metric at all are one theory.
  C2 (traceless form). With det g held fixed the same step gives
     g = h / |det h|^(1/4) and L = 2 |det h|^(1/4) - U: homogeneous of
     degree two in the gradients, invariant under volume preserving
     maps, and blind to a constant in U.
  C3 (matter). About h = c eta + s, s the gradient bilinear of the
     fluctuations, the full form is V0 + eta.s/2 - dV - T.Pi.T/(4 V0)
     + O(1/V0^2), T the flat energy momentum tensor and
     Pi = (eta eta + eta eta - eta eta)/2 (arXiv:1610.08521 Eqs. 3.5,
     3.8). The traceless form is 2c + eta.s/2 - U - (s.s - (eta.s)^2/4)
     /(4c): the free carrier plus one contact term, of the traceless
     part alone in the second case.
  C4 (the tensor term at one loop). For a traceless transverse h the
     induced quadratic action is W2 = [w0 + kappa q^2 + ...] h.h. With
     I1 = <d phi . d phi>, I2 = int p^4 G^2 and rho the vacuum term,
     the regulator gives I2 = 3 I1 and I1 = -4 rho, so that
         full form:        w0 = -rho_total / 4, the vacuum term times
                           the second order part of sqrt g (g = 1 + h),
         traceless form:  w0 = I1/16 - I2/48 = 0 identically (det g = 1),
     and kappa = (1/4) (1/6) sum_j c_j M_j^2 ln M_j^2 / (32 pi^2) per
     real component: the heat kernel weight a1 = 1/6 is the stiffness
     of the tensor mode (arXiv:2011.06541 Eqs. 2.27-2.31).
  C5 (states). Every classical configuration has h positive
     semidefinite. The clock gives h = diag(h_00, 0, 0, 0), rank one; a
     texture gives spatial rank three where it sits. The traceless
     spatial part of h is twice the Q_ij of gravity_tensor.py.
  C6 (ledger). The venting ledger of interference/gravity.py sums the
     same a1 over its rows. With one proper time cutoff the induced
     vacuum term and stiffness obey Lambda G = -(pi/2) str(1) / Sigma^2,
     free of the cutoff.
  C7 (the form binds nothing). For a field of amplitude sigma turning
     at rate Omega with dV = m^2 sigma^2, the full form is
     L = V0 [sqrt(1 + 2 w x)/(1 + x) - 1], x = m^2 sigma^2/V0,
     w = Omega^2/m^2: negative at every x > 0 for Omega < m, so the
     form with a mass term holds no lump. Its contact term is the
     energy +rho^2/(8 V0) for matter of density rho, +3 rho^2/(16 c)
     in the traceless form: repulsive. The binding is the coupling's
     (self_bound_vacuum.py), not the metric's.

RESULT. Written without a background metric the carrier eliminates to
g = h/U with a vanishing energy momentum tensor, and the traceless
form is L = 2|det h|^(1/4) - U. At one loop under the Pauli-Villars
regulator the q = 0 tensor term is -rho/4 in the full form and zero
in the traceless form, and the q^2 stiffness is the heat kernel
weight a1 = 1/6 per real component. The form binds nothing, and both
contact energies are repulsive.

Usage: python3 tests/probes/composite_metric.py     (about half a minute
standalone: C6 solves the engine web once)
"""

import contextlib
import io
import math
import os
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import couplings, gravity, masses, root  # noqa: E402
from substrate3d import texture_state  # noqa: E402
from gravity_tensor import neutral_tensor  # noqa: E402
import carrier_sakharov as cs  # noqa: E402

ETA = np.diag([1.0, -1.0, -1.0, -1.0])
PV4 = ((1.0, 0.0), (-3.0, 1.0), (3.0, 2.0), (-1.0, 3.0))   # (c_j, s_j / L^2)
MEASURE = 1.0 / (8 * math.pi**2)                           # d^4p/(2pi)^4 = p^3 dp MEASURE


def gram(grad):
    """h_mu nu for grad[a, mu] = d_mu phi_a (complex components)."""
    grad = np.asarray(grad, dtype=complex)
    return 2 * (grad.conj().T @ grad).real


def lagrangian(g, h, U):
    return math.sqrt(abs(np.linalg.det(g))) * (0.5 * np.trace(np.linalg.inv(g) @ h) - U)


def stress(g, h, U):
    return h - g * (0.5 * np.trace(np.linalg.inv(g) @ h) - U)


def fd_metric_gradient(f, g, step=1e-6):
    out = np.zeros((4, 4))
    for i in range(4):
        for j in range(i, 4):
            d = np.zeros((4, 4))
            d[i, j] = d[j, i] = step
            out[i, j] = out[j, i] = (f(g + d) - f(g - d)) / (2 * step)
    return out


def elimination(rng):
    """C1 and C2 on a random Lorentzian signature h."""
    A = rng.normal(size=(4, 4))
    h, U = A.T @ ETA @ A, 0.7 + rng.random()
    g = h / U
    t_full = np.max(abs(stress(g, h, U)))
    l_err = abs(lagrangian(g, h, U) - math.sqrt(abs(np.linalg.det(h))) / U)
    grad = np.max(abs(fd_metric_gradient(lambda x: lagrangian(x, h, U), g)))
    off = np.max(abs(fd_metric_gradient(lambda x: lagrangian(x, h, U), 1.3 * g)))

    root4 = abs(np.linalg.det(h))**0.25
    gu = h / root4
    trace_free = h - 0.25 * gu * np.trace(np.linalg.inv(gu) @ h)

    def unimodular(x):
        x = x / abs(np.linalg.det(x))**0.25
        return 0.5 * np.trace(np.linalg.inv(x) @ h) - U
    grad_u = np.max(abs(fd_metric_gradient(unimodular, gu)))
    return dict(stress=t_full, lagrangian=l_err, stationary=grad, control=off,
                unit_det=abs(abs(np.linalg.det(gu)) - 1),
                trace_free=np.max(abs(trace_free)),
                value=abs(unimodular(gu) - (2 * root4 - U)), stationary_tf=grad_u)


def expansion(rng, V0):
    """C3: exact actions against their second order forms at scale V0."""
    s = gram(rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4)))
    dV = 0.4
    tr, ss = np.trace(ETA @ s), np.trace(ETA @ s @ ETA @ s)
    T = s - ETA * (0.5 * tr - dV)
    TT, Ttr = np.trace(ETA @ T @ ETA @ T), np.trace(ETA @ T)
    exact = V0 * (V0 / (V0 + dV)) * math.sqrt(abs(np.linalg.det(ETA + s / V0)))
    series = V0 + 0.5 * tr - dV - (TT - 0.5 * Ttr**2) / (4 * V0)
    cev = (V0 + 0.5 * tr - dV - (ss - 0.5 * tr**2) / (4 * V0)
           - 0.5 * dV * tr / V0 + dV**2 / V0)
    exact_tf = 2 * V0 * abs(np.linalg.det(ETA + s / V0))**0.25
    series_tf = 2 * V0 + 0.5 * tr - (ss - 0.25 * tr**2) / (4 * V0)
    return exact - series, abs(series - cev), exact_tf - series_tf


def pv_sum(direct, tail, m2, lam, pv, split=None):
    """sum_j c_j int d^4p/(2pi)^4 f(p^2, M_j^2). Tail(x), x = M^2/p^2, is f
    minus the terms the regulator conditions cancel, used above the split."""
    split = split or 3 * lam
    masses2 = [(c, m2 + s * lam**2) for c, s in pv]
    lo = quad(lambda p: p**3 * sum(c * direct(p * p, M) for c, M in masses2),
              0, split, limit=200, epsabs=0, epsrel=1e-11)[0]
    hi = quad(lambda p: p**3 * sum(c * tail(p * p, M) for c, M in masses2),
              split, np.inf, limit=200, epsabs=0, epsrel=1e-11)[0]
    return MEASURE * (lo + hi)


def log_tail(x):
    if x > 1e-2:
        return math.log1p(x) - x + x * x / 2
    return x**3 * (1 / 3 - x / 4 + x * x / 5 - x**3 / 6 + x**4 / 7)


def one_loop(m2=1.0, lam=3.0):
    """C4: tadpole, q = 0 and q^2 terms of the tensor kernel, per real field."""
    I1 = pv_sum(lambda q, M: q / (q + M),
                lambda q, M: -(M / q)**3 / (1 + M / q), m2, lam, PV4)
    I2 = pv_sum(lambda q, M: q * q / (q + M)**2,
                lambda q, M: -(4 * (M / q)**3 + 3 * (M / q)**4) / (1 + M / q)**2,
                m2, lam, PV4)
    rho = 0.5 * pv_sum(lambda q, M: math.log(q + M),
                       lambda q, M: log_tail(M / q), m2, lam, PV4)
    rho_closed = sum(c * (m2 + s * lam**2)**2 * math.log(m2 + s * lam**2)
                     for c, s in PV4) / (64 * math.pi**2)
    w0_full = I1 / 8 - I2 / 48          # <L> = 0 under the regulator
    w0_tf = I1 / 16 - I2 / 48

    # The q^2 term needs two regulator conditions only: the set of
    # carrier_sakharov.py, whose closed form a_li is the comparison.
    def x3(q, M):
        x = M / q
        return (6 * x**2 + 8 * x**3 + 3 * x**4) / (1 + x)**3 / q

    def x4(q, M):
        x = M / q
        return (10 * x**2 + 20 * x**3 + 15 * x**4 + 4 * x**5) / (1 + x)**4 / q
    J3 = pv_sum(lambda q, M: q * q / (q + M)**3, x3, m2, lam, cs.PV)
    J4 = pv_sum(lambda q, M: q**3 / (q + M)**4, x4, m2, lam, cs.PV)
    kappa = -0.25 * (-J3 / 24 + J4 / 96)
    stiffness = cs.a_li([m2], lam, "pv")           # (1/6) sum c M^2 ln M^2 / 32 pi^2
    return dict(I1=I1, I2=I2, rho=rho, rho_closed=rho_closed, w0_full=w0_full,
                w0_trace_free=w0_tf, kappa=kappa, stiffness=stiffness)


def state_metrics(rng):
    """C5: h on a random configuration, the clock and a texture."""
    random_min = min(np.linalg.eigvalsh(gram(
        rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4)))).min()
        for _ in range(50))
    mu, phi = 3 - 1 / math.sqrt(2), np.exp(2j * np.pi * np.arange(3) / 3)
    clock = np.zeros((3, 4), dtype=complex)
    clock[:, 0] = -1j * mu * phi                   # d_t phi = -i mu phi, uniform
    clock_eigs = np.linalg.eigvalsh(gram(clock))

    N, L = 24, 12.0
    psi, dx, R, _ = texture_state(N, L, amp=0.7, size=1.5)
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    K = np.meshgrid(k, k, k, indexing="ij")
    grad = np.stack([np.fft.ifftn(1j * Ki * np.fft.fftn(psi, axes=(1, 2, 3)),
                                  axes=(1, 2, 3)) for Ki in K], axis=1)   # [a, i, x, y, z]
    h = 2 * np.einsum("aixyz,ajxyz->xyzij", grad.conj(), grad).real
    eig = np.linalg.eigvalsh(h)
    core = R < 2.0
    rank3 = np.mean(eig[core][:, 0] > 1e-3 * eig[core][:, 2])
    far = eig[R > 5.5][:, 2].max() / eig[..., 2].max()
    at = (N // 2 + 2, N // 2 + 1, N // 2 + 3)
    hx = h[at]
    q_err = np.max(abs(hx - np.eye(3) * np.trace(hx) / 3
                       - 2 * neutral_tensor(grad[(slice(None), slice(None)) + at])))
    return dict(random_min=random_min, clock_eigs=clock_eigs, texture_rank3=rank3,
                texture_far=far, q_identity=q_err, texture_min=eig.min())


def contact_term():
    """C7: the rotating uniform field in the full form, and both contact energies."""
    x = np.geomspace(1e-6, 1e6, 400)
    worst = max(float(np.max(np.sqrt(1 + 2 * w * x) / (1 + x) - 1))
                for w in (0.05, 0.5, 0.9, 0.999))
    V0, rho = 1e3, 1.0                               # dust: T_00 = s_00 = rho
    s = np.zeros((4, 4))
    s[0, 0] = rho
    full = V0 * math.sqrt(abs(np.linalg.det(ETA + s / V0))) - V0 - 0.5 * rho
    free = 2 * V0 * abs(np.linalg.det(ETA + s / V0))**0.25 - 2 * V0 - 0.5 * rho
    return worst, -full * 8 * V0 / rho**2, -free * 16 * V0 / (3 * rho**2)


def ledger():
    """C6: the counted ledger's rows and the cutoff free ratio."""
    name = "Sigma_vent(UV)"
    if name not in root.WEB.state:                 # standalone: solve the web once
        with contextlib.redirect_stdout(io.StringIO()):
            R = root.derive()
            m = masses.derive(R)
            gravity.derive(R, m, couplings.derive(R, m))
    sigma = root.WEB.state[name]
    matter = sum(t.value(root.WEB.state) for t in root.WEB[name].terms
                 if not t.path[0].startswith("bridge"))
    d10, d11, n7 = root.d10, root.d11, root.n7
    bosons = (n7 + root.N_bridge + 2 * d11 * (d11 + 1)
              + d10 * (2 * n7 - (d11**2 - 1)))
    fermions = (d10 * d11 * (d11 * d10**2 + d10 + 1)       # rows of gravity.py
                + d10 * (root.n26 - 3 * (d11**2 - 1)))     # with the two right handed neutrinos
    str1 = bosons - fermions
    return dict(str1=str1, sigma=sigma, sigma_matter=matter,
                lambda_G=-(math.pi / 2) * str1 / sigma**2)


def run(report=print):
    rng = np.random.default_rng(11)
    el = elimination(rng)
    assert max(el["stress"], el["lagrangian"], el["unit_det"], el["trace_free"],
               el["value"]) < 1e-10
    assert el["stationary"] < 1e-6 and el["stationary_tf"] < 1e-6 < el["control"]

    r1, c1, t1 = expansion(np.random.default_rng(5), 2000.0)
    r2, c2, t2 = expansion(np.random.default_rng(5), 4000.0)
    assert max(c1, c2) < 1e-9
    assert abs(r1 / r2 - 4) < 0.05 and abs(t1 / t2 - 4) < 0.05

    lp = one_loop()
    assert abs(lp["I2"] / lp["I1"] - 3) < 1e-6
    assert abs(lp["I1"] / (-4 * lp["rho"]) - 1) < 1e-6
    assert abs(lp["rho"] / lp["rho_closed"] - 1) < 1e-6
    assert abs(lp["w0_full"] / (-lp["rho"] / 4) - 1) < 1e-6
    assert abs(lp["w0_trace_free"]) < 1e-6 * abs(lp["w0_full"])
    assert lp["kappa"] > 0 and abs(4 * lp["kappa"] / lp["stiffness"] - 1) < 1e-6

    st = state_metrics(rng)
    assert st["random_min"] > -1e-10 and st["texture_min"] > -1e-10
    assert st["clock_eigs"][2] < 1e-12 < st["clock_eigs"][3]
    assert st["texture_rank3"] > 0.9 and st["texture_far"] < 1e-2
    assert st["q_identity"] < 1e-10

    led = ledger()
    assert led["sigma"] > 0 > led["sigma_matter"] and led["str1"] > 0

    bound, contact_full, contact_free = contact_term()
    assert bound < 0 and abs(contact_full - 1) < 2e-3 and abs(contact_free - 1) < 2e-3

    report(f"C1 elimination: |T| = {el['stress']:.1e}, |L - sqrt|det h|/U| = "
           f"{el['lagrangian']:.1e}, dL/dg = {el['stationary']:.1e} at g = h/U "
           f"({el['control']:.2f} off it).")
    report(f"C2 traceless: |det g| - 1 = {el['unit_det']:.1e}, traceless stress "
           f"{el['trace_free']:.1e}, L = 2|det h|^(1/4) - U to {el['value']:.1e}.")
    report(f"C3 expansion: residual at V0 = 2000, 4000: full {r1:.2e}, {r2:.2e} "
           f"(ratio {r1 / r2:.2f}). Traceless {t1:.2e}, {t2:.2e} (ratio {t1 / t2:.2f}); "
           f"T.Pi.T form = the Eq. 3.5 form to {max(c1, c2):.1e}.")
    report(f"C4 one loop (PV, m = 1, L = 3): I2/I1 = {lp['I2'] / lp['I1']:.7f}, "
           f"I1/(-4 rho) = {lp['I1'] / (-4 * lp['rho']):.7f}, rho = {lp['rho']:.6f} "
           f"(closed form {lp['rho_closed']:.6f}).")
    report(f"   q = 0 tensor term: full form {lp['w0_full']:.6f} = -rho/4. Traceless "
           f"form {lp['w0_trace_free']:.1e}.")
    report(f"   q^2 tensor term: kappa = {lp['kappa']:.6e}, 4 kappa / [(1/6) sum c M^2 "
           f"ln M^2 / 32 pi^2] = {4 * lp['kappa'] / lp['stiffness']:.7f}.")
    report(f"C5 states: min eigenvalue of h over random configurations "
           f"{st['random_min']:.2e}. Clock eigenvalues {np.round(st['clock_eigs'], 6)}; "
           f"texture core rank three on {100 * st['texture_rank3']:.0f}% of points, "
           f"far field {st['texture_far']:.1e} of the peak. Traceless part = 2 Q_ij to "
           f"{st['q_identity']:.1e}.")
    report(f"C6 counted ledger: Sigma = {led['sigma']:.4f} (matter rows alone "
           f"{led['sigma_matter']:+.4f}), str(1) = {led['str1']}, "
           f"Lambda G = -(pi/2) str(1)/Sigma^2 = {led['lambda_G']:+.4f}.")
    report(f"C7 rotating uniform field, max of L/V0 over amplitude for Omega < m: {bound:.1e} "
           f"(never positive). Contact energy over rho^2/(8 V0): {contact_full:.5f}, "
           f"traceless over 3 rho^2/(16 c): {contact_free:.5f}.")
    report("RESULT: written without a background metric the carrier is the theory "
           "of an actionless metric, and that metric is the gradient bilinear "
           "over U. The q = 0 tensor term is the vacuum term alone: zero for every "
           "vacuum term in the traceless form, zero in the full form when the "
           "vacuum term cancels. The tensor stiffness per component is the weight "
           "a1 = 1/6 the ledger sums. Classical states supply h_00 (the clock) and "
           "spatial entries where a texture sits. The Lorentzian vacuum value of h "
           "is the regulated fluctuation term.")
    return dict(elimination=el, expansion_ratio=r1 / r2, expansion_ratio_tf=t1 / t2,
                expansion_form=max(c1, c2), loop=lp, states=st, ledger=led,
                contact=(bound, contact_full, contact_free))


if __name__ == "__main__":
    run()
