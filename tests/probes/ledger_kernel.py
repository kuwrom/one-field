"""The ledger's rows as tensor kernels: scalar, fermion, vector, and their sum.

composite_pole.py computes the kernel of the composite metric for the
carrier's scalar components by direct diagrams. This probe extends it to
every row of the venting ledger of interference/gravity.py, through the
second order heat trace, and sums the rows.

PREMISES. P1 each row couples minimally to the one metric g (the
actionless metric of composite_metric.py), with the endomorphism of its
own Laplacian: scalar U = xi R, Dirac U = R/4, vector U = Ricci with two
ghosts. P2 Euclidean D = 4, one loop, the proper time cutoff
rho = theta(t - t0) of sakharov_normalisation.py FACT 1. P3 massless
rows (the symmetric phase that the ledger counts).

THE HEAT TRACE. For Delta = -D^2 + U on a bundle with curvature Omega,
to second order in curvature,
    Tr e^{-t Delta} = (4 pi t)^-2 int sqrt g tr{ 1 + t (R/6 - U)
        + t^2 [ Ric f_Ric Ric + R f_R R + R f_RU U + U f_U U
                + Omega f_Om Omega ] },
    f(x) = int_0^1 d xi exp(-x xi (1 - xi)),     x = -t D^2,
    f_Ric = 1/(6x) + (f - 1)/x^2,
    f_R   = f/32 + (f - 1)/(8x) - 1/(48x) - (f - 1)/(8x^2),
    f_RU  = -f/4 - (f - 1)/(2x),   f_U = f/2,   f_Om = -(f - 1)/(2x)
(Barvinsky and Vilkovisky; in this form Codello and Zanusso
arXiv:1203.2034). Each form factor is checked here against an
independent computation before it is used.

STATEMENTS (the numbers below check each one).

  R1 (form factors). The heat trace form reproduces the direct kernel
     of composite_pole.py in all ten components at every momentum
     (f_Ric, f_R). It reproduces the direct kernel of a scalar with
     (xi/2) R phi^2 (f_RU, f_U). And f_Om = (1/4) int (1 - 2 xi)^2
     exp(-x xi (1 - xi)) is the charged scalar vacuum polarization,
     whose transversality is checked.
  R2 (rows). At low momentum each row is the linearized
     Einstein-Hilbert operator, F0/F2 -> -2, with stiffness per state
        scalar +1/6 - xi,  Weyl +1/12,  gauge vector -1/3,  Proca -1/6,
     the weights of the ledger.
  R3 (the ledger). The sum over the engine's rows has stiffness
     Sigma = 27.0000 (the ledger's total)
     and F0/F2 -> -2. On the Euclidean axis its spin-2 form factor is
     positive and its spin-0 form factor negative from q^2 t0 = 1e-3
     to 1e3, with no zero. So the massless pole is the only one. The
     amount is minimally coupled (boundary.xi_of). A single scalar row
     has no zero either.
  R4 (real time). For a massive component continued to q^2 = -s, the
     ratio of the kernel to its low momentum value is 1.00 at
     s t0 = 0.01 and 1.29 at s t0 = 1. Its real part changes sign
     between s t0 = 10 and 12 and then grows exponentially. The continuation holds below the cutoff
     scale. Above it the kernel is not that of a local field theory.
  R5 (the next order). The composite metric's own loop adds the
     graviton plus ghost coefficient of R, which depends on gauge and
     cutoff off shell: -8/3 on shell (Bastianelli, Bonezzi and Melis
     arXiv:2206.13287 Eq. 43), -23/3 in de Donder gauge off shell
     (Eq. 35; Solodukhin arXiv:1104.3712), -11/3 and -35/8 for the
     type I cutoffs of Codello, Percacci and Rahmede arXiv:0805.2909
     (Table I), in units where one real scalar is 1/6. Added to
     Sigma = 27.00 these give 19.3 to 24.3. The stiffness stays
     positive and moves by 10 to 30 percent, and K9 of
     composite_pole.py keeps the structure. No number in the lane
     depends on Sigma (interference/gravity.py), so the next order
     changes only a diagnostic ratio.

RESULT. Computed. The second order heat trace matches the direct kernels
to 1e-10, and the engine's rows sum to Sigma = 27.00 with no zero of
either form factor from q^2 t0 = 1e-3 to 1e3, so the massless pole is
the only one. The graviton's own loop keeps Sigma positive in all four
listed conventions.

Usage: python3 tests/probes/ledger_kernel.py     (a few seconds)
"""

import math
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import exp1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root  # noqa: E402
import composite_pole as cp  # noqa: E402

XG, WG, C16 = cp.XG, cp.WG, cp.C16
E0 = np.array([1.0, 0.0, 0.0, 0.0])
# trace content of one field: (tr 1, tr U / R, Ric-Ric extras, R-R extras, sign)
# extras are the coefficients of (f_RU, f_U, f_Om)
SCALAR = (1.0, 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0)
DIRAC = (4.0, 1.0, (0.0, 0.0, -2.0), (1.0, 0.25, 0.5), -1.0)   # tr Om^2 = -Riem^2/2
VECTOR = (4.0, 1.0, (0.0, 1.0, -4.0), (1.0, 0.0, 1.0), 1.0)    # tr Om^2 = -Riem^2


def scalar(xi):
    return (1.0, xi, (0.0, 0.0, 0.0), (xi, xi * xi, 0.0), 1.0)


def curvature(h, q):
    """Linearized R and Ricci. The quadratic forms of sqrt g and sqrt g R."""
    h = 0.5 * (h + h.T)
    tr, hq = np.trace(h), h @ q
    q2, qhq = q @ q, q @ hq
    ricci = 0.5 * (q2 * h + np.outer(q, q) * tr - np.outer(q, hq) - np.outer(hq, q))
    einstein = q2 * np.trace(h @ h) - q2 * tr * tr + 2 * qhq * tr - 2 * (hq @ hq)
    return (q2 * tr - qhq, ricci, -0.25 * np.trace(h @ h) + 0.125 * tr * tr,
            -0.25 * einstein, float(q2))


def form_factor_integrals(J, q2):
    """int dt rho t^-1 f_X(t q^2) for the five form factors."""
    X = XG * (1 - XG) * q2
    j00, j10, j20 = (float(WG @ J(a, 0, X)) for a in (0, 1, 2))
    j10_0, j20_0 = float(J(1, 0, 0.0)), float(J(2, 0, 0.0))
    d1, d2 = (j10 - j10_0) / q2, (j20 - j20_0) / q2**2
    return dict(Ric=j10_0 / (6 * q2) + d2,
                R=j00 / 32 + d1 / 8 - j10_0 / (48 * q2) - d2 / 8,
                RU=-j00 / 4 - d1 / 2, U=j00 / 2, Om=-d1 / 2,
                eh=j10_0, vac=j20_0)


def heat_kernel(h, q, J, content, vacuum=True):
    """Quadratic action of one field from the second order heat trace."""
    tr1, u, a_x, b_x, sign = content
    R1, ricci, vol2, eh2, q2 = curvature(h, q)
    I = form_factor_integrals(J, q2)
    a = tr1 * I["Ric"] + a_x[0] * I["RU"] + a_x[1] * I["U"] + a_x[2] * I["Om"]
    b = tr1 * I["R"] + b_x[0] * I["RU"] + b_x[1] * I["U"] + b_x[2] * I["Om"]
    out = (tr1 / 6 - u) * eh2 * I["eh"] + a * np.sum(ricci * ricci) + b * R1 * R1
    if vacuum:
        out = out + tr1 * vol2 * I["vac"]
    return -0.5 * sign * C16 * out


def direct_nonminimal(h, q, J, xi):
    """Direct diagrams for a scalar with (xi/2) R phi^2: h:V -> h:V - xi R1."""
    h = 0.5 * (h + h.T)
    trh = np.trace(h)
    A = h - 0.5 * trh * np.eye(4)
    Aq = A @ q
    q2, qAq, b0 = q @ q, q @ Aq, Aq @ Aq
    trA, trA2 = np.trace(A), np.trace(A @ A)
    X = XG * (1 - XG) * q2
    R1, _, _, eh2, _ = curvature(h, q)
    c0 = -XG * (1 - XG) * qAq - xi * R1
    b2 = (1 - 2 * XG)**2 * b0
    integrand = (6 * (trA * trA + 2 * trA2) / 24 * J(2, 0, X)
                 + 2 * (b2 / 4 + c0 * trA / 2) * J(1, 0, X)
                 - 0.5 * trh * trA * J(1, 1, X) + c0 * c0 * J(0, 0, X)
                 - c0 * trh * J(0, 1, X) + 0.25 * trh * trh * J(0, 2, X))
    I1 = 2 * C16 * float(J(2, 0, 0.0))
    return ((I1 / 4) * (0.5 * np.trace(h @ h) - 0.25 * trh * trh)
            - 0.25 * C16 * float(WG @ integrand)
            + 0.5 * xi * eh2 * C16 * float(J(1, 0, 0.0)))


def row(kind, xi=0.0):
    """Kernel of one physical field, vacuum term removed."""
    def k(h, q, J):
        if kind == "scalar":
            return heat_kernel(h, q, J, scalar(xi), vacuum=False)
        if kind == "dirac":
            return heat_kernel(h, q, J, DIRAC, vacuum=False)
        ghosts = 2 if kind == "vector" else 1          # gauge vector or Proca
        return (heat_kernel(h, q, J, VECTOR, vacuum=False)
                - ghosts * heat_kernel(h, q, J, SCALAR, vacuum=False))
    return k


def form_factors(kernel, q, J):
    theta = np.eye(4) - np.outer(q, q) / (q @ q)
    w, v = np.linalg.eigh(theta)
    tt = (np.outer(v[:, -1], v[:, -2]) + np.outer(v[:, -2], v[:, -1])) / math.sqrt(2)
    return kernel(tt, q, J) / cp.spin2(tt, theta), kernel(theta, q, J) / cp.spin0(theta, theta)


def ledger_rows():
    """The engine's rows as (field kind, number of fields, xi): gravity.py counts."""
    d10, d11, n7 = root.d10, root.d11, root.n7
    states = d10 * d11 * (d11 * d10**2 + d10 + 1)                  # 90 fermion states
    return [("vector", d11 * (d11 + 1), 0.0),                      # 24 gauge states
            ("dirac", states / 4, 0.0),
            ("vector", 2 * n7 - (d11**2 - 1), 0.0),                # 12 coset states
            ("scalar", n7, 0.0),
            ("scalar", root.N_bridge, 0.0)]                         # the amount, minimally coupled


def timelike(s, t0, m2=1.0, eps=1e-10):
    """R4: 4 F2/q^2 over its q -> 0 value at q^2 = -s, one massive component."""
    u, w = np.polynomial.legendre.leggauss(500)
    u, w = 0.5 * (u + 1), 0.5 * w
    x, w = 0.5 * (1 - np.cos(math.pi * u)), w * 0.5 * math.pi * np.sin(math.pi * u)

    def E(k, z):
        e = exp1(z)
        for j in range(1, k):
            e = (np.exp(-z) - z * e) / j
        return e
    D = m2 - x * (1 - x) * s - 1j * eps
    F2 = -0.125 * C16 * np.sum(w * (E(3, t0 * D) - E(3, t0 * m2 + 0j))) / t0**2
    stiffness = E(2, t0 * m2 + 0j).real / (192 * math.pi**2 * t0)
    return 4 * F2 / (-s) / stiffness


def run(report=print):
    rng = np.random.default_rng(1)
    # R1: the form factors against independent computations.
    check = {}
    for name, (J, _, _) in {"pv": cp.pauli_villars(1.0, 3.0),
                            "pt": cp.proper_time(1.0, 1e-2)}.items():
        minimal = nonminimal = 0.0
        for size in (0.1, 0.7, 2.0, 6.0):
            q = rng.normal(size=4)
            q *= size / np.linalg.norm(q)
            for _ in range(3):
                h = rng.normal(size=(4, 4))
                h = h + h.T
                minimal = max(minimal, abs(heat_kernel(h, q, J, SCALAR) / cp.gamma2(h, q, J) - 1))
                for xi in (0.1, 1 / 6, -0.3):
                    nonminimal = max(nonminimal, abs(heat_kernel(h, q, J, scalar(xi))
                                                     / direct_nonminimal(h, q, J, xi) - 1))
        check[name] = (minimal, nonminimal)
        assert minimal < 1e-10 and nonminimal < 1e-10
    Jm, _, _ = cp.proper_time(1.0, 1e-2)
    X = XG * (1 - XG) * 4.0
    lhs = 2 * (float(WG @ Jm(1, 0, X)) - float(Jm(1, 0, 0.0)))
    polarization = float(WG @ ((1 - 2 * XG)**2 * Jm(0, 0, X)))
    omega_err = max(abs(lhs + 4.0 * polarization),
                    abs(polarization - 4 * form_factor_integrals(Jm, 4.0)["Om"]))
    assert omega_err < 1e-10

    # R2: the rows, massless, proper time cutoff t0 = 1.
    J, _, _ = cp.proper_time(1e-10, 1.0)
    unit = 1.0 / (32 * math.pi**2)
    q = 3e-3 * E0
    xi_b = float(root.xi_bridge)
    weights = {}
    for label, kind, xi, states in (("scalar", "scalar", 0.0, 1), ("bridge scalar", "scalar", xi_b, 1),
                                    ("Dirac", "dirac", 0.0, 4), ("gauge vector", "vector", 0.0, 2),
                                    ("Proca", "proca", 0.0, 3)):
        F2, F0 = form_factors(row(kind, xi), q, J)
        weights[label] = (4 * F2 / (q @ q) / unit / states, F0 / F2)
    want = {"scalar": 1 / 6, "bridge scalar": 1 / 6 - xi_b, "Dirac": 1 / 12,
            "gauge vector": -1 / 3, "Proca": -1 / 6}
    for label, (wgt, ratio) in weights.items():
        assert abs(wgt - want[label]) < 2e-4 and abs(ratio + 2) < 5e-3

    # R3: the ledger.
    rows = ledger_rows()

    def ledger(h, q, J):
        return sum(n * row(kind, xi)(h, q, J) for kind, n, xi in rows)
    F2, F0 = form_factors(ledger, q, J)
    sigma_kernel = 4 * F2 / (q @ q) / unit
    sigma_rows = sum(n * want[{"scalar": "scalar" if xi == 0 else "bridge scalar",
                               "dirac": "Dirac", "vector": "gauge vector"}[kind]]
                     * {"scalar": 1, "dirac": 4, "vector": 2}[kind] for kind, n, xi in rows)
    assert abs(sigma_rows - 27.0) < 1e-3 and abs(sigma_kernel - sigma_rows) < 2e-2
    assert abs(F0 / F2 + 2) < 1e-2

    def channel(i):
        return lambda ly: form_factors(ledger, math.sqrt(math.exp(ly)) * E0, J)[i]
    grid = np.linspace(math.log(1e-3), math.log(1e3), 120)
    zeros, signs = [], []
    for i in (0, 1):
        v = [channel(i)(x) for x in grid]
        signs.append(v)
        zeros.append([math.exp(brentq(channel(i), a, b))
                      for a, b, va, vb in zip(grid[:-1], grid[1:], v[:-1], v[1:]) if va * vb < 0])
    assert zeros == [[], []] and all(f2 > 0 > f0 for f2, f0 in zip(*signs))
    one_scalar = [form_factors(row("scalar"), math.sqrt(y) * E0, J) for y in np.geomspace(1e-3, 1e3, 40)]
    assert all(f2 > 0 > f0 for f2, f0 in one_scalar)

    # R4: real time, one massive component.
    real_time = {st: timelike(st / 1e-2, 1e-2) for st in (0.01, 1.0, 10.0, 12.0, 30.0)}
    assert abs(real_time[0.01].real - 1) < 1e-2 and abs(real_time[1.0].real - 1.285) < 5e-3
    assert real_time[10.0].real > 0 > real_time[12.0].real and real_time[30.0].real < -20
    growth = math.log(abs(timelike(4000.0, 1e-2).real / real_time[30.0].real)) / 10.0

    # R5: the next order, the graviton's own coefficient of R.
    graviton = {"on shell": -8 / 3, "de Donder off shell": -23 / 3,
                "type Ia cutoff": -11 / 3, "type Ib cutoff": -35 / 8}
    next_order = {k: sigma_rows + v for k, v in graviton.items()}
    assert min(next_order.values()) > 0

    report(f"R1 heat trace against direct diagrams: minimal scalar {check['pv'][0]:.1e} (PV), "
           f"{check['pt'][0]:.1e} (proper time). With xi R phi^2 {check['pv'][1]:.1e}, "
           f"{check['pt'][1]:.1e}. Charged scalar polarization {omega_err:.1e}.")
    for label, (wgt, ratio) in weights.items():
        report(f"R2 {label:13s}: stiffness per state {wgt:+.5f} (ledger {want[label]:+.5f}), "
               f"F0/F2 = {ratio:+.4f}.")
    report(f"R3 ledger: Sigma from the kernel {sigma_kernel:.4f}, from the rows {sigma_rows:.4f}; "
           f"F0/F2 = {F0 / F2:+.4f}.")
    report("   Euclidean axis: F2 > 0 > F0 from q^2 t0 = 1e-3 to 1e3, no zero; a single scalar row has none.")
    report("R4 real time, 4F2/q^2 over its low momentum value at s t0 = "
           + ", ".join(f"{st}: {v.real:+.3f}{v.imag:+.3f}i" for st, v in real_time.items())
           + f". The real part changes sign between s t0 = 10 and 12. Growth rate "
           f"{growth:.3f} per unit s t0 from 30 to 40.")
    report("R5 next order, Sigma + graviton: "
           + ", ".join(f"{k} {v:.2f}" for k, v in next_order.items())
           + ". Positive in every convention, the structure fixed by locality.")
    report("RESULT: every row of the ledger is the linearized Einstein-Hilbert operator at "
           "low momentum with the ledger's weight, and their sum has stiffness Sigma. "
           "On the Euclidean side at every momentum computed, and below the cutoff scale in real "
           "time, the massless pole is the only one.")
    return dict(check=check, omega_error=omega_err, weights=weights, sigma_kernel=sigma_kernel,
                sigma_rows=sigma_rows, ledger_ratio=F0 / F2, zeros_spin2=zeros[0],
                zeros_spin0=zeros[1], real_time=real_time, growth=growth,
                next_order=next_order)


if __name__ == "__main__":
    run()
