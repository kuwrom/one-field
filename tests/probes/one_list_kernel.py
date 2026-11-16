"""One action on one list: every boundary row as one kernel of the composite metric.

THE QUESTION. boundary_content.py builds one list of boundary rows. In
it the bridge carries its searchlight label, and both the G2 running
and the stiffness see that label. The list sums per state weights: a
row of mass m enters with E2(m^2 t0) (its P3). The stiffness of this
labelled list is -0.175 with the regulator at the boundary and 3.698 at
the self-consistent regulator. ledger_kernel.py computes the counted
rows (the rows the ledger counts) as tensor kernels of one composite
metric. Written as one metric free action that contains every row of a
list, the induced kernel is a single object. Does it carry the weights
the sum assumes? Up to what momentum does it keep the Einstein sign?
The engine's one list reads the bridge as amount: light, minimally
coupled and with no share of the running. Its kernel is
tests/probes/horizon_amount.py A4.

PREMISES. P1 the metric free form of each row. Scalar rows: the
composite metric g = h/U of the gradient bilinear, L = sqrt|det h|/U
(composite_metric.py C1; Carone, Erlich and Vaman arXiv:1610.08521), the
bridge with xi R phi^2 in the curvature of that metric. Vector rows: the
Hermitian determinant sqrt|det(h + i sqrt(U) F)|/U of
determinant_vectors.py, Maxwell on the composite metric at quadratic
order (V2, V5), so the heat trace vector row applies. Fermion rows: the
composite vierbein (Wetterich arXiv:hep-th/0307145), whose pole carries
half the scalar weight per state (Carone, Erlich and Vaman
arXiv:1812.08201: N Dirac flavors give M_P^2 = N m^2/(24 pi^2 |eps|),
Eq. 3.32, N scalars N m^2/(48 pi^2 |eps|), Eq. 2.13 at D = 4 - eps).
Heavy rows: the same forms with a mass term. The one action is the sum
of these over a list, and its kernel is the sum of the rows'
kernels. P2 the second order heat trace of ledger_kernel.py (R1 checks
it against direct diagrams), Euclidean D = 4, one loop of the rows,
under the proper time cutoff rho = theta(t - t0) of
sakharov_normalisation.py FACT 1. Vector rows carry the Ricci
endomorphism and two ghosts (contact c = +1). A mass enters as a weight
on the proper time, Tr exp(-t(Delta + m^2)) = exp(-m^2 t)
Tr exp(-t Delta): a heavy row is its massless row read with
rho(t) = theta(t - t0) exp(-m^2 t), the weight of
composite_pole.proper_time(m2, t0), whose q^2 term carries
int_t0 dt t^-2 exp(-m^2 t) = E2(m^2 t0)/t0. P3 the lists.
The counted list is ledger_kernel.py's rows, all massless. Field for
field it is the counted list of boundary_content.py. The labelled list is
boundary_content.rows('boundary'): 36 vector states and 90 fermion
states massless, the 182 bridge channels at mass M_b with
xi = 1/(48 pi), no G2 scalars. t0 is in units of M_b^-2, and the
boundary identification M_b^2 t0 = Sigma/(2 pi) (root.py,
boundary_content.py P1) places the regulator. The carrier's own six
components (at most +1, composite_pole.py K5) and the composite metric's
own loop (-8/3, ledger_kernel.py R5, beyond one loop of the rows) are
left out of the labelled list, as in boundary_content.py B5. Together
they lower the stiffness.

STATEMENTS (the numbers below check each one).

  K1 (heavy rows). A scalar row of mass m with xi R phi^2 is, at low
     momentum, the linearized Einstein-Hilbert operator, F0/F2 = -2 to
     4.3e-6, with stiffness (1/6 - xi) E2(m^2 t0) per real component:
     the kernel over that weight is 1 to 1.0e-6 at m^2 t0 = 0.01, 0.3,
     1, 3 (E2 = 0.9497, 0.4691, 0.1485, 0.0106) for xi = 0, 1/(48 pi),
     0.1, -0.3, and the heavy row's heat trace matches the direct
     diagrams of ledger_kernel.py to 2.3e-13. This is boundary_content.py
     P3, read off the kernel.
  K2 (the counted list). The summed kernel has stiffness 27.0000 and
     F0/F2 = -2.0000, and its spin-2 form factor keeps the Einstein
     sign from q^2 t0 = 1e-3 to 1e3, with no zero. The amount is
     minimally coupled (boundary.xi_of). The zero of ledger_kernel.py R3
     belongs to a row at xi = 1/(48 pi). With
     M_b^2 t0 = Sigma/(2 pi) = 4.297 the regulator scale t0^(-1/2) is
     0.482 M_b, inside the boundary. Barcelo, Liberati and Visser
     (arXiv:gr-qc/0505065v4 Sec. 7.9, Sec. 7.8 of the Living Reviews
     text) take the scale where the geometric picture disappears to be
     the cutoff itself. Read that way, the counted list breaks down at
     its regulator, 0.482 M_b. It departs from the Einstein value by 10
     and 50 percent at q^2 t0 = 0.072 and 0.805, 0.129 and 0.433 M_b,
     inside the boundary.
  K3 (the labelled list). Its kernel has the stiffness of
     boundary_content.sigma, -0.1749 at t0 = 1 and +3.6978 at the
     self-consistent t0 = 0.5885, both to 1.1e-5, with F0/F2 = -2 to
     1.0e-4. At t0 = 1 the stiffness is negative and F2 stays negative
     up to q^2 t0 = 1e3, with no zero. The massless spin-2 pole between
     conserved sources has negative residue, so it is a ghost. The
     static potential between two positive masses is the attractive
     Newton value times -1.000, so it is repulsive. Its spin-0 form factor changes sign
     at q^2 t0 = 0.043, 0.208 M_b. At the self-consistent t0 it is
     attractive (+1.000), and the spin-2 form factor vanishes at
     q^2 t0 = 0.2295, q = 0.624 M_b. It keeps the Einstein sign only
     below 0.624 M_b, under the boundary and under the regulator scale
     1.304 M_b. It departs from the Einstein value by 10 and 50 percent
     at 0.137 and 0.380 M_b. So on the labelled list the composite metric
     loses the Einstein sign inside the very boundary that defines the
     list.
  K4 (the carrier alone). The six real components have F2 > 0 from
     q^2 t0 = 1e-3 to 1e4, no zero, falling monotonically below the
     Einstein value: by 10 percent at q^2 t0 = 0.259 and 50 percent at
     3.317 massless, 0.315 and 3.522 at m^2 t0 = 0.01, 0.455 and 4.214
     at 0.1. The scale rises with the mass and tends to the massless
     value as m^2 t0 -> 0 (shifts of 0.7 and 0.2 percent at
     m^2 t0 = 1e-4, 21.4 and 6.2 percent at 0.01): the mass enters only
     through m^2 t0. In the limit m^2 t0 -> 0 the departure sits at
     q = 0.509 and 1.821 t0^(-1/2), and a heavy carrier departs at
     q^2 t0 of order one too (10 percent at 0.737, 0.994 and 1.068 for
     m^2 t0 = 1, 10 and 100), below its mass (0.86, 0.32 and 0.10 of
     it). At every mass the carrier's breakdown computation returns the
     regulator scale and nothing beyond it. So taking the cutoff to be
     that scale fixes no t0.

RESULT. Written as one metric free action, a list induces one kernel.
That kernel carries the weights boundary_content.py sums: at low
momentum a heavy row is the Einstein operator with stiffness
(1/6 - xi) E2(m^2 t0). Under the ledger's sharp cutoff the labelled
list has two cases. With the regulator at the boundary its stiffness
is -0.175, the massless spin-2 pole is a ghost and the static potential
repels. At the self-consistent regulator its stiffness is 3.698 and it
attracts, but its spin-2 form factor vanishes at 0.624 M_b, below the
boundary and the regulator. The counted list keeps the Einstein sign
at every momentum computed. Both lists depart from the Einstein value
inside the boundary. The carrier alone keeps the Einstein sign over the
whole range and departs from the Einstein value at q^2 t0 of order one
at every mass. The breakdown computation returns t0 and no scale beyond
it.

Usage: python3 tests/probes/one_list_kernel.py     (about a second)
"""

import math
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B, root as r  # noqa: E402
import boundary_content as bc  # noqa: E402
import composite_pole as cp  # noqa: E402
import ledger_kernel as lk  # noqa: E402

XI = float(r.xi_bridge)
UNIT = 1.0 / (32 * math.pi**2)           # 4 F2/q^2 = Sigma UNIT/t0 (ledger_kernel.py R2)
Y_LOW = 1e-7                             # the low momentum point, q^2 t0
E0 = np.array([1.0, 0.0, 0.0, 0.0])
HEAVY = (0.01, 0.3, 1.0, 3.0)            # m^2 t0 of K1
XIS = (0.0, XI, 0.1, -0.3)               # minimal, the bridge's, two more
CARRIER = (0.0, 1e-4, 0.01, 0.1)         # m^2 t0 of K4
CARRIER_HEAVY = (1.0, 10.0, 100.0)       # heavy carriers, K4
LEVELS = (0.1, 0.5)                      # departures from the Einstein value, K4
FIELD = {"vector": ("vector", 2), "weyl": ("dirac", 4), "scalar": ("scalar", 1), "bridge": ("scalar", 1)}


def kernel_rows(table):
    """boundary.py rows as (field, fields, xi, m^2/M_b^2): light and amount massless, boundary at M_b.

    The bridge carries the label's xi only when it carries the label (boundary.xi_of). The amount is minimally coupled.
    """
    out = []
    for _, states, kind, _, reading in table:
        if reading in ("light", "amount", "boundary"):
            field, per = FIELD[kind]
            out.append((field, states / per, B.xi_of(reading, XI) if kind == "bridge" else 0.0,
                        1.0 if reading == "boundary" else 0.0))
    return out


def content(table):
    """Fields per (field, xi): what a list counts, however its rows are split."""
    out = {}
    for field, n, xi, _ in table:
        out[field, round(xi, 12)] = out.get((field, round(xi, 12)), 0) + n
    return out


def kernel(table, t0):
    """The one action's kernel: every row's heat trace kernel under its own weight exp(-m^2 t), summed."""
    parts = [(n, lk.row(field, xi), cp.proper_time(m2, t0)[0]) for field, n, xi, m2 in table]
    return lambda h, q, _J=None: sum(n * k(h, q, J) for n, k, J in parts)


def form_factors(table, t0, y):
    """F2, F0 of the summed kernel at Euclidean q^2 t0 = y."""
    return lk.form_factors(kernel(table, t0), math.sqrt(y / t0) * E0, None)


def stiffness(table, t0):
    """Sigma, from 4 F2/q^2 = Sigma UNIT/t0 at q^2 t0 = Y_LOW, and F0/F2 there."""
    F2, F0 = form_factors(table, t0, Y_LOW)
    return 4 * F2 * t0 * t0 / (Y_LOW * UNIT), F0 / F2


def static_sign(table, t0):
    """Static potential between two conserved sources over the attractive Newton value at |G|."""
    q = math.sqrt(Y_LOW / t0) * E0
    F2, F0 = lk.form_factors(kernel(table, t0), q, None)
    theta = np.eye(4) - np.outer(q, q) / (q @ q)
    T = np.diag([0.0, 0.0, 0.0, 1.0])                  # a static source, u.q = 0 (composite_pole.py K5)
    potential = -0.125 * (cp.spin2(T, theta) / F2 + cp.spin0(T, theta) / F0)
    G = 1.0 / (16 * math.pi * abs(4 * F2 / (q @ q)))
    return potential / (-4 * math.pi * G / (q @ q))


def axis(table, t0, lo=1e-3, hi=1e3, n=120):
    """F2 on a log grid of q^2 t0, and its first zero (None if it has none).

    The Feynman parameter nodes hold F2 to 1e-3 up to q^2 t0 = 1e3 and
    lose the endpoints beyond; kernel_sign.py has the closed form.
    """
    K = kernel(table, t0)

    def f(ly):
        return lk.form_factors(K, math.sqrt(math.exp(ly) / t0) * E0, None)[0]
    grid = np.linspace(math.log(lo), math.log(hi), n)
    v = [f(x) for x in grid]
    for a, b, va, vb in zip(grid[:-1], grid[1:], v[:-1], v[1:]):
        if va * vb < 0:
            return v, math.exp(brentq(f, a, b))
    return v, None


def departures(table, t0):
    """q^2 t0 where the summed F2 falls 10 and 50 percent below its Einstein value Sigma UNIT y/(4 t0^2)."""
    s = stiffness(table, t0)[0]
    f = lambda y: form_factors(table, t0, y)[0] / (s * UNIT * y / (4 * t0 * t0))
    return [brentq(lambda y: f(y) - (1 - lv), 1e-5, 1e2) for lv in LEVELS]


def first_zero(fn, lo=1e-4, hi=1e3, n=200):
    grid = np.linspace(math.log(lo), math.log(hi), n)
    v = [fn(x) for x in grid]
    for a, b, va, vb in zip(grid[:-1], grid[1:], v[:-1], v[1:]):
        if va * vb < 0:
            return math.exp(brentq(fn, a, b))
    return None


def einstein_ratio(m2t0, y):
    """The carrier's F2 over its Einstein value (N/6) E2(m^2 t0) UNIT q^2/4, at t0 = 1."""
    table = [("scalar", cp.N_COMPONENTS, 0.0, m2t0)]
    return form_factors(table, 1.0, y)[0] / (cp.N_COMPONENTS / 6 * float(expn(2, m2t0)) * UNIT * y / 4)


def run(report=print):
    out = {}

    # K1: a heavy scalar row against the weight (1/6 - xi) E2(m^2 t0) of boundary_content.py P3.
    rng = np.random.default_rng(2)
    heavy, direct = {}, 0.0
    for z in HEAVY:
        J = cp.proper_time(z, 1.0)[0]
        ratio, shape = {}, {}
        for xi in XIS:
            s, shape[xi] = stiffness([("scalar", 1, xi, z)], 1.0)
            ratio[xi] = s / ((1 / 6 - xi) * float(expn(2, z)))
            for size in (0.1, 0.7, 2.0):
                q = rng.normal(size=4)
                q *= size / np.linalg.norm(q)
                h = rng.normal(size=(4, 4))
                h = h + h.T
                direct = max(direct, abs(lk.heat_kernel(h, q, J, lk.scalar(xi))
                                         / lk.direct_nonminimal(h, q, J, xi) - 1))
        heavy[z] = dict(E2=float(expn(2, z)), bridge=ratio[XI], worst=max(abs(v - 1) for v in ratio.values()),
                        shape=max(abs(v + 2) for v in shape.values()))
    assert direct < 1e-10
    assert all(v["worst"] < 1e-5 and v["shape"] < 1e-4 for v in heavy.values())
    out.update(heavy=heavy, direct=direct)
    for z, v in heavy.items():
        report(f"K1 m^2 t0 = {z:4}: E2 = {v['E2']:.4f}. Kernel over (1/6 - xi) E2 = {v['bridge']:.7f} "
               f"(bridge xi), worst |ratio - 1| = {v['worst']:.1e}, |F0/F2 + 2| <= {v['shape']:.1e} "
               "over xi = 0, 1/(48 pi), 0.1, -0.3")
    report(f"   the heavy row's heat trace against direct diagrams: {direct:.1e}")

    # K2: the counted list. The boundary identification M_b^2 t0 = Sigma/(2 pi).
    booked = [(field, n, xi, 0.0) for field, n, xi in lk.ledger_rows()]
    assert content(booked) == content(kernel_rows(bc.rows(bridge="amount", scalars="light")))
    sig_b, shape_b = stiffness(booked, 1.0)
    values_b, zero_b = axis(booked, 1.0)
    t0_b = sig_b / (2 * math.pi)
    reg_b = t0_b**-0.5
    assert abs(sig_b - 27.0) < 1e-3 and abs(shape_b + 2) < 1e-3
    dep_b = departures(booked, 1.0)
    qdep_b = [math.sqrt(y / t0_b) for y in dep_b]
    assert zero_b is None and all(v > 0 for v in values_b) and reg_b < 1.0 and max(qdep_b) < 1.0
    out.update(sigma_booked=sig_b, shape_booked=shape_b, zero_booked=zero_b, t0_booked=t0_b,
               regulator_booked=reg_b, departure_booked=dep_b, q_departure_booked=qdep_b)
    report(f"K2 counted list: Sigma = {sig_b:.4f} from the kernel, F0/F2 = {shape_b:+.4f}; F2 > 0 from "
           f"q^2 t0 = 1e-3 to 1e3, no zero")
    report(f"   M_b^2 t0 = Sigma/(2 pi) = {t0_b:.3f}: regulator t0^-1/2 = {reg_b:.3f} M_b < M_b; "
           f"10%/50% below Einstein at q^2 t0 = {dep_b[0]:.3f}/{dep_b[1]:.3f}, "
           f"q = {qdep_b[0]:.3f}/{qdep_b[1]:.3f} M_b")

    # K3: the labelled list, with the regulator at the boundary and at the self-consistent t0.
    table = bc.rows(bridge="boundary")
    one = kernel_rows(table)
    counts = content(one)
    assert counts == {("vector", 0.0): 18, ("dirac", 0.0): 22.5, ("scalar", round(XI, 12)): 182}
    t0_s = bc.fixed_points(table, 1.0)[0] / (2 * math.pi)
    k3 = {}
    for name, t0 in (("boundary", 1.0), ("self", t0_s)):
        s, shape = stiffness(one, t0)
        values, zero = axis(one, t0)
        k3[name] = dict(t0=t0, sigma=s, want=bc.sigma(table, t0, 1.0), shape=shape,
                        sign=static_sign(one, t0), zero=zero, negative=all(v < 0 for v in values))
    at_b, at_s = k3["boundary"], k3["self"]
    agree = max(abs(v["sigma"] - v["want"]) for v in k3.values())
    shape_one = max(abs(v["shape"] + 2) for v in k3.values())
    assert agree < 1e-3 and shape_one < 1e-3
    assert at_b["sigma"] < 0 and at_b["negative"] and at_b["zero"] is None and abs(at_b["sign"] + 1) < 1e-3
    assert at_s["sigma"] > 0 and abs(at_s["sign"] - 1) < 1e-3
    qz_s = math.sqrt(at_s["zero"] / t0_s)
    f0_zero = first_zero(lambda ly: form_factors(one, 1.0, math.exp(ly))[1])
    dep_s = departures(one, t0_s)
    qdep_s = [math.sqrt(y / t0_s) for y in dep_s]
    assert qz_s < 1.0 < t0_s**-0.5 and f0_zero is not None and max(qdep_s) < 1.0
    out.update(labelled_list=k3, agree=agree, shape_labelled=shape_one, t0_self=t0_s, regulator_self=t0_s**-0.5,
               zero_labelled=at_s["zero"], q_zero_labelled=qz_s, f0_zero_boundary=f0_zero,
               departure_labelled=dep_s, q_departure_labelled=qdep_s)
    report(f"K3 labelled list ({2 * counts['vector', 0.0]:.0f} vector and {4 * counts['dirac', 0.0]:.0f} fermion states "
           f"massless, {counts['scalar', round(XI, 12)]:.0f} bridge channels at M_b, no G2 scalars)")
    report(f"   t0 = 1 M_b^-2: Sigma = {at_b['sigma']:+.4f} (boundary_content {at_b['want']:+.4f}); "
           f"F2 < 0 from q^2 t0 = 1e-3 to 1e3: {at_b['negative']}. Static potential over attractive "
           f"Newton {at_b['sign']:+.3f}; F0 changes sign at q^2 t0 = {f0_zero:.4f}, q = {math.sqrt(f0_zero):.3f} M_b")
    report(f"   self-consistent t0 = {t0_s:.4f} M_b^-2: Sigma = {at_s['sigma']:+.4f} (boundary_content "
           f"{at_s['want']:+.4f}), static potential {at_s['sign']:+.3f}; F2 vanishes at q^2 t0 = "
           f"{at_s['zero']:.4f}, q = {qz_s:.3f} M_b. Regulator t0^-1/2 = {t0_s**-0.5:.3f} M_b; 10%/50% below "
           f"Einstein at q = {qdep_s[0]:.3f}/{qdep_s[1]:.3f} M_b")
    report(f"   kernel against boundary_content.sigma to {agree:.1e}; |F0/F2 + 2| <= {shape_one:.1e}")

    # K4: the carrier's six components alone. Where F2 leaves its Einstein value.
    grid = np.geomspace(1e-3, 1e4, 60)
    k4 = {}
    for z in CARRIER:
        ratios = [einstein_ratio(z, y) for y in grid]
        assert all(v > 0 for v in ratios) and all(a > b for a, b in zip(ratios[:-1], ratios[1:]))
        k4[z] = [brentq(lambda y: einstein_ratio(z, y) - (1 - lv), 1e-4, 1e4) for lv in LEVELS]
    shift = {z: [d / k4[0.0][i] - 1 for i, d in enumerate(k4[z])] for z in CARRIER}
    assert all(0 <= v < 1e-2 for v in shift[1e-4])
    assert all(k4[a][i] < k4[b][i] for a, b in zip(CARRIER[:-1], CARRIER[1:]) for i in (0, 1))
    heavy_dep = {z: [brentq(lambda y: einstein_ratio(z, y) - (1 - lv), 1e-4, 1e4) for lv in LEVELS]
                 for z in CARRIER_HEAVY}
    assert all(0.5 < v[0] < 1.5 for v in heavy_dep.values())
    out.update(departure=k4, shift=shift, departure_heavy=heavy_dep)
    for z in CARRIER:
        report(f"K4 carrier, m^2 t0 = {z:6}: F2 > 0, falling below Einstein monotonically from q^2 t0 = 1e-3 "
               f"to 1e4; 10% at q^2 t0 = {k4[z][0]:.3f}, 50% at {k4[z][1]:.3f} (shift from massless "
               f"{100 * shift[z][0]:.1f}%, {100 * shift[z][1]:.1f}%)")
    report(f"   massless: 10% at q = {math.sqrt(k4[0.0][0]):.3f} t0^-1/2, 50% at {math.sqrt(k4[0.0][1]):.3f} "
           "t0^-1/2: the regulator scale. Heavy, 10% at q^2 t0 = "
           + ", ".join(f"{v[0]:.3f} (m^2 t0 = {z:g})" for z, v in heavy_dep.items()))
    report("RESULT: the one action's kernel carries the weights boundary_content.py sums. Under the ledger's sharp "
           "cutoff the labelled list's spin-2 pole is a ghost with the regulator at the boundary, and at its "
           f"self-consistent regulator its spin-2 form factor vanishes at {qz_s:.3f} M_b, inside the boundary "
           "(counted list: no zero). The carrier departs from Einstein at the regulator scale at every mass.")
    return out


if __name__ == "__main__":
    run()
