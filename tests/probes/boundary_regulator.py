"""The labelled bridge's weight at the boundary under the regulator shapes the sources write down.

THE QUESTION. Read with its searchlight label as running charge, the
182 bridge rows of interference/boundary.py can leave the running only
at the boundary's mass (bridge="boundary"). Gate S then passes only if
they enter the stiffness at the boundary with at least 0.370 of a
massless row's weight (0.462 with the composite metric's own loop,
boundary_content.py B5). The engine's one list reads the bridge as
amount: light, with no share of the running, and unaffected by any
regulator's shape (tests/probes/horizon_amount.py). The ledger's sharp
proper time cutoff gives E2(1) = 0.148. Over monotone regulators the
weight can lie anywhere from 0.148 to 1, so the regulator's shape
decides the sign. What do the sources fix, and what does each shape
they write down give?

PREMISES. P1 the labelled list and its weights (boundary.py): the light
rows give Sigma_light = 13.5 - 18 c, the bridge gives 182 (1/6 - xi) =
29.13 at a massless row's weight, and with G anchored,
t0 = Sigma/(2 pi) M_b^-2. P2 one normalization for every shape, with
one rho shared by every row. A row of mass m enters the Einstein term
with Int dt t^-2 rho(t) exp(-m^2 t), t0 is fixed by the massless row
(Int dt t^-2 rho = 1/t0), and the row's weight is w(x), x = m^2 t0. Its
running threshold is T(x) = Int dt t^-1 rho(t) exp(-m^2 t). A cutoff
function chi(D^2/Lambda^2) is written in one loop units as
chi(u) = Int ds s^-1 rho(s) exp(-s u), so w(x) = F(x f2)/f2 with
F(y) = Int_0^inf chi(u + y) du, f2 = F(0), and T(x) = chi(x f2). P3 the
anchor (root.py): alpha_G2 = 1/(24 pi) at M_b, run by the light rows'
b0 = 32/3, puts M_Pl at 22 ppm. P4 the composite metric's own loop,
-8/3 (ledger_kernel.py R5).

STATEMENTS (the numbers below check each one).

  R1 (what the sources write down). Each source fixes a condition or a
     count and writes one example. The sharp cutoff theta(t - t0), the
     ledger's own ("nothing particularly special", Visser
     gr-qc/0204062; the m -> inf optimum of the proper time flows,
     Mazza and Zappala hep-th/0106230): 0.148. Path integral duality
     exp(-t0/t), the one shape a stated principle fixes for every row
     and for the running (Padmanabhan hep-th/9608182 Eq. 4;
     Srinivasan, Sriramkumar and Padmanabhan gr-qc/9710104
     Eqs. 147-150): 0.280. The edge of the proper time window
     P(2, t/t0) (Litim and Pawlowski hep-th/0202188), which is two
     equal Pauli-Villars fields, a count no Pauli-Villars source uses:
     1 - ln 2 = 0.307. The flow family stays at or below it for m >= 1
     (Bonanno and Reuter hep-th/0410191). Three equal Pauli-Villars
     fields, the minimum count, with the masses set equal to simplify
     (Akama, Prog. Theor. Phys. 60 (1978) 1900; Li and Vaman
     arXiv:2512.06132 Eq. 2.41): 0.235. Four (Eq. 2.51): 0.207. Frolov
     and Fursaev's restrictions at the Demers-Lafrance-Myers masses,
     one resolution (hep-th/9802010 Eqs. 4.18-4.25): 0.249. The bosonic
     entropy (Dong, Khalkhali and van Suijlekom arXiv:1903.09624), with
     its amplitude fixed by its log: 0.181. Band pass cutoff functions
     at chi(0) = 1: eigenvalue counting theta(1 - u): 0. The heat
     kernel exp(-u): exp(-1) = 0.368 (a near miss). The Gaussian in
     D^2: 0.210. The exponential in |D| (Lizzi arXiv:1805.00411
     Eq. 2.3): 0.587. The fermionic entropy (Chamseddine, Connes and
     van Suijlekom arXiv:1809.02944): 0.426. Outside P2, Pauli-Villars
     applied per species gives 0.057.
  R2 (the datum's scheme). If one regulator serves both the running and
     the stiffness, every row adds its regulated threshold to
     1/alpha_G2. A light row adds the constant kappa =
     lim (T(x) + ln x), which is -gamma for the sharp cutoff and
     -2 gamma for duality. The bridge adds (13/(12 pi)) T(1). With the
     datum at the regulator, the light rows alone would move M_Pl by
     exp(-kappa/2), x1.33 to x1.78 across the monotone written shapes.
     The bridge would move it by a further exp(13 T(1)/64), x1.040 to
     x1.047. The anchor's residual is 22 ppm. The anchor takes the
     datum to be the light rows' one loop coupling at M_b, and that
     coupling absorbs every threshold constant at the boundary. So the
     running depends on no shape, and the shape enters only the
     stiffness.
  R3 (the gates). For the labelled bridge every monotone written shape
     fails S, with and without the loop. The best, the window edge,
     leaves Sigma = 4.44 at the boundary, 1.85 short of 2 pi. Of the
     band pass examples at chi(0) = 1, the exponential in |D| and the
     fermionic entropy pass. R5 treats their amplitude. Gate M passes
     under every shape except the sharp cutoff (t0* = 0.589). With the
     loop, M also fails for the bosonic entropy, the Frolov-Fursaev
     masses and three and four Pauli-Villars fields, which already fail
     S. So the labelled bridge's verdict rests on S.
  R4 (what a pass needs). With one rho for every row, the sharp cutoff
     is the floor: w(1) = E[E2(z)] >= E2(E[z]) = E2(1) (Jensen). A pass
     needs the regulator spread over a wide range of proper times, and
     the size of that spread is a parameter the sources leave free.
     Three Pauli-Villars fields at masses squared 1, 100 and 10^4, Li
     and Vaman's count, give 0.465 and pass with the loop. Two levels
     pass at a ratio r >= 13.4 (87.3 with the loop). A one scale onset
     passes if it is slow: P(a, t) for a <= 1.650, outside the proper
     time window, or (t/(1+t))^p for p <= 2.433. A band pass amplitude
     passes below a ceiling: 2.28 for the exponential in |D| (1.63
     with the loop), 1.19 for the fermionic entropy (0.89), 0.994 for
     the heat kernel (0.773), 0.630 for counting (0.538). A finite
     vacuum term (rho below t^2 at small t) excludes the window edge
     and the two level split, which start as t^2. It leaves duality,
     0.280, as the best monotone written shape. (t/(1+t))^p starts as
     t^p, so it keeps the vacuum term finite and still passes.
  R5 (the amplitude). A band pass shape has no log: T(0) = chi(0).
     Take the amplitude from Akama and Terazawa's compositeness
     condition, with no loop of the G2 field's own
     (boundary_content.py B6). The matter rows then induce
     1/alpha_G2 = (b_m/4 pi) chi(0) with b_m = 4, so
     chi(0) = 4 pi/(b_m alpha_G2) = 24 pi^2 = 236.9 in one loop units
     (Chamseddine, Connes and Marcolli's f0 = pi/(8 alpha_G2) = 3 pi^2,
     hep-th/0610241 Eq. 4.10). That is 104 times the ceiling of the
     exponential in |D|. At this amplitude the heat kernel, the
     exponential in |D| and the fermionic entropy give the bridge less
     than 1e-8, while a band pass shape with a second proper time keeps
     it (0.606). With the G2 field's own loop, the induced coupling is
     -(b0/4 pi) chi(0) < 0. B6 treats the coupling as a datum, which
     leaves chi(0) free.

RESULT. No regulator that a source writes down, read by its own
source, gives the bridge rows the weight gate S needs at the boundary.
The monotone examples fail. The best is the window edge at 0.307, and
path integral duality, the one shape a stated principle fixes for every
row and for the running, gives 0.280. The band pass examples pass at
chi(0) = 1, but the spectral action ties chi(0) to the coupling at
24 pi^2, and there the one scale shapes lose the bridge. The sources
fix conditions and counts, not masses, widths or amplitudes, and those
decide the sign. Within the sources' own families, a wide spread of
proper times passes. Read as amount beyond the horizon, the engine's
reading, the bridge is light and the spread does not enter
(tests/probes/horizon_amount.py). The labelled bridge is the reading
whose sign rests on a spread no source fixes.

Usage: python3 tests/probes/boundary_regulator.py
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erfc, exp1, expn, gammainc, kv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

XI = float(r.xi_bridge)
TWO_PI = 2 * math.pi
GAMMA = 0.5772156649015329
ZETA3 = 1.2020569031595942
LOOP = -8.0 / 3.0                                           # ledger_kernel.py R5, on shell
CONTACTS = B.CONTACTS
ONE = B.rows(r.d10, r.d11, r.n7, r.n26, bridge="boundary")  # the labelled bridge at the boundary
LIGHT = {c: B.sigma(B.rows(r.d10, r.d11, r.n7, r.n26, bridge="absent"), 0.0, c, XI) for c in CONTACTS}
BRIDGE = r.n7 * r.n26 * (1 / 6 - XI)                        # the bridge at a massless row's weight
NEED = (TWO_PI - LIGHT[1.0]) / BRIDGE                       # 0.370
NEED_LOOP = (TWO_PI - LIGHT[1.0] - LOOP) / BRIDGE           # 0.462
INF = math.inf


def proper_time(rho):
    """w(x), T(x) and the light row constant kappa = lim (T(x) + ln x) of a shape rho (P2)."""
    tau = 1.0 / quad(lambda t: rho(t) / t**2, 0, INF, limit=200)[0]
    w = lambda x: tau * quad(lambda t: rho(t) * math.exp(-x * t / tau) / t**2, 0, INF, limit=200)[0]
    T = lambda x: quad(lambda t: rho(t) * math.exp(-x * t / tau) / t, 0, INF, limit=200)[0]
    kappa = (quad(lambda t: rho(t) / t, 0, tau, limit=200)[0]
             + quad(lambda t: (rho(t) - 1) / t, tau, INF, limit=200)[0] - GAMMA)
    return w, T, kappa


def pauli_villars(masses2):
    """rho = 1 + sum c_i exp(-a_i t), one field per condition sum c = -1, sum c a^k = 0 (k < n)."""
    a = np.array(masses2, float)
    c = np.linalg.solve(np.vstack([a**k for k in range(len(a))]), np.eye(len(a))[0] * -1.0)
    return lambda t: float(np.sum(c * np.expm1(-a * t)))


def entropy(fn, c):
    """Cutoff function chi(u) = c fn(sqrt u), integrated in v = sqrt(u)."""
    F = lambda y: c * quad(lambda v: 2 * v * fn(v), math.sqrt(y), INF, limit=200)[0]
    f2 = F(0.0)
    return (lambda x: F(x * f2) / f2), (lambda x: c * fn(math.sqrt(x * f2))), f2


def h_fermi(v):
    """Fermionic entropy function (Chamseddine, Connes, van Suijlekom)."""
    return v * math.exp(-v) / (1 + math.exp(-v)) + math.log1p(math.exp(-v))


def k_bose(v):
    """Bosonic entropy function (Dong, Khalkhali, van Suijlekom)."""
    return v * math.exp(-v) / -math.expm1(-v) - math.log(-math.expm1(-v))


def band(name, c):
    """Band pass cutoff functions at chi(0) = c: w(x), T(x)."""
    if name == "heat kernel":
        return (lambda x: math.exp(-c * x)), (lambda x: c * math.exp(-c * x))
    if name == "exponential in |D|":
        return (lambda x: (1 + math.sqrt(2 * c * x)) * math.exp(-math.sqrt(2 * c * x))), \
               (lambda x: c * math.exp(-math.sqrt(2 * c * x)))
    if name == "eigenvalue counting":
        return (lambda x: max(0.0, 1 - c * x)), (lambda x: c if c * x < 1 else 0.0)
    if name == "Gaussian in D^2":
        s = c * math.sqrt(math.pi) / 2
        return (lambda x: float(erfc(s * x))), (lambda x: c * math.exp(-(s * x)**2))
    return entropy(h_fermi, c / math.log(2))[:2]                     # fermionic entropy


SHARP = ((lambda x: float(expn(2, x))), (lambda x: float(exp1(x))), -GAMMA)
WRITTEN = {                                                          # name: (w, T, kappa), kind
    "sharp": (SHARP, "monotone"),
    "duality": (proper_time(lambda t: math.exp(-1 / t)), "monotone"),
    "window edge P(2)": (proper_time(lambda t: float(gammainc(2, t))), "monotone"),
    "three equal PV": (proper_time(lambda t: float(gammainc(3, t))), "monotone"),
    "four equal PV": (proper_time(lambda t: float(gammainc(4, t))), "monotone"),
    "Frolov-Fursaev masses": (proper_time(lambda t: (1 - math.exp(-t))**3 * (1 + math.exp(-t))), "monotone"),
    "bosonic entropy": ((*entropy(k_bose, 2.0)[:2], 2 - math.log(12 * ZETA3)), "monotone"),
    **{k: ((*band(k, 1.0), None), "band") for k in ("eigenvalue counting", "heat kernel", "Gaussian in D^2",
                                                    "exponential in |D|", "fermionic entropy")},
}


def anchored(w, c, extra=0.0):
    """t0 solved from Sigma = Sigma_light + extra + bridge w(Sigma/2 pi); w falls with x, so one root."""
    f = lambda s: LIGHT[c] + extra + BRIDGE * w(s / TWO_PI) - s
    return brentq(f, 1e-9, LIGHT[c] + extra + BRIDGE + 1.0) / TWO_PI if f(1e-9) > 0 else float("nan")


def gates(w):
    sig = {c: LIGHT[c] + BRIDGE * w(1.0) for c in CONTACTS}
    t0, t0l = ({c: anchored(w, c, e) for c in CONTACTS} for e in (0.0, LOOP))
    return dict(sigma=sig, t0=t0, t0_loop=t0l, S=all(v >= TWO_PI for v in sig.values()),
                S_loop=all(v + LOOP >= TWO_PI for v in sig.values()),
                M=all(0.6 < v < 15.0 for v in t0.values()), M_loop=all(0.6 < v < 15.0 for v in t0l.values()))


def run(report=print):
    out = {}

    # R1: what the sources write down.
    table = {k: dict(w1=v[0][0](1.0), T1=v[0][1](1.0), kappa=v[0][2], kind=v[1]) for k, v in WRITTEN.items()}
    expect = {"sharp": 0.148496, "duality": 0.279732, "window edge P(2)": 0.306853, "three equal PV": 0.234721,
              "four equal PV": 0.207456, "Frolov-Fursaev masses": 0.248941, "bosonic entropy": 0.180610,
              "eigenvalue counting": 0.0, "heat kernel": 0.367879, "Gaussian in D^2": 0.210091,
              "exponential in |D|": 0.586936, "fermionic entropy": 0.425603}
    assert all(abs(table[k]["w1"] - v) < 2e-6 for k, v in expect.items())
    assert abs(table["window edge P(2)"]["w1"] - (1 - math.log(2))) < 1e-9
    assert abs(entropy(h_fermi, 1 / math.log(2))[2] - 9 * ZETA3 / (2 * math.log(2))) < 1e-9
    flow = {m: proper_time(lambda t, m=m: float(gammainc(m + 1, m * t)))[0](1.0) for m in (1.0, 1.5, 2.0, 5.0, 40.0)}
    y = 0.5                                                          # per species, three fields, x/(n - 1)
    per_species = 2 * ((1 - y) + y * math.log(y) - (1 - y)**2 / 2)
    assert all(v <= 1 - math.log(2) + 1e-9 for v in flow.values()) and abs(per_species - 0.05685) < 1e-4
    out.update(table=table, flow=flow, per_species=per_species)
    for k, v in table.items():
        report(f"R1 {k:22s} w(1) = {v['w1']:.4f}  T(1) = {v['T1']:.4f}  ({v['kind']})")
    report("   proper time flow m = 1, 1.5, 2, 5, 40: " + ", ".join(f"{v:.4f}" for v in flow.values())
           + f". Pauli-Villars per species (outside P2): {per_species:.4f}")

    # R2: the datum's scheme.
    b_bridge = abs(float(ONE[-1][3]))                                # 13/3, the bridge's share if light
    mono = {k: v for k, v in table.items() if v["kind"] == "monotone"}
    light_f = {k: math.exp(-v["kappa"] / 2) for k, v in mono.items()}
    bridge_f = {k: math.exp(b_bridge * v["T1"] / (2 * float(r.b0_G2))) for k, v in mono.items()}
    x = 1e-14
    assert abs(SHARP[1](x) + math.log(x) + GAMMA) < 1e-9 and abs(table["duality"]["kappa"] + 2 * GAMMA) < 1e-6
    assert abs(WRITTEN["bosonic entropy"][0][1](x) + math.log(x) - table["bosonic entropy"]["kappa"]) < 1e-5
    assert abs(b_bridge / (2 * float(r.b0_G2)) - 13 / 64) < 1e-12
    assert 1.33 < min(light_f.values()) and max(light_f.values()) < 1.79
    assert 1.039 < min(bridge_f.values()) and max(bridge_f.values()) < 1.048
    out.update(light_factor=light_f, bridge_factor=bridge_f)
    report(f"R2 datum at the regulator: light rows M_Pl x{min(light_f.values()):.3f} to x{max(light_f.values()):.3f}, "
           f"the bridge a further x{min(bridge_f.values()):.3f} to x{max(bridge_f.values()):.3f} (anchor 22 ppm). "
           "The datum is the light rows' one loop coupling at M_b, and the shape enters only the stiffness")

    # R3: the gates for the labelled bridge.
    assert abs(B.sigma(ONE, 1.0, 1.0, XI, heavy=0.3) - (LIGHT[1.0] + BRIDGE * 0.3)) < 1e-12
    g = {k: gates(v[0][0]) for k, v in WRITTEN.items()}
    assert all(not g[k]["S"] and not g[k]["S_loop"] for k in mono)
    assert {k for k, v in g.items() if v["S"]} == {"exponential in |D|", "fermionic entropy"}
    assert {k for k, v in g.items() if not v["M"]} == {"sharp"}
    m_loop_fail = {k for k, v in g.items() if not v["M_loop"]}
    assert m_loop_fail == {"sharp", "bosonic entropy", "Frolov-Fursaev masses", "three equal PV", "four equal PV"}
    best = g["window edge P(2)"]["sigma"][1.0]
    assert abs(best - 4.4375) < 1e-3 and abs(g["sharp"]["t0"][1.0] - 0.5885) < 1e-4
    out.update(gates=g, m_loop_fail=m_loop_fail)
    for k, v in g.items():
        report(f"R3 {k:22s} Sigma(M_b) = {v['sigma'][1.0]:+7.3f}  S {'PASS' if v['S'] else 'FAIL'}"
               f" ({'PASS' if v['S_loop'] else 'FAIL'} with the loop)  t0* = {v['t0'][1.0]:.3f}"
               f" ({v['t0_loop'][1.0]:.3f})  M {'PASS' if v['M'] else 'FAIL'} ({'PASS' if v['M_loop'] else 'FAIL'})")

    # R4: the floor and what a pass needs.
    rng = np.random.default_rng(5)
    floor = 1.0
    for _ in range(2000):
        s, p = np.exp(rng.uniform(-5, 5, 4)), rng.dirichlet(np.ones(4))
        t0 = 1 / np.sum(p / s)
        floor = min(floor, t0 * float(np.sum(p * np.array([expn(2, z) for z in s / t0]) / s)))
    pv3 = proper_time(pauli_villars([1.0, 100.0, 1e4]))[0](1.0)
    two = lambda q, need: proper_time(pauli_villars([1.0, q]))[0](1.0) - need
    onset = lambda p: proper_time(lambda t: (t / (1 + t))**p)[0](1.0)
    cross = dict(
        pv_ratio=brentq(two, 1.5, 1e4, args=(NEED,)), pv_ratio_loop=brentq(two, 1.5, 1e4, args=(NEED_LOOP,)),
        gamma_a=brentq(lambda a: proper_time(lambda t: float(gammainc(a, t)))[0](1.0) - NEED, 1.2, 1.99),
        onset_p=brentq(lambda q: onset(q) - NEED, 2.05, 3.0),
        heat=-math.log(NEED), heat_loop=-math.log(NEED_LOOP), counting=1 - NEED, counting_loop=1 - NEED_LOOP)
    for k in ("exponential in |D|", "fermionic entropy"):
        cross[k] = brentq(lambda c: band(k, c)[0](1.0) - NEED, 0.05, 20)
        cross[k + " loop"] = brentq(lambda c: band(k, c)[0](1.0) - NEED_LOOP, 0.05, 20)
    p2, slow = (lambda t: float(gammainc(2, t))), (lambda t: (t / (1 + t))**2.2)
    assert floor >= expn(2, 1.0) - 1e-12 and abs(pv3 - 0.46524) < 1e-4 and abs(onset(2.2) - 0.38566) < 1e-4
    assert abs(cross["pv_ratio"] - 13.38) < 0.02 and abs(cross["pv_ratio_loop"] - 87.32) < 0.05
    assert abs(cross["gamma_a"] - 1.6498) < 1e-3 and abs(cross["onset_p"] - 2.4332) < 1e-3
    assert abs(cross["exponential in |D|"] - 2.283) < 2e-3 and abs(cross["fermionic entropy"] - 1.190) < 2e-3
    assert abs(p2(1e-5) / 1e-10 - 0.5) < 1e-4 and slow(1e-6) / 1e-12 < 0.1   # t^2/2 against t^2.2
    out.update(jensen_floor=floor, pv3_spread=pv3, onset_22=onset(2.2), cross=cross, need=NEED, need_loop=NEED_LOOP)
    report(f"R4 need w(1) >= {NEED:.3f} ({NEED_LOOP:.3f}). 2000 mixtures stay above E2(1) (least {floor:.4f}). PV 1:100:1e4 {pv3:.3f}. "
           f"two PV levels r >= {cross['pv_ratio']:.2f} ({cross['pv_ratio_loop']:.2f}). P(a, t) a <= {cross['gamma_a']:.3f}. "
           f"(t/(1+t))^p p <= {cross['onset_p']:.3f}")
    report("   band pass ceilings: " + ", ".join(f"{k} {cross[k]:.3f}" for k in (
        "exponential in |D|", "exponential in |D| loop", "fermionic entropy", "fermionic entropy loop",
        "heat", "heat_loop", "counting", "counting_loop")))

    # R5: the amplitude from a compositeness condition.
    b_m = float(4 * r.n_dirac_G2 * r.T_7 / 3)
    chi0 = 4 * math.pi / (b_m * r.alpha_G2_Pl)
    f0 = math.pi / (8 * r.alpha_G2_Pl)
    at_chi0 = {k: band(k, chi0)[0](1.0) for k in ("heat kernel", "exponential in |D|", "fermionic entropy")}
    A, lam = chi0 - 0.5, 1e6                                         # chi = A exp(-lam u) + 0.5 exp(-u)
    F = lambda yy: A * math.exp(-lam * yy) / lam + 0.5 * math.exp(-yy)
    two_time = F(F(0.0)) / F(0.0)
    assert b_m == 4 and abs(chi0 - 24 * math.pi**2) < 1e-9 and abs(chi0 / f0 - 8) < 1e-12
    assert max(at_chi0.values()) < 1e-8 and abs(two_time - 0.6061) < 1e-3 and two_time > NEED_LOOP
    assert -float(r.b0_G2) / (4 * math.pi) * chi0 < 0
    out.update(b_m=b_m, chi0=chi0, f0_ccm=f0, band_at_chi0=at_chi0, two_time=two_time)
    report(f"R5 chi(0) = 4 pi/(b_m alpha_G2) = {chi0:.1f} = 24 pi^2 (CCM f0 = {f0:.2f}). "
           + ", ".join(f"{k} {v:.1e}" for k, v in at_chi0.items()) + f". Two proper times {two_time:.3f}")
    report("RESULT: no regulator a source writes down gives the labelled bridge the weight S needs at the boundary when read by "
           f"its own source (best monotone {table['window edge P(2)']['w1']:.3f}, duality {table['duality']['w1']:.3f}). "
           "The sources leave the spread of proper times free, and a wide spread passes.")
    return out


if __name__ == "__main__":
    run()
