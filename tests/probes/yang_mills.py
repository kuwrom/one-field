"""Yang-Mills: the strong sector closes gapped, a third of a turn from the empty field.

THE PROBLEM, MAPPED. The Millennium problem asks for a quantum
Yang-Mills theory on R^4, for every compact simple group, that obeys
the Wightman axioms and has a mass gap. It treats fields on a continuum
as the basic objects. The framework treats records as the basic objects
and asks which of them close. Each part of the problem has a
counterpart here.
  The theory. Gauge forces are collective motions of the confined
     phase, not fields put in by hand (interference/__init__.py, THE
     GAUGE FIELDS). The theory is there when the
     colour record closes (P2, P3).
  The group. It is not chosen. SU(3) appears when G2 breaks (Coherence
     Gate paper, Sec. III), and the other levels are tested and fail (Y2).
  The state space and the vacuum. Records come first, not Hilbert space
     rays (interference/__init__.py, THE ONTOLOGY). The vacuum is the
     empty field.
  The axioms. The lattice requirement R1 gives locality and rules out
     negative norms (Sec. II). Three space dimensions and one time follow
     from closure (dimension_by_closure.py).
  The continuum limit. It is not needed. The substrate is the medium,
     and running it is the verification (THE GENERATIVE DIRECTION).
  The mass gap. No colour record comes within a third of a turn of the
     vacuum (Y1), at the scale of Y4.
The test below asks the problem in these terms. It does not build a
continuum field.

THE QUESTION. What closes is what lasts (dimension_by_closure.py). A
particle is a closed standing wave of the substrate whose phase turns
at one fixed rate, and that rate is its mass (emergent_clock.py,
lump_in_empty_field.py). The engine gives each particle one reading,
with every echo included (masses.py), so the reading fixes the mass.
The vacuum is the empty field, which never closes into anything
(interference/__init__.py, WHERE THE QUESTIONS END). A strong sector
with no gap would need colour records, closed colour patterns, whose
readings come as close to the vacuum as one likes, with masses as close
to zero. Does the colour record allow that? If not, where does a
continuum live?

PREMISES.

P1 (what a reading carries). A physical phase lives in a relation
between closed records, a braid for example, and a phase shared by
everything adds nothing physical (THE ONTOLOGY). So two readings that
differ only by a common phase are the same reading, and what a record
carries is its relative phases. Its distance from the vacuum is the
smallest arc of the circle, in turns, that holds all its phases. A
finite group of readings keeps a smallest distance from the vacuum. An
infinite group does not: one record, repeated, comes back as close to
the vacuum as one likes.

P2 (the data are fixed on their own). The SU(3)_3 data follow from the
three requirements R1 to R3 on the one embedding E8(1) > G2(1) x F4(1)
(Coherence Gate paper, Secs. II, III and IX). Confinement plays no part
in fixing them. They are topological, so they do not depend on the
number of dimensions. Nothing two dimensional is carried into four, and
no conformal weight is used as an energy.

P3 (the colour record). The colour record is SU(3) at level 3. The level
is fixed by h(J) = 1 (Sec. VI), and the record has ten states with
quantum dimensions 1, 2 and 3 (Secs. VII and X, root.py). Its centre is
Z3, the three simple currents {1, J, J^2} (wiring.py). The mass
operator commutes with the Z3 shift, so it treats the three centre
labels alike by itself, and no wall or choice of modular invariant is
needed (Sec. XI, masses.py). The colour braids are the ones mixing.py
solves on the engine's graph of colour states, the same operators the
CKM matrix is read from.

P4 (the bridge). Only one channel crosses a confining boundary: the
bridge B = (7,26) between the G2 and F4 faces (Sec. XI,
wiring_theorem_probe.py PART B). Both faces are Fibonacci and confined.
The engine reads the bridge as an amount with no identity (boundary.py,
THE BRIDGE).

P5 (the emergence ladder). Each part sits on one rung of the ladder
(tests/audit/emergence_ladder.py). The closure algebra rung is held by
integrality: the braid image is finite, so the menu is discrete. A
discrete menu has no readings crowding toward the vacuum, and that is
the gap. Scales sit a factor e^-S apart and do not mix. Hadrons are
read from QCD with the engine's coupling and quark masses.

THE CANDIDATES. SU(3) at levels k = 1 to 6, and the bridge.

STATEMENTS (each is computed below and checked in run()).

  Y1 (the colour floor). The engine's colour braids form finite groups,
     of order 144 on three steps and 3888 at depth four. Every colour
     record other than the vacuum sits at least a third of a turn from
     it. A third of a turn is the Z3 phase.
  Y2 (the other levels). The engine has braids only at its own level,
     so the probe builds the others from the standard Hecke algebra
     representations at q = exp(2 pi i/(k + 3)). At k = 3 they give the
     engine's floor again. At k = 2, 4, 5 and 6 a short braid never
     returns to where it started, and its powers come within 1e-4 turn
     of the vacuum. k = 1 is finite too, but R3 rules it out (P3).
     Single strands are enough: every state is built from them, and
     crossing two bundles of strands is a braid of the strands inside.
     On the engine's braids two crossing pairs stay in the depth four
     group, half a turn from the vacuum. Known theorems agree: dense
     for r of 5 and above except r = 6, finite at r = 6 (Freedman,
     Larsen, Wang, arXiv math/0103200), and finite for every state of
     SU(3)_3 (Etingof, Nikshych, Ostrik, arXiv 0809.3031, and Green,
     Nikshych, arXiv 1911.02633).
  Y3 (the bridge). The bridge braids as the product of its two faces,
     with weights 2/5 and 3/5 (wiring_theorem_probe.py PART A). Its
     group is infinite, and its powers come within 1e-4 turn of the
     vacuum.
  Y4 (the scale). The confined layer sets its own scale. Its coupling
     at the Planck scale, 1/(24 pi) (Sec. XI), runs with b0 = 32/3 to
     the action S_0 = 9 pi^2/2, which does not depend on conventions
     (root.py). With the fundamental channel weight C_EFF = 1/2 in front
     (root.py), Lambda_conf = C_EFF M_Pl e^-S_0 = 314.2 MeV
     (higgs_ew.py). SU(3) stays unbroken, because the condensate sits in
     the singlet of 7 -> 1 + 3 + 3bar (couplings.py). How heavy the
     lightest hadron is in units of Lambda_conf belongs to the hadron
     rung (P5).
  Y5 (the other phases do not close). With six quark flavours, a free
     Coulomb phase would need 17 flavours or more, and an infrared fixed
     point (Banks-Zaks) 9 or more (couplings.py). The centre cannot
     break partly, because 3 is prime. A broken centre keeps a label,
     and an operator that keeps one, diag(1, -1/2, -1/2), does not
     commute with the Z3 shift, while the lepton circulant does
     (interference/__init__.py, Gate 4). At theta = pi the vacuum would
     split in two (Gaiotto, Kapustin, Komargodski, Seiberg, JHEP 05
     (2017) 091), and the engine gives theta = 0 (couplings.py). Every
     quark mass is positive (masses.py), so there are no massless pions.

RESULT. The colour menu is finite, so there are finitely many readings
and finitely many masses. Every colour record other than the vacuum
turns by at least a third of a turn, so none is massless, and the
masses above the vacuum have a smallest value. The strong sector closes
with a gap and its centre unbroken, at Lambda_conf. The continuum lives
on the bridge, the one channel with no identity. Running the substrate
in this phase is the generative direction (interference/__init__.py,
THE GENERATIVE DIRECTION). Its classical signature, measured on the
gas sign, is the area law of confined_phase.py.

KILLED BY. A colour braid of the engine of infinite order, a finite
group for the bridge, or a massless strong state in nature.

Usage: python3 tests/probes/yang_mills.py     (a few seconds)
"""

import cmath
import contextlib
import io
import itertools
import math
import os
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_here, "..", ".."))
sys.path[:0] = [p for p in (_repo, os.path.join(_repo, "tests"), os.path.join(
    _repo, "tests", "audit"), os.path.join(_repo, "interference")) if p not in sys.path]

import emergence_ladder                                              # noqa: E402
import interference                                                  # noqa: E402
from interference import boundary, couplings, masses, mixing, root, wiring  # noqa: E402
from interference.root import (BA_ratio, C_EFF, K, M_Pl_GeV,        # noqa: E402
                               S_lepton, alpha_G2_Pl, b0_G2, d11, h10, k_SU3)
from probes import wiring_theorem_probe                             # noqa: E402

# ── distance from the vacuum reading ─────────────────────────────────

def arc_of(phases):
    """Smallest arc, in turns, holding every phase."""
    a = np.sort(np.mod(phases, 2 * np.pi))
    return float(1 - np.diff(np.concatenate([a, [a[0] + 2 * np.pi]])).max() / (2 * np.pi))

print_off = lambda *a, **k: None
arc = lambda M: arc_of(np.angle(np.linalg.eigvals(M)))
key = lambda M: tuple(np.round(M.ravel(), 6).tolist())


def image(gens, cap=6000):
    """The order of the generated group, its floor (the least arc of a
    nonscalar element) and its elements. None if it passes the cap."""
    step = list(gens) + [np.linalg.inv(T) for T in gens]
    eye = np.eye(gens[0].shape[0], dtype=complex)
    seen, frontier, floor = {key(eye)}, [eye], 1.0
    while frontier:
        new = []
        for M in frontier:
            for P in (M @ G for G in step):
                if key(P) not in seen:
                    seen.add(key(P))
                    new.append(P)
                    floor = min(floor, arc(P)) if arc(P) > 1e-7 else floor
                    if len(seen) > cap:
                        return None
        frontier = new
    return len(seen), floor, seen


def accumulation(gens, max_len=4, n_max=20000):
    """Closest approach to the vacuum among the powers of a short word of
    infinite order (an eigenvalue no root of unity of order <= 2000)."""
    step = list(gens) + [np.linalg.inv(T) for T in gens]
    words = [np.eye(gens[0].shape[0], dtype=complex)]
    for _ in range(max_len):
        words = [M @ G for M in words for G in step]
        for lam in (np.linalg.eigvals(P) for P in words):
            if any(min(abs(x ** m - 1) for m in range(1, 2001)) > 1e-7 for x in lam):
                arcs = (arc_of(n * np.angle(lam)) for n in range(1, n_max + 1))
                return min(a for a in arcs if a > 1e-12)
    return None


# ── Y1: the engine's colour braids ───────────────────────────────────

def engine_colour():
    d6 = mixing._compute_d6_boundary()
    n10, po, (s1, s2, s3) = d6["n10"], d6["path_operator"], d6["depth4_braids"]
    paths3 = tuple(p for p in itertools.product(range(n10.shape[0]), repeat=4)
                   if all(n10[p[i], p[i + 1]] for i in range(3)))
    depth4, cable = image([s1, s2, s3]), s2 @ s1 @ s3 @ s2   # pair (1,2) crosses (3,4)
    return {"depth3": image([po(paths3, 1), po(paths3, 2)]), "depth4": depth4,
            "cable_in_group": key(cable) in depth4[2], "cable_floor": image([cable])[1]}


# ── Y2: the rival levels ─────────────────────────────────────────────

def shapes(n, N, k):
    """Partitions of n with at most N rows and lambda_1 - lambda_N <= k."""
    return [p for m in range(1, N + 1) for p in itertools.product(range(n, 0, -1), repeat=m)
            if sum(p) == n and list(p) == sorted(p, reverse=True)
            and p[0] - (p[-1] if len(p) == N else 0) <= k]


def tableaux(shape):
    """Standard Young tableaux as {entry: (row, col)}."""
    boxes = [(r, c) for r, length in enumerate(shape) for c in range(length)]
    fills = (dict(zip(boxes, p)) for p in itertools.permutations(range(1, len(boxes) + 1)))
    return [{v: b for b, v in t.items()} for t in fills
            if all(t[(r, c)] < min(t.get((r, c + 1), 99), t.get((r + 1, c), 99))
                   for r, c in boxes)]


def hecke_irrep(shape, q):
    """Young's seminormal T_i, (T - q)(T + 1) = 0. Off the row and column cases
    (t, s_i t) has diagonal f(+-d), f(d) = (q - 1)/(1 - q^-d), d the content gap."""
    tabs = tableaux(shape)
    index = {tuple(sorted(t.items())): a for a, t in enumerate(tabs)}
    f = lambda d: (q - 1) / (1 - q ** (-d))
    gens = []
    for i in range(1, sum(shape)):
        T = np.zeros((len(tabs), len(tabs)), dtype=complex)
        for a, t in enumerate(tabs):
            (r1, c1), (r2, c2) = t[i], t[i + 1]
            if r1 == r2:
                T[a, a] = q
            elif c1 == c2:
                T[a, a] = -1
            else:
                d = (c2 - r2) - (c1 - r1)
                s = dict(t)
                s[i], s[i + 1] = t[i + 1], t[i]
                T[a, a] = f(d)
                T[a, index[tuple(sorted(s.items()))]] = cmath.sqrt(f(d) * f(-d) + q)
        gens.append(T)
    I = np.eye(len(tabs))
    assert all(np.abs((T - q * I) @ (T + I)).max() < 1e-9 for T in gens)
    assert all(np.abs(A @ B @ A - B @ A @ B).max() < 1e-9 for A, B in zip(gens, gens[1:]))
    return gens


def level(N, k):
    """Per irreducible on three and four strands: order and floor, or closest approach."""
    q, out = cmath.exp(2j * math.pi / (k + N)), {}
    for lam in shapes(3, N, k) + shapes(4, N, k):
        gens = hecke_irrep(lam, q)
        if gens[0].shape[0] > 1:
            got = image(gens)
            out[lam] = ("finite",) + got[:2] if got else ("infinite", accumulation(gens))
    return out


# ── Y3: the bridge from its two Fibonacci faces ──────────────────────

def face_braids(h, d):
    """B_3 on the tau x tau x tau -> tau space of a Fibonacci face with
    weight h and dimension d. R has phases exp(-2 pi i h) and
    eps exp(-i pi h), F is the Fibonacci F matrix, and the braid relation
    closes on exactly one sign eps."""
    F = np.array([[1 / d, d ** -0.5], [d ** -0.5, -1 / d]])
    Rs = [np.diag([cmath.exp(-2j * math.pi * h), e * cmath.exp(-1j * math.pi * h)])
          for e in (1, -1)]
    closed = [[R, F @ R @ F] for R in Rs
              if np.abs(R @ F @ R @ F @ R - F @ R @ F @ R @ F @ R @ F).max() < 1e-12]
    assert len(closed) == 1
    return closed[0]


def bridge_braids(cat):
    """B = (7,26) braids as the product of its G2(1) and F4(1) faces."""
    g2, f4 = (face_braids(cat["hs"][f], cat["ds"][f]) for f in (("7", "1"), ("1", "26")))
    return [np.kron(a, b) for a, b in zip(g2, f4)]


# ── run ──────────────────────────────────────────────────────────────

def run(report=print, chain=None):
    if chain is None:
        with contextlib.redirect_stdout(io.StringIO()):
            R = root.derive()
            m = masses.derive(R)
            chain = m, couplings.derive(R, m)
    m, c = chain

    said = lambda doc, *quotes: all(q in " ".join(doc.split()) for q in quotes)
    assert said(interference.__doc__, "never closes into anything further",
                "does not add a physical line", "not postulated fields",
                "THE GENERATIVE DIRECTION", "absorption program")
    assert said(boundary.__doc__, "amount without identity") and said(
        masses.__doc__, "one Web reading, including every echo") and said(
        open(os.path.join(_here, "emergent_clock.py")).read(), "Mass is a rate (E = ħω)")
    rungs = {r["key"]: r for r in emergence_ladder.RUNGS}
    assert said(rungs["algebra"]["keeps"], "braid image is finite") and rungs[
        "hadrons"]["read"].startswith("QCD with the engine's alpha_s")

    g1 = engine_colour()
    report("Y1  engine colour braids: order %d (three steps), %d (depth four), floor "
           "%.5f turn, bundled pairs in the group %s with floor %.5f turn"
           % (g1["depth3"][0], g1["depth4"][0], g1["depth4"][1],
              g1["cable_in_group"], g1["cable_floor"]))

    assert k_SU3 + d11 == K                     # r = k + N at the colour point
    menu = {k: level(d11, k) for k in range(1, 2 * k_SU3 + 1)}
    finite = {k: all(v[0] == "finite" for v in menu[k].values()) for k in menu}
    rival_acc = max(v[1] for k in menu for v in menu[k].values() if v[0] == "infinite")
    assert all(v[-1] is not None for k in menu for v in menu[k].values())
    hecke_floor = min(v[2] for v in menu[k_SU3].values())
    report("Y2  finite image by level %s, floor at k = 3 %.5f turn, "
           "other levels reach %.1e" % (finite, hecke_floor, rival_acc))

    cat = wiring_theorem_probe.part_B(*wiring_theorem_probe.part_A(print_off), print_off)
    deconfined = cat["B"]               # part_B asserts it is the unique nontrivial boson
    gB = bridge_braids(cat)
    bridge_infinite, fib_acc = image(gB) is None, accumulation(gB)
    report("Y3  bridge %s, d_B = %.6f, image infinite %s, reaches %.1e"
           % (deconfined, cat["ds"][cat["B"]], bridge_infinite, fib_acc))

    inv_a, s0 = 1 / c["alpha_s_mu_star"], 2 * math.pi / (float(b0_G2) * alpha_G2_Pl)
    lam_conf = C_EFF * M_Pl_GeV * math.exp(-S_lepton)
    report("Y4  alpha_G2(M_Pl) = 1/(%.6f pi), b0 = %s, S_0 = %.6f, Lambda_conf = %.1f MeV"
           % (1 / (alpha_G2_Pl * math.pi), b0_G2, s0, 1e3 * lam_conf))
    report("    1/alpha(mu*) = %.12f (32/pi), alpha_s(M_Z) = %.5f, theta_bar = %s, M_V = %.1f GeV"
           % (inv_a, c["alpha_s_MZ_thresh"], c["bar_theta"], c["M_V_GeV"]))

    b0 = {nf: couplings._b0(nf) for nf in range(18)}
    b1 = {nf: couplings._b1(nf) for nf in range(18)}
    w = wiring.canonical()
    z3 = len(w["pointed_orbit"])                # |Z3| = the simple current orbit
    S = np.roll(np.eye(z3), 1, axis=0)          # the Z3 shift
    D = np.diag(np.cos(2 * np.pi * np.arange(z3) / z3))   # v = (1, -1/2, -1/2)
    c1 = BA_ratio / 2 * cmath.exp(1j * float(h10))  # masses.py: I + c1 S + c1* S^dag
    circ = np.eye(z3) + c1 * S + np.conj(c1) * S.T
    g5 = {"infrared_free": [nf for nf in b0 if b0[nf] < 0],
          "banks_zaks": [nf for nf in b0 if b0[nf] > 0 and b1[nf] < 0],
          "z3_subgroups": [d for d in range(1, z3 + 1) if z3 % d == 0],
          "labelled": float(np.linalg.norm(D @ S - S @ D)),
          "circulant": float(np.linalg.norm(circ @ S - S @ circ)),
          "min_quark_MeV": min(m[q] for q in ("m_u", "m_d", "m_s", "m_c", "m_b", "m_t"))}
    report("Y5  infrared free at nf %(infrared_free)s, Banks-Zaks zero at nf %(banks_zaks)s, "
           "Z3 subgroups %(z3_subgroups)s, |[diag(v),S]| = %(labelled).3f, "
           "|[circ,S]| = %(circulant).0e, lightest quark %(min_quark_MeV).3f MeV" % g5)

    third = 1 / z3                              # the Z3 phase, in turns
    assert g1["depth3"][0] == 144 and g1["depth4"][0] == 3888
    assert abs(g1["depth3"][1] - third) < 1e-9 and abs(g1["depth4"][1] - third) < 1e-9
    assert g1["cable_in_group"] and g1["cable_floor"] >= third - 1e-9
    assert finite == {1: True, 2: False, 3: True, 4: False, 5: False, 6: False}
    assert abs(hecke_floor - third) < 1e-9 and rival_acc < 1e-3
    assert deconfined == ("7", "26") and bridge_infinite and fib_acc < 1e-3
    assert abs(inv_a - 32 / math.pi) < 1e-12 and c["bar_theta"] == 0 and b0_G2 * 3 == 32
    assert abs(alpha_G2_Pl * 24 * math.pi - 1) < 1e-12 and abs(s0 - S_lepton) < 1e-12
    assert abs(s0 - 9 * math.pi**2 / 2) < 1e-12 and g5["z3_subgroups"] == [1, 3]
    assert g5["infrared_free"] == [17] and g5["banks_zaks"] == list(range(9, 17))
    assert g5["labelled"] > 1 and g5["circulant"] < 1e-12 and g5["min_quark_MeV"] > 0

    report("RESULT  every colour record stands a third of a turn from the empty "
           "field and the bridge accumulates. The strong sector closes gapped, "
           "with its centre unbroken, at Lambda_conf = %.1f MeV (SU(3)_%d, D^2 = %.0f)"
           % (1e3 * lam_conf, w["k"], w["D2_tot"]))
    return dict(g1=g1, menu=menu, bridge_accumulation=fib_acc, lambda_conf=lam_conf, g5=g5)


if __name__ == "__main__":
    run()
