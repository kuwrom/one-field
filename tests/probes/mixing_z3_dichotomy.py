"""Why the mixing sector cannot close in its original shape, and where
its premise fails.

The mixing sector asks for CKM and PMNS as RELATIVE DIAGONALIZERS of four
derived sector transports (u, d, e, nu), each a response operator on a
three-dimensional generation carrier.  This probe records four exact facts
about that question.  It consults no mixing data until the last section,
which only counts how many independent numbers the CKM scorecard actually
contains.

  1. DICHOTOMY.  The engine's closure axiom is forgetting: the confined
     sector admits only Z3-invariant operators, [Delta, S] = 0 (the
     closure theorem, test_coverage.py).  Every such operator is a
     circulant, hence diagonal in the Z3 Fourier frame.  Two circulants
     therefore share one eigenframe: the relative diagonalizer of any two
     sector transports built this way is the identity, exactly.  A
     Z3-equivariant engine predicts V = 1 and cannot predict anything else.

  2. RIGIDITY.  The one escape route is the Ocneanu algebra of quantum
     symmetries, whose M3(C) block carries the Delta(27) clock/shift pair.
     Delta(27) contains four Z3 subgroups.  Their eigenbases are pairwise
     MUTUALLY UNBIASED, exactly: every overlap is 1/sqrt(3).  So the
     relative diagonalizer of two Delta(27) frames is either the identity
     (same subgroup) or exactly democratic (different subgroup).  There is
     no interpolation, and no breaking INSIDE Delta(27) can reach a
     hierarchy.  The escape is not a leading-order answer awaiting a
     correction; it is rigid.

  3. PREMISE.  The engine's own mass sector contradicts the question.  The
     circulant amplitude is DERIVED, |B/A| = sqrt(d10), equivalently the
     Koide value Q0 = 1/3 + d10/6 = d10/d11.  Any three masses can be
     written as A + B cos(theta + 2 pi k / 3); what carries content is that
     amplitude.  It holds for (e, mu, tau) and for the two ISOSPIN-CROSSING
     triplets the engine uses to derive m_b and m_s, and it fails for the
     same-charge triplets (u, c, t) and (d, s, b) in every sign branch.
     The engine's Z3 does not sit on per-sector carriers, so the four
     sector transports the mixing sector asks to diagonalize are not objects the engine
     builds.

  4. COUNTING.  The 9 |V_ij| of a unitary 3x3 matrix are a function of
     exactly four real parameters.  The CKM scorecard's 15 entries are
     therefore four independent measurements against four engine numbers,
     not fifteen.

This probe records those facts only.  It proposes no replacement map.
"""

import os
import sys
import itertools
import math
import cmath

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
for _path in (os.path.dirname(_here),
              os.path.join(os.path.dirname(os.path.dirname(_here)),
                           "interference"),
              os.path.join(_here, "interference")):
    if _path not in sys.path:
        sys.path.insert(0, _path)
import _bootstrap  # noqa: E402,F401

from root import d10, d11, Q0  # noqa: E402

_SHIFT = np.roll(np.eye(3), 1, axis=0)
_OMEGA = cmath.exp(2j * math.pi / 3.0)
_CLOCK = np.diag([1.0 + 0j, _OMEGA, _OMEGA**2])


def _fourier_frame():
    """The Z3 Fourier frame: the eigenframe of every circulant."""
    return np.array([[_OMEGA**(i * j) for j in range(3)]
                     for i in range(3)]) / math.sqrt(float(d11))


def same_z3_dichotomy(trials=400, seed=0):
    """Two circulants share one eigenframe: V = 1 exactly."""
    rng = np.random.default_rng(seed)
    frame = _fourier_frame()
    worst_offdiag = 0.0
    worst_mixing = 0.0
    for _ in range(trials):
        pair = []
        for _sector in range(2):
            c = rng.normal(size=3) + 1j * rng.normal(size=3)
            yukawa = (c[0] * np.eye(3) + c[1] * _SHIFT
                      + c[2] * _SHIFT.conjugate().T)
            assert np.allclose(yukawa @ _SHIFT, _SHIFT @ yukawa)
            hermitian = yukawa @ yukawa.conjugate().T
            diagonalized = frame.conjugate().T @ hermitian @ frame
            worst_offdiag = max(worst_offdiag, float(np.abs(
                diagonalized - np.diag(np.diag(diagonalized))).max()))
            pair.append(hermitian)
        # relative diagonalizer of the two sector transports
        worst_mixing = max(worst_mixing,
                           float(np.abs(pair[0] @ pair[1]
                                        - pair[1] @ pair[0]).max()))
    return {"circulant_offdiagonal": worst_offdiag,
            "sector_commutator": worst_mixing}


def delta27_rigidity():
    """The four Z3 subgroups of Delta(27) are exactly mutually unbiased."""
    generators = {
        "S": _SHIFT.astype(complex),
        "C": _CLOCK,
        "CS": _CLOCK @ _SHIFT,
        "CS2": _CLOCK @ _SHIFT @ _SHIFT,
    }
    frames = {}
    for name, generator in generators.items():
        _values, vectors = np.linalg.eig(generator)
        frames[name] = vectors / np.linalg.norm(vectors, axis=0)
    unbiased = 1.0 / math.sqrt(float(d11))
    deviation = 0.0
    for (_a, left), (_b, right) in itertools.combinations(frames.items(), 2):
        overlaps = np.abs(left.conjugate().T @ right)
        deviation = max(deviation, float(np.abs(overlaps - unbiased).max()))
    # the democratic answer the escape produces, and its Jarlskog
    frame = _fourier_frame()
    jarlskog = float((frame[0, 0] * frame[1, 1]
                      * np.conj(frame[0, 1]) * np.conj(frame[1, 0])).imag)
    return {"unbiased_deviation": deviation,
            "overlap": unbiased,
            "democratic_jarlskog": abs(jarlskog),
            "closed_form": 1.0 / (2.0 * d11 * math.sqrt(float(d11)))}


def _circulant_amplitude(masses, signs):
    """|B/A| for sqrt(m_k) = A(1 + (B/A) cos(theta + 2 pi k / 3))."""
    roots = [sign * math.sqrt(mass) for sign, mass in zip(signs, masses)]
    if sum(roots) < 0:
        roots = [-root for root in roots]
    mean = sum(roots) / 3.0
    residues = [root / mean - 1.0 for root in roots]
    harmonic = (2.0 / 3.0) * sum(residues[k] * _OMEGA**(-k) for k in range(3))
    koide = sum(root * root for root in roots) / sum(roots)**2
    return abs(harmonic), cmath.phase(harmonic), koide


def engine_triplet_amplitudes(engine_masses):
    """Which triplets carry the DERIVED amplitude |B/A| = sqrt(d10)?

    engine_masses: dict of the engine's own mass readings, plus the raw
    strange root the Koide solve actually uses (m_s before the bridge).
    """
    derived = math.sqrt(float(d10))
    triplets = {
        "leptons (e,mu,tau)": (("e", "mu", "tau"), (1, 1, 1)),
        "crossing (c,b,t)": (("c", "b", "t"), (1, 1, 1)),
        "crossing (s,c,b)": (("s_raw", "c", "b"), (1, -1, -1)),
        "up-type (u,c,t)": (("u", "c", "t"), None),
        "down-type (d,s,b)": (("d", "s", "b"), None),
    }
    out = {}
    for name, (keys, signs) in triplets.items():
        masses = [engine_masses[key] for key in keys]
        if signs is not None:
            amplitude, theta, koide = _circulant_amplitude(masses, signs)
        else:
            # give the same-charge triplets every sign branch; keep the best
            best = None
            for branch in itertools.product((1, -1), repeat=3):
                if branch[0] < 0:
                    continue
                candidate = _circulant_amplitude(masses, branch)
                if best is None or abs(candidate[0] - derived) < abs(
                        best[0] - derived):
                    best = candidate
            amplitude, theta, koide = best
        out[name] = {"amplitude": amplitude,
                     "theta": theta,
                     "koide": koide,
                     "deviation": abs(amplitude / derived - 1.0)}
    return out


def modular_transition_menu():
    """Every canonical transition amplitude of SU(3)_3 is 0, 1 or 2.

    Two facts, both exact.  (a) The down channel (1,1) is the FIXED POINT
    of the Z3 simple current (3,0), and a simple-current fixed point has
    S_fb = 0 for every b of nonzero triality.  The engine's
    up/down assignment therefore puts the two sectors in exactly
    non-communicating superselection sectors: the up-to-down modular
    amplitude is not small, it is zero.  (b) The whole normalized menu
    S_ab/S_0b out of the up channel takes only the values 0, 1, 2,
    because SU(3)_3 has integer quantum dimensions.  A mixing hierarchy
    needs 0.22, 0.04, 0.004.  No transport in this theory is small.
    """
    import mixing  # local: the probe must not import it at module load

    wzw = mixing._compute_wzw()
    primaries, index, smat = wzw["PRIMARIES"], wzw["IDX"], wzw["S_mat"]
    triality = {p: (p[0] + 2 * p[1]) % 3 for p in primaries}
    up, down = (1, 0), (1, 1)

    charged_overlap = max(
        abs(smat[index[down], index[p]])
        for p in primaries if triality[p] != 0)
    menu = sorted({round(abs(smat[index[up], index[p]])
                         / abs(smat[0, index[p]]), 9)
                   for p in primaries})
    up_down = abs(smat[index[up], index[down]]) / abs(smat[0, index[down]])
    # the orientation-reversed alternative: down = antifundamental
    anti = abs(smat[index[up], index[(0, 1)]]) / abs(smat[0, index[(0, 1)]])
    return {"fixed_point_charged_overlap": charged_overlap,
            "menu": menu,
            "up_to_down": up_down,
            "up_to_antifundamental": anti}


def braiding_phase():
    """pi/d11^2 has a shape: it is the braiding monodromy, read three
    ways that are ONE datum.

    The R-matrix eigenvalue of two fundamentals in the symmetric channel,
    the diagonal modular element arg S_ff, and the identity part of the
    D6* connection q^(2/3) at q = exp(i pi / 6) all equal pi/9.  They are
    NOT three independent witnesses: T_ll = exp(2 pi i (h_l - c/24)), so
    all three are functions of the same conformal weights h_3 = 1/9 and
    h_6 = 5/18, and pi/9 = 2 pi (h_6 - 2 h_3) exactly.  Their agreement
    is a consistency check of the modular data, a good one, and not
    triple confirmation.  It is the one
    CKM coordinate with a derivation shape rather than an
    identification, with the caveat stated in the mixing sector's
    derived_coordinate: what the engine scores against PDG is eta_bar,
    the Wolfenstein APEX COORDINATE, which is magnitude-like.
    """
    import mixing

    wzw = mixing._compute_wzw()
    primaries, index, smat = wzw["PRIMARIES"], wzw["IDX"], wzw["S_mat"]

    def weight(lam):
        return (lam[0]**2 + lam[1]**2 + lam[0] * lam[1]
                + 3 * lam[0] + 3 * lam[1]) / 3.0 / float(d10 * d11)

    fund = (1, 0)
    r_symmetric = cmath.exp(1j * math.pi
                            * (weight((2, 0)) - 2 * weight(fund)))
    return {
        "r_matrix_symmetric": cmath.phase(r_symmetric),
        "arg_S_ff": cmath.phase(smat[index[fund], index[fund]]),
        "connection_identity": cmath.phase(
            cmath.exp(1j * math.pi / float(d10 * d11))**(2.0 / 3.0)),
        "closed_form": math.pi / float(d11**2),
    }


def ckm_scorecard_rank():
    """The 9 |V_ij| are a function of exactly four real parameters."""
    def magnitudes(lam, a_wolf, rho, eta):
        s12 = lam
        s23 = a_wolf * lam**2
        c12 = math.sqrt(1 - s12**2)
        c23 = math.sqrt(1 - s23**2)
        apex = complex(rho, eta)
        z13 = apex * s23 * s12 * c23 / (c12 * (1.0 - apex * s23**2))
        s13 = abs(z13)
        c13 = math.sqrt(1 - s13**2)
        phase = cmath.exp(1j * cmath.phase(z13))
        return np.abs(np.array([
            [c12 * c13, s12 * c13, s13],
            [-s12 * c23 - c12 * s23 * s13 * phase,
             c12 * c23 - s12 * s23 * s13 * phase, s23 * c13],
            [s12 * s23 - c12 * c23 * s13 * phase,
             -c12 * s23 - s12 * c23 * s13 * phase, c23 * c13],
        ])).ravel()

    # a GENERIC point: the rank is a property of the parameterization,
    # not of any measured value, so this probe stays blind
    point = np.array([0.2, 0.8, 0.15, 0.35])
    step = 1e-6
    jacobian = np.array([
        (magnitudes(*(point + step * np.eye(4)[i])) - magnitudes(*point)) / step
        for i in range(4)])
    return int(np.linalg.matrix_rank(jacobian, tol=1e-6))


def run(report=print, engine_masses=None):
    report("MIXING SECTOR Z3 DICHOTOMY: what the question can return")
    report("=" * 62)

    dichotomy = same_z3_dichotomy()
    assert dichotomy["circulant_offdiagonal"] < 1e-9
    assert dichotomy["sector_commutator"] < 1e-9
    report("  1. Selector theorem [Delta,S]=0 -> circulant -> one shared")
    report("     Fourier frame.  Two sector transports commute exactly;")
    report(f"     max sector commutator = "
           f"{dichotomy['sector_commutator']:.1e}.")
    report("     A Z3-equivariant engine predicts V = 1 and nothing else.")

    rigidity = delta27_rigidity()
    assert rigidity["unbiased_deviation"] < 1e-12
    assert abs(rigidity["democratic_jarlskog"]
               - rigidity["closed_form"]) < 1e-12
    report("  2. The four Z3 subgroups of Delta(27) are EXACTLY mutually")
    report(f"     unbiased (max deviation from 1/sqrt(3): "
           f"{rigidity['unbiased_deviation']:.1e}).  The Ocneanu escape")
    report(f"     returns |V_ij| = 1/sqrt(3), |J| = "
           f"{rigidity['democratic_jarlskog']:.6f} = 1/(2 d11 sqrt(d11)),")
    report("     with no dial between it and the identity.")

    triplets = None
    if engine_masses is not None:
        triplets = engine_triplet_amplitudes(engine_masses)
        crossing = [v for k, v in triplets.items() if "type" not in k]
        same_charge = [v for k, v in triplets.items() if "type" in k]
        assert all(v["deviation"] < 0.01 for v in crossing)
        assert all(v["deviation"] > 0.05 for v in same_charge)
        report(f"  3. Derived circulant amplitude |B/A| = sqrt(d10) = "
               f"{math.sqrt(float(d10)):.5f}  (Koide Q0 = {float(Q0)}):")
        for name, value in triplets.items():
            report(f"       {name:<20} |B/A| = {value['amplitude']:.5f}"
                   f"   dev {value['deviation']*100:6.2f}%")
        report("     It holds on the ISOSPIN-CROSSING triplets the engine")
        report("     uses to derive m_b and m_s, and fails on the")
        report("     same-charge triplets in every sign branch.  The engine")
        report("     builds no per-sector Z3 carrier for the mixing sector to diagonalize.")

    menu = modular_transition_menu()
    assert menu["fixed_point_charged_overlap"] < 1e-12
    assert menu["up_to_down"] < 1e-12
    assert all(abs(value - round(value)) < 1e-9 for value in menu["menu"])
    assert set(round(value) for value in menu["menu"]) <= {0, 1, 2}
    report("  4. The down channel (1,1) is the fixed point of the Z3 simple")
    report("     current (3,0), so its overlap with every charged sector")
    report(f"     vanishes: max |S_(1,1),b| = "
           f"{menu['fixed_point_charged_overlap']:.1e}.  The up/down")
    report("     assignment puts them in non-communicating")
    report("     superselection sectors: the amplitude is zero, not small.")
    report(f"     The whole normalized menu out of the up channel is "
           f"{[int(round(v)) for v in menu['menu']]}.")
    report(f"     Orientation-reversed alternative (down = antifundamental)"
           f" gives {menu['up_to_antifundamental']:.4f}, maximal.")
    report("     No transport amplitude in this theory is small; a CKM")
    report("     hierarchy needs 0.22, 0.04, 0.004.")

    phase = braiding_phase()
    for key in ("r_matrix_symmetric", "arg_S_ff", "connection_identity"):
        assert abs(phase[key] - phase["closed_form"]) < 1e-12
    report(f"  5. pi/d11^2 = {phase['closed_form']:.6f}, one datum read "
           f"three ways (all")
    report("     functions of the same conformal weights, not "
           "independent witnesses):")
    report("     the R-matrix eigenvalue of two fundamentals in the")
    report("     symmetric channel, arg S_ff, and the identity part of the")
    report("     D6* connection q^(2/3).  eta_bar is a BRAIDING phase, the")
    report("     one CKM coordinate with a derivation shape.")

    rank = ckm_scorecard_rank()
    assert rank == 4
    report(f"  6. The 9 |V_ij| have parameter rank {rank}: the 15-entry CKM")
    report("     scorecard is four independent measurements, not fifteen.")

    return {
        "sector_commutator": dichotomy["sector_commutator"],
        "circulant_offdiagonal": dichotomy["circulant_offdiagonal"],
        "unbiased_deviation": rigidity["unbiased_deviation"],
        "democratic_jarlskog": rigidity["democratic_jarlskog"],
        "triplet_amplitudes": triplets,
        "modular_menu": menu,
        "braiding_phase": phase,
        "ckm_parameter_rank": rank,
    }


if __name__ == "__main__":
    import contextlib
    import io

    import root as _root
    import masses as _masses

    with contextlib.redirect_stdout(io.StringIO()):
        _R = _root.derive()
        _M = _masses.derive(_R)
    _reading = {k: _M[f"m_{k}"] for k in
                ("e", "mu", "tau", "u", "d", "s", "c", "b", "t")}
    _reading["s_raw"] = _M["m_s"] / _M["bridge"]
    run(engine_masses=_reading)
