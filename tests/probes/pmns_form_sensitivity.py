"""Local sensitivity audit for the committed PMNS rotation ansatz.

It runs the declared 144-row Cartesian scan around the committed
ansatz. It then audits nearby omitted options and convention choices.  The
scan measures local sensitivity of the ansatz.  The frozen NuFit
marginals are used only for comparison, once, last.

Declared axes (144 = 3*3*2*2*2*2): rotation planes, order, sign(phi),
theta reading {h10, lambda}, and phi scale {total, per mode}, with the
TBM convention and delta = acos(-lambda) fixed.  The phi scale pair and
the form axes are the sensitivity axes of the scan.

RESULT. Of the 144 declared rows, 2 lie within 1 sigma of the NuFIT
marginals, and no theta sign or label permutation turns a miss into a
hit. The committed R13 rotation phase arccos(-lambda) = 103.02 deg lies
inside the 0.5 sigma solar band. It is not delta_CP: the engine carries
the sign -1 on this phase, whose rephasing invariant is delta_CP =
275.74 deg.
"""

import cmath
import itertools
import math

import numpy as np

import _bootstrap  # noqa: E402,F401

from root import d10, K, h10, h11, BA_ratio, NUFIT_PMNS  # noqa: E402


PLANES = {"12": (0, 1), "13": (0, 2), "23": (1, 2)}


def _rot(plane, angle, phase=None):
    """Unitary rotation in one labeled plane."""
    M = np.eye(3, dtype=complex)
    i, j = PLANES[plane]
    c, s = math.cos(angle), math.sin(angle)
    e = 1.0 if phase is None else cmath.exp(1j * phase)
    M[i, i] = M[j, j] = c
    M[i, j] = s * e
    M[j, i] = -s * (e.conjugate() if phase is not None else 1.0)
    return M


def _observables(U):
    s13 = abs(U[0, 2])
    c13_sq = 1.0 - s13**2
    return {
        "12": abs(U[0, 1])**2 / c13_sq,
        "23": abs(U[1, 2])**2 / c13_sq,
        "13": s13**2,
    }


def _pulls(U):
    obs = _observables(U)
    return {k: (obs[k] - NUFIT_PMNS[f"sin2_{k}"][0])
            / NUFIT_PMNS[f"sin2_{k}"][1] for k in obs}


def _hit(pulls, tolerance):
    return max(abs(x) for x in pulls.values()) <= tolerance


def _scan(U_TBM, r, delta, theta_readings, theta_signs=(1,)):
    rows = []
    for p_phi, p_theta in itertools.product(PLANES, repeat=2):
        for order in ("phi_first", "theta_first"):
            for sign_phi in (+1, -1):
                for theta_name, theta in theta_readings:
                    for phi_name, phi in (("total", r),
                                          ("per_mode", r / math.sqrt(float(d10)))):
                        for sign_theta in theta_signs:
                            R_phi = _rot(p_phi, sign_phi * phi)
                            R_theta = _rot(p_theta, sign_theta * theta, delta)
                            U_e = (R_phi @ R_theta if order == "phi_first"
                                   else R_theta @ R_phi)
                            U = U_e.conj().T @ U_TBM
                            key = (p_phi, p_theta, order, sign_phi,
                                   theta_name, phi_name, sign_theta)
                            rows.append((key, U, _pulls(U)))
    return rows


def _fingerprint(U):
    obs = _observables(U)
    return tuple(round(obs[k], 13) for k in ("12", "23", "13"))


def run(report=print):
    from mixing import _compute_wzw, _derive_pmns

    # Engine inputs: the circulant mass ratio, exact echoed lambda, and
    # the engine's own TBM block.  No reference data enters this stage.
    gaps = {k: 1.0 + float(BA_ratio) * math.cos(float(h10)
                                                + 2 * math.pi * k / 3)
            for k in range(3)}
    mmu_over_me = (gaps[2] / gaps[1])**2
    r = math.sqrt(1.0 / mmu_over_me)
    lam = float(h10 * (1 + h11 / K**2))
    delta = math.acos(-lam)
    pm = _derive_pmns(_compute_wzw(), lam, 1.0, mmu_over_me)
    U_TBM = pm["U_TBM"]

    declared_theta = (("h10", float(h10)), ("lambda", lam))
    declared = _scan(U_TBM, r, delta, declared_theta)
    assert len(declared) == 144

    committed = ("23", "13", "phi_first", -1,
                 "h10", "per_mode", 1)
    committed_row = next(row for row in declared if row[0] == committed)
    for angle in ("12", "23", "13"):
        assert abs(committed_row[2][angle] - pm["pulls"][angle]) < 1e-9

    declared_hits = {tol: [row for row in declared if _hit(row[2], tol)]
                   for tol in (1.0, 0.5)}
    unique_fingerprints = {_fingerprint(U) for _, U, _ in declared}

    # Same plane rows are reducible to a single SU(2) plane up to phase
    # conventions.  Removing them changes 144 -> 96 but not the hits.
    distinct_planes = [row for row in declared if row[0][0] != row[0][1]]
    distinct_hits = {tol: sum(_hit(row[2], tol) for row in distinct_planes)
                     for tol in (1.0, 0.5)}

    # The engine source names arctan(lambda) as the nearby theta fork.
    # arcsin(lambda) is included as a second standard angle map.
    extended_theta = declared_theta + (("atan_lambda", math.atan(lam)),
                                     ("asin_lambda", math.asin(lam)))
    extended = _scan(U_TBM, r, delta, extended_theta)
    extended_hits = {tol: sum(_hit(row[2], tol) for row in extended)
                     for tol in (1.0, 0.5)}

    # Compare the three nearby theta readings at the committed form.
    # The declared scan calls theta = lambda the echoed fork. The engine
    # source names atan(lambda) instead. Both definitions are kept explicit.
    def committed_with_theta(theta):
        phi = r / math.sqrt(float(d10))
        U_e = _rot("23", -phi) @ _rot("13", theta, delta)
        return U_e.conj().T @ U_TBM

    fork_comparison = {}
    for name, theta in (("h10", float(h10)),
                        ("lambda", lam),
                        ("atan_lambda", math.atan(lam))):
        pulls = _pulls(committed_with_theta(theta))
        fork_comparison[name] = {
            "pulls": pulls,
            "max_abs_pull": max(abs(x) for x in pulls.values()),
            "diagonal_score": sum(x*x for x in pulls.values()),
        }

    # delta -> -delta complex conjugates the construction because TBM
    # and the phi rotation are real. The three angle observables agree.
    conjugate = _scan(U_TBM, r, -delta, declared_theta)
    assert all(_fingerprint(a[1]) == _fingerprint(b[1])
               for a, b in zip(declared, conjugate))

    # Adding the theta sign doubles the raw rows. Record whether a negative
    # theta row turns a nonhit of the declared scan into a hit, rather
    # than only adding conventions.
    signed = _scan(U_TBM, r, delta, declared_theta, theta_signs=(1, -1))
    declared_hit_keys = {row[0][:-1] for row in declared_hits[1.0]}
    theta_sign_promotions = {
        row[0][:-1] for row in signed
        if row[0][-1] == -1 and _hit(row[2], 1.0)
        and row[0][:-1] not in declared_hit_keys
    }

    # Check all physical label permutations as a diagnostic, without
    # counting them as grammar generated alternatives.
    perms = tuple(itertools.permutations(range(3)))
    permutation_promotions = set()
    for key, U, pulls in declared:
        if _hit(pulls, 1.0):
            continue
        if any(_hit(_pulls(U[np.ix_(rp, cp)]), 1.0)
               for rp in perms for cp in perms):
            permutation_promotions.add(key[:-1])

    # delta materially controls the solar angle.  Fit its exact a+b cosδ
    # dependence for the committed form and solve the 0.5σ marginal band.
    def committed_U(d):
        phi = r / math.sqrt(float(d10))
        U_e = _rot("23", -phi) @ _rot("13", float(h10), d)
        return U_e.conj().T @ U_TBM

    y0 = _observables(committed_U(0.0))["12"]
    ypi = _observables(committed_U(math.pi))["12"]
    solar_offset, solar_cosine = (y0 + ypi) / 2, (y0 - ypi) / 2
    assert abs(_observables(committed_U(delta))["12"]
               - (solar_offset + solar_cosine * math.cos(delta))) < 1e-12
    ref, sig = NUFIT_PMNS["sin2_12"]
    cos_bounds = sorted(((ref - 0.5 * sig - solar_offset) / solar_cosine,
                         (ref + 0.5 * sig - solar_offset) / solar_cosine))
    principal_delta_band = tuple(math.degrees(x) for x in
                                 (math.acos(cos_bounds[1]),
                                  math.acos(cos_bounds[0])))

    report("PMNS FORM SENSITIVITY (local audit around the committed ansatz)")
    report("=" * 70)
    report(f"  declared menu: {len(declared)} rows, "
           f"{len(unique_fingerprints)} unique angle fingerprints")
    report(f"  marginal box hits: {len(declared_hits[1.0])} at 1σ, "
           f"{len(declared_hits[0.5])} at 0.5σ")
    report(f"  distinct plane diagnostic: {len(distinct_planes)} rows, "
           f"{distinct_hits[1.0]} / {distinct_hits[0.5]} hits")
    report(f"  add atan/asin theta maps: {len(extended)} rows, "
           f"{extended_hits[1.0]} / {extended_hits[0.5]} hits")
    report("  committed form theta forks (max pull / diagonal score):")
    for name in ("h10", "lambda", "atan_lambda"):
        row = fork_comparison[name]
        report(f"    {name:<11s} {row['max_abs_pull']:.3f}σ / "
               f"{row['diagonal_score']:.5f}")
    report(f"  theta sign promotions: {len(theta_sign_promotions)}; "
           f"label permutation promotions: {len(permutation_promotions)}")
    report(f"  delta controls sin²theta12 as {solar_offset:.8f} "
           f"+ {solar_cosine:.8f} cos(delta)")
    report(f"  committed 0.5σ solar band: "
           f"{principal_delta_band[0]:.2f}°–{principal_delta_band[1]:.2f}°; "
           f"delta = {math.degrees(delta):.2f}°")

    return {
        "declared_n_rows": len(declared),
        "declared_unique_fingerprints": len(unique_fingerprints),
        "declared_n_within_1sigma": len(declared_hits[1.0]),
        "declared_n_within_half_sigma": len(declared_hits[0.5]),
        "committed": committed,
        "committed_pulls": committed_row[2],
        "distinct_plane_n_rows": len(distinct_planes),
        "distinct_plane_hits": distinct_hits,
        "extended_n_rows": len(extended),
        "extended_hits": extended_hits,
        "fork_comparison": fork_comparison,
        "theta_sign_promotions": sorted(theta_sign_promotions),
        "permutation_promotions": sorted(permutation_promotions),
        "delta_conjugation_invariant": True,
        "solar_offset": solar_offset,
        "solar_cosine": solar_cosine,
        "principal_delta_band_deg": principal_delta_band,
        "delta_deg": math.degrees(delta),
    }


if __name__ == "__main__":
    run()
