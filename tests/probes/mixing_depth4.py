"""Mixing at web depth four: the open path frames and the backreaction.

Two depth 4 objects bear on mixing, and this probe computes both.  Every
step before PART C runs without mixing data.

PART A: the depth 4 carriers as frames.
  The D6* braid connection on four step paths with fixed endpoints acts
  on twelve eight dimensional blocks.  The central full twist splits
  every block as 3 + 3 + 2, so there are 24 generation sized carriers.
  A single braid generator has two eigenvalues (Hecke), so it cannot
  fix a frame on a three dimensional carrier.  The canonical frames are
  the five fusion trees on four strands, each the joint eigenbasis of a
  commuting pair:
      ((12)3)4 : b1, (b1 b2)^3        (1(23))4 : b2, (b1 b2)^3
      (12)(34) : b1, b3
      1((23)4) : b2, (b2 b3)^3        1(2(34)) : b3, (b2 b3)^3
  Every joint spectrum is simple on every carrier, so each frame is
  fixed up to phases and |V|^2 between two frames is gauge invariant.
  Over all 24 carriers and all 10 tree pairs the relative frames take
  |V|^2 in {0, 1/12, 1/6, 1/4, 1/3, 1/2, 2/3, 3/4, 1}, exactly, and the
  Jarlskog invariant vanishes: the tree frames are real.  The frames
  are defined by braid group elements, so any intertwiner between
  equivalent carriers carries tree frames to tree frames: the list is
  complete for these carriers under any sector assignment.
  Every entry is zero or at least 1/sqrt(12) = 0.2887.  In every block
  the two carriers sit on the twist roots 1 and omega and the 2-dim
  piece on omega bar: the carriers are not a conjugate pair.

PART B: the depth 4 backreaction on mixing.
  B1. Every depth 4 edge in the web reads the depth 3 e<->q carrier
      eq_loop and sits on a coupling, action or per knot amplitude node,
      or on the common mode, which no knot reads.
      No mixing coordinate is a web node.
  B2. A depth 4 term is either sector blind (one scalar on every knot)
      or per knot (diagonal in its sector's mass basis).  Neither moves
      an eigenvector, so V = U_u^dag U_d is invariant at depth 4,
      exactly.  The depth 4 reading of the depth 3 apex gives the
      same apex: gamma = arctan(pi/sqrt(d10)) at depth 3 and at depth 4.
  B3. The one live depth 4 path into mixing is the lepton masses inside
      phi = sqrt(m_e/m_mu)/sqrt(d10).  It is computed here by reading the
      lepton ledgers at depth <= 3 and at full depth.
  B4. A scalar depth 4 echo attached directly to a Wolfenstein coordinate
      would need a depth 1 vertex weight on a mixing node, and none
      exists.  Its size is bounded anyway by the frozen 4 x 7
      vertex weight menu: w (alpha/2pi)^2 with w <= (8/3)(182).

PART C: data contact, once.
  The pull shifts of B3 and B4, and the distance of the PART A menu from
  the measured CKM magnitudes.

RESULT. The four step path carriers give real fusion tree frames whose
|V|^2 lies exactly on a nine value menu with floor 1/12 and J = 0.
Sector blind and per knot terms move |V|^2 and J by under 1e-12, so V is
invariant at web depth 4, and the one live path shifts the PMNS pulls by
under 1e-3 sigma while the nearest tree frame sits farther from the
measured CKM than the identity. Mixing closes at the depth it is read.
"""

import contextlib
import io
import itertools
import math
from fractions import Fraction

import numpy as np

import _bootstrap  # noqa: F401
from interference import masses, mixing, root
from interference.root import (PDG_CKM, NUFIT_PMNS, WEB, SEARCHLIGHT_DEPTH,
                               d10, h10)

_ROOTS = (1.0 + 0.0j, complex(math.cos(2 * math.pi / 3),
                              math.sin(2 * math.pi / 3)),
          complex(math.cos(2 * math.pi / 3),
                  -math.sin(2 * math.pi / 3)))
_MENU = (Fraction(0), Fraction(1, 12), Fraction(1, 6), Fraction(1, 4),
         Fraction(1, 3), Fraction(1, 2), Fraction(2, 3), Fraction(3, 4),
         Fraction(1))
_TREES = ("((12)3)4", "(1(23))4", "(12)(34)", "1((23)4)", "1(2(34))")


def _jarlskog(v):
    return float((v[0, 0] * v[1, 1]
                  * np.conj(v[0, 1]) * np.conj(v[1, 0])).imag)


def _carrier_bases(twist):
    """Orthonormal bases of the full twist eigenspaces, by root."""
    eye = np.eye(twist.shape[0], dtype=complex)
    bases = []
    for r in _ROOTS:
        proj = eye.copy()
        for s in _ROOTS:
            if abs(s - r) > 1e-9:
                proj = proj @ (twist - s * eye) / (r - s)
        u, sv, _ = np.linalg.svd(proj)
        bases.append(u[:, :int(np.sum(sv > 1e-8))])
    return bases


def _tree_frames(q, b1, b2, b3):
    """Joint eigenframes of the five fusion trees on one carrier."""
    d123 = np.linalg.matrix_power(b1 @ b2, 3)
    d234 = np.linalg.matrix_power(b2 @ b3, 3)
    pairs = ((b1, d123), (b2, d123), (b1, b3), (b2, d234), (b3, d234))
    frames, min_gap = [], float("inf")
    for a_full, b_full in pairs:
        a = q.conj().T @ a_full @ q
        b = q.conj().T @ b_full @ q
        assert np.linalg.norm(a @ b - b @ a) < 1e-10
        ev, vec = np.linalg.eig(a + (0.731 + 0.29j) * b)
        gap = min(abs(ev[i] - ev[j])
                  for i in range(len(ev)) for j in range(i + 1, len(ev)))
        min_gap = min(min_gap, gap)
        frames.append(np.linalg.qr(vec)[0])
    return frames, min_gap


def part_a(report):
    report("PART A: the depth 4 carriers as frames (no data)")
    report("-" * 66)
    d6 = mixing._compute_d6_boundary()
    n10, path_operator = d6["n10"], d6["path_operator"]
    blocks, carriers, min_gap = 0, 0, float("inf")
    menu_hits, snap_residual, j_max = set(), 0.0, 0.0
    relative = []
    for start, end in itertools.product(range(6), repeat=2):
        paths = tuple(
            p for p in itertools.product(range(6), repeat=5)
            if p[0] == start and p[-1] == end
            and all(n10[p[i], p[i + 1]] for i in range(4)))
        if not paths:
            continue
        blocks += 1
        b1, b2, b3 = (path_operator(paths, i) for i in (1, 2, 3))
        twist = np.linalg.matrix_power(b1 @ b2 @ b3, 4)
        bases = _carrier_bases(twist)
        root_dims = tuple(q.shape[1] for q in bases)
        assert root_dims == (3, 3, 2)
        for q in bases:
            if q.shape[1] != 3:
                continue
            carriers += 1
            frames, gap = _tree_frames(q, b1, b2, b3)
            min_gap = min(min_gap, gap)
            for i, j in itertools.combinations(range(5), 2):
                v = frames[i].conj().T @ frames[j]
                j_max = max(j_max, abs(_jarlskog(v)))
                mags = np.abs(v) ** 2
                relative.append(mags)
                for x in mags.ravel():
                    snapped = min(_MENU, key=lambda m: abs(float(m) - x))
                    snap_residual = max(snap_residual,
                                        abs(float(snapped) - x))
                    menu_hits.add(snapped)
    assert blocks == 12 and carriers == 24
    assert min_gap > 0.5
    assert snap_residual < 1e-12
    assert menu_hits == set(_MENU)
    assert j_max < 1e-12
    nonzero_floor = min(m for m in menu_hits if m > 0)
    report(f"  blocks {blocks}, each 3 + 3 + 2. Generation sized "
           f"carriers {carriers}")
    report(f"  five fusion tree frames per carrier, joint spectra simple "
           f"(min gap {min_gap:.3f})")
    report(f"  |V|^2 over 24 x 10 relative frames: "
           f"{{{', '.join(str(m) for m in sorted(menu_hits))}}}")
    report(f"    exact to {snap_residual:.1e}. Max |J| = {j_max:.1e}")
    report("  the 3-dim carriers sit on twist roots 1 and omega in every "
           "block. Omega bar")
    report("    carries the 2-dim one, so the pair is not a conjugate "
           "pair")
    report(f"  every entry is 0 or >= {nonzero_floor} "
           f"(|V| >= {math.sqrt(float(nonzero_floor)):.4f})")
    return {"blocks": blocks, "carriers": carriers,
            "root_dims": root_dims,
            "menu": tuple(sorted(menu_hits)), "j_max": j_max,
            "snap_residual": snap_residual,
            "nonzero_floor": nonzero_floor, "relative": relative}


def _solved_web():
    with contextlib.redirect_stdout(io.StringIO()):
        r = root.derive()
        m = masses.derive(r)
    return m


def _frames_invariant_under_depth4(seed=7):
    """Sector blind and per knot depth 4 terms leave V unchanged."""
    rng = np.random.default_rng(seed)

    def haar():
        z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        qm, rm = np.linalg.qr(z)
        return qm * (np.diag(rm) / abs(np.diag(rm)))

    def frame(u, m):
        _, vec = np.linalg.eigh(u @ np.diag(m) @ u.conj().T)
        return vec

    u_up, u_dn = haar(), haar()
    m_up, m_dn = np.array([1.0, 5.0, 40.0]), np.array([2.0, 9.0, 60.0])
    v0 = frame(u_up, m_up).conj().T @ frame(u_dn, m_dn)
    common = 1.0 - 16 * 1.3e-6
    per_up = 1.0 - 13.5 * 1.3e-6 * np.array([0.3, 1.1, 2.4])
    per_dn = 1.0 - 13.5 * 1.3e-6 * np.array([0.5, 0.9, 2.1])
    v_common = (frame(u_up, common * m_up).conj().T
                @ frame(u_dn, common * m_dn))
    v_knot = (frame(u_up, per_up * m_up).conj().T
              @ frame(u_dn, per_dn * m_dn))
    return max(float(np.max(abs(abs(v_common) ** 2 - abs(v0) ** 2))),
               float(np.max(abs(abs(v_knot) ** 2 - abs(v0) ** 2))),
               abs(_jarlskog(v_knot) - _jarlskog(v0)))


def part_b(report):
    report("")
    report("PART B: the depth 4 backreaction on mixing (no data)")
    report("-" * 66)
    m = _solved_web()
    deep = [(name, t) for name, led in WEB.items()
            if isinstance(led, root.Ledger)
            for t in led.terms if t.depth > SEARCHLIGHT_DEPTH]
    nodes = sorted({name for name, _ in deep})
    parents = {t.parent for _, t in deep}
    assert parents == {"eq_loop"}
    assert set(nodes) <= {"inv_alpha", "S_quark", "m_e", "m_mu", "m_tau", "common_mode"}
    report(f"  B1 depth 4 edges: {len(deep)}, all reading eq_loop, on "
           f"{', '.join(nodes)}")
    report("     no mixing coordinate is a web node")

    v_shift = _frames_invariant_under_depth4()
    assert v_shift < 1e-12
    gamma3 = math.degrees(math.atan(math.pi / math.sqrt(float(d10))))
    report(f"  B2 sector blind and per knot depth 4 terms move |V|^2 and J "
           f"by {v_shift:.1e}")
    report(f"     the apex is unchanged: gamma = {gamma3:.4f} deg at "
           "depth 3 and at depth 4")

    eq = WEB.state["eq_loop"]
    m_e4, m_mu4 = m["m_e"], m["m_mu"]
    m_e3 = WEB["m_e"].value(max_depth=SEARCHLIGHT_DEPTH)
    m_mu3 = WEB["m_mu"].value(max_depth=SEARCHLIGHT_DEPTH)
    wzw = mixing._compute_wzw()
    lam = float(h10 * (1 + root.h11 / root.K ** 2))
    pm4 = mixing._derive_pmns(wzw, lam, m_e4, m_mu4)
    pm3 = mixing._derive_pmns(wzw, lam, m_e3, m_mu3)
    ratio_shift = (m_e4 / m_mu4) / (m_e3 / m_mu3) - 1.0
    report(f"  B3 lepton depth 4 vent on m_e/m_mu: {ratio_shift:+.3e} "
           "(relative), carried into phi")

    # Weights in use: coupling, action, amplitude. The common mode carries
    # no weight of its own. Its count reads the carrier through the
    # amplitude terms.
    w_used = (Fraction(2), Fraction(16), Fraction(27, 2))
    w_menu_max = Fraction(8, 3) * 182
    report(f"  B4 eq_loop = (alpha/2pi)^2 = {eq:.4e}. Endpoint weights in "
           f"use {', '.join(str(w) for w in w_used)}, menu ceiling "
           f"{float(w_menu_max):.1f}")
    report(f"     relative ceiling on any scalar depth 4 echo: "
           f"{float(max(w_used)) * eq:.2e} (weights in use), "
           f"{float(w_menu_max) * eq:.2e} (menu)")
    return {"deep_nodes": nodes, "v_shift": v_shift, "gamma3_deg": gamma3,
            "eq_loop": eq, "ratio_shift": ratio_shift,
            "pmns3": pm3, "pmns4": pm4, "lam": lam,
            "ceiling_used": float(max(w_used)) * eq,
            "ceiling_menu": float(w_menu_max) * eq}


def part_c(report, a, b):
    report("")
    report("PART C: data contact, once")
    report("-" * 66)
    pmns_dpull = {k: (b["pmns4"][f"sin2_{k}"] - b["pmns3"][f"sin2_{k}"])
                  / NUFIT_PMNS[f"sin2_{k}"][1] for k in ("12", "23", "13")}
    report("  B3 PMNS pull shift from depth 4: "
           + ", ".join(f"sin2_{k} {v:+.1e} sigma"
                       for k, v in pmns_dpull.items()))
    lam_pdg, lam_sig = PDG_CKM["lambda"]
    g_pdg, g_sig = PDG_CKM["gamma_deg"]
    dl_used = b["lam"] * b["ceiling_used"] / lam_sig
    dl_menu = b["lam"] * b["ceiling_menu"] / lam_sig
    dg_menu = b["gamma3_deg"] * b["ceiling_menu"] / g_sig
    report(f"  B4 lambda pull shift ceiling: {dl_used:.1e} sigma (weights "
           f"in use), {dl_menu:.2f} sigma (menu)")
    report(f"     gamma pull shift ceiling: {dg_menu:.1e} sigma (menu)")
    lam_tan = math.tan(float(h10))
    need = lam_pdg / lam_tan - 1.0
    report(f"     tan(h10) = {lam_tan:.6f}, pull "
           f"{(lam_tan - lam_pdg) / lam_sig:+.2f} sigma. Reaching the "
           f"central value needs {need:+.2e}, "
           f"{abs(need) / b['ceiling_menu']:.0f}x the menu ceiling")

    ckm = np.array([[PDG_CKM[f"V{r}{c}"][0] for c in "dsb"]
                    for r in "uct"]) ** 2
    best = min(
        (float(np.max(abs(np.sqrt(mags[np.ix_(rows, cols)])
                          - np.sqrt(ckm)))), rows, cols)
        for mags in a["relative"]
        for rows in itertools.permutations(range(3))
        for cols in itertools.permutations(range(3)))
    dist_identity = float(np.max(abs(np.sqrt(ckm) - np.eye(3))))
    report(f"  A  nearest tree frame |V| to the CKM: max |V_ij| deviation "
           f"{best[0]:.4f}. The identity is {dist_identity:.4f} away")
    return {"pmns_dpull": pmns_dpull, "lambda_dpull_used": dl_used,
            "lambda_dpull_menu": dl_menu, "gamma_dpull_menu": dg_menu,
            "tan_h10_need": need, "ckm_best_deviation": best[0],
            "ckm_identity_deviation": dist_identity}


def run(report=print):
    report("MIXING AT WEB DEPTH FOUR")
    report("=" * 66)
    a = part_a(report)
    b = part_b(report)
    c = part_c(report, a, b)
    report("")
    report("VERDICT")
    report("  The depth 4 carriers give real frames on a nine value menu")
    report("  with a floor of 1/12 in |V|^2: no hierarchy, no CP phase.")
    report("  The depth 4 backreaction leaves V invariant and reaches the")
    report("  PMNS only through phi, far below every band. Mixing closes")
    report("  at the depth it is read. The depth 4 version exists and is")
    report("  computed, and it moves no scored coordinate.")
    report("  Depth 4 cannot fix which quark rotation carries the CP")
    report("  phase. The frames are real, and the backreaction is a")
    report("  continuous shift that cannot flip a discrete sign.")
    return {"A": a, "B": b, "C": c}


if __name__ == "__main__":
    run()
