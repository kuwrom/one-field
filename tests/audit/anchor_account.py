"""
anchor_account.py: the account of the anchor relation, with every factor between m_e and G entered once.

The lepton masses are web ledgers with registered echo terms
(interference/masses.py).  The anchor inverts the electron's closure in
one expression (root.py, _F_e_of).  g_path_trials.py,
gravity_resolution.py and g_lane_exhaustion.py each write that
expression out again.  The comparison with G passes through a constant
printed to six digits.  This module holds the relation as one object
and checks what a hand written product cannot.

THE ACCOUNT.  With e = alpha/(2 pi), Delta_k = 1 + (B/A) cos(theta + 2 pi k/3):
    m_k  = M_b C_EFF exp(-S0) Delta_k^2 (1 + c2 e) (1 - w4 e^2 Delta_k)
    v_EW = M_b exp(-S0 + C2 - dS + w16 e^2)
    G    = 1/(c M_b^2)
    1/alpha(0): the root of x^3 = 16 A [(1 - 1/(2 pi)) x^2 - 1/(2 pi^2)],
                A = 24 pi - (b0/2 pi)(S0 - C2), the lever arm's endpoint
The last line is root.py's own.  Running from the boundary over
S0 - C2 gives 1/alpha_G2(mu*) = 32/pi, "Computed here, not assumed".
The EM base is that endpoint divided by 16.
Each symbol on the right is an Entry: its value, the web depth it
enters at, the repository's status word for it, and where it is fixed.

WHAT THE OBJECT ENFORCES (the numbers below check each one).
  A1 (nothing left out).  Evaluated at the entries' values, the
     account gives the engine's M_Pl, m_mu, m_tau and v_EW to 1e-12, and
     the three hand written copies of the inversion agree with it.
  A2 (readers are computed).  An entry's readers are the observables
     that move when it is shifted.  They are found by finite
     differences on the account, never declared.
  A3 (what fixes each entry).  For each entry, the largest shift that
     keeps every reader other than G within two standard deviations.
  A4 (where G's residual can sit).  For each entry, the shift that
     takes G's pull to zero and the largest pull it then leaves on the
     other observables.  An entry is a home if that pull stays below
     two.  There is one home: c, the boundary reading, which only G
     reads.  S0, the exponent every mass shares, moves no mass ratio.
     But S0 is also the length of the lever arm whose endpoint is the
     EM base.  The shift that zeroes G's pull moves 1/alpha(0) by 6e-4,
     which is twenty three thousand standard deviations.  Every entry
     specific to one sector is excluded by G_F or by the lepton ratios,
     by tens to thousands of standard deviations.
  A5 (second order entries enter once).  The account holds three
     entries of order e^2: w4 on each lepton, w16 on v_EW, and the square
     of the depth 2 term, which is absent.  G_F reads the electron's
     terms and v_EW's terms together.  It fixes w16 minus the electron's
     terms to 0.38 at two standard deviations, not each one alone.  A
     term common to every mass node is a factor on M_b, which is a
     shift of c.
  A6 (comparison constants from primaries).  sqrt(hbar c^5/G) from the
     exact SI constants and G = 6.67430e-11 is 1.2208901e19 GeV. The
     engine's comparison constant is printed as 1.22089e19, lower by
     0.105 ppm.  Computed from primaries the engine's G is
     6.6739371e-11 and its pull is -2.42.  Through the printed constant
     they are 6.673936e-11 and -2.43 (root.py G_PRED).  The anchor's own
     uncertainty is 6e-10 from m_e.

RESULT.  The 22 ppm has one possible home in the account: the boundary
reading c.  G is its only reader.  No other observable can test a term
placed there, no term specific to one sector can carry it, and the
fine structure constant holds the shared exponent fixed.  The account
carries the common mode as the engine reads it, G M_b^2 = 3 Q_l - 1
(root.py, THE COMMON MODE).  This is the count of the lepton circulant
on the depth four web, 1 - 7.50 (alpha/2pi)^2, fixed with no measured G
(tests/probes/common_mode_count.py).  It is computed from the account's
own theta, B/A and w4, so each of these is read by G as well as by the
lepton ratios.  The count carries none of the residual and moves G by
-10.1 ppm.  The entry c is the boundary reading beside the count.  Its
value is 1, and only G reads it.
If c is read as the whole second moment, c = 2/B^2
(tests/probes/boundary_whole.py B5), it moves with the harmonic, and
the lepton ratio holds it to 9e-10.  The residual is G's distance from
the recommended value, which depends on the rule that combines the
measurements (tests/audit/boundary_reading_scan.py).

Usage: python3 tests/audit/anchor_account.py     (about 20 s, runs the chain)
"""

import contextlib
import io
import math
import os
import sys
from dataclasses import dataclass

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import masses  # noqa: E402
from interference import root as r  # noqa: E402
import g_measurements as gm  # noqa: E402
import g_path_trials as gpt  # noqa: E402
import gravity_resolution as gr  # noqa: E402

HBAR, C_LIGHT, E_CHARGE = 1.054571817e-34, 299792458.0, 1.602176634e-19   # exact SI
OBSERVABLES = ("G", "v_EW", "m_mu", "m_tau", "1/alpha")
REL_M_E = 3.1e-10                                # CODATA relative uncertainty of m_e


@dataclass(frozen=True)
class Entry:
    name: str
    value: float
    depth: int          # web depth at which the entry enters
    status: str         # the repository's word for it
    source: str         # where it is fixed
    scope: str          # the nodes that read it in the account
    step: float         # finite difference step, in the entry's own units


class Account:
    """The anchor relation as a list of entries and the observables computed from them."""

    def __init__(self):
        a = 1.0 / r.inv_alpha_phys
        self.e = a / (2 * math.pi)
        self.m_e = r.M_E_ANCHOR_MEV
        rows = [
            Entry("S0", r.S_lepton, 1, "closure action, one loop (identification G tests)", "root.py S_lepton",
                  "every node", 1e-7),
            Entry("c", 1.0, 0, "identification by necessity, the boundary reading beside the count",
                  "root.py, the boundary reading", "G", 1e-7),
            Entry("C_EFF", r.C_EFF, 1, "the closure's common mode read whole (frame_coefficients F1), tested by G_F", "root.py C_EFF", "leptons", 1e-8),
            Entry("theta", float(r.h10), 1, "derived (h10)", "root.py h10", "leptons", 1e-9),
            Entry("B/A", r.BA_ratio, 1, "derived (sqrt d10)", "root.py BA_ratio", "leptons", 1e-9),
            Entry("c2", -1.0, 2, "the bridge vertex -2 pi R(h_bridge) (frame_coefficients F2), tested by G_F", "root.py WEB['lepton_EM']", "leptons", 1e-5),
            Entry("w4", r.d11**3 / r.d10, 4, "C_EFF x dim J3(O) (frame_coefficients F3), parent eq_loop", "masses.py lepton ledgers", "leptons", 1e-3),
            Entry("C2", float(r.C2_26), 1, "FORCED, Casimir vent", "root.py WEB['S_quark']", "v_EW", 1e-7),
            Entry("dS", float(r.delta_S), 1, "FORCED, E8 current trace", "root.py delta_S", "v_EW", 1e-7),
            Entry("w16", float(r.charge_trace * r.C2_26), 4, "FORCED, parent eq_loop", "root.py WEB['S_quark']",
                  "v_EW", 1e-3),
        ]
        names = [row.name for row in rows]
        assert len(set(names)) == len(names), "an entry is entered twice"
        self.entries = {row.name: row for row in rows}

    def values(self, **shift):
        unknown = set(shift) - set(self.entries)
        assert not unknown, f"not in the account: {sorted(unknown)}"
        return {k: row.value + shift.get(k, 0.0) for k, row in self.entries.items()}

    def observe(self, **shift):
        """G (in units of 1/M_Pl(CODATA)^2), v_EW in GeV, m_mu and m_tau in MeV, from the account."""
        x, e = self.values(**shift), self.e
        gap = lambda k: 1.0 + x["B/A"] * math.cos(x["theta"] + 2 * math.pi * k / 3)
        closure = lambda k: (x["C_EFF"] * math.exp(-x["S0"]) * gap(k)**2 * (1 + x["c2"] * e)
                             * (1 - x["w4"] * e * e * gap(k)))
        M_b = self.m_e / closure(1)                                    # MeV: the anchor inversion
        amp = [math.sqrt(closure(k)) for k in range(3)]                # the common mode's count, P v^2 P
        q0 = (sum(amp) / 3) ** 2
        count = (sum(closure(k) for k in range(3)) / 3 - q0) / q0
        lever = 1.0 / r.alpha_G2_Pl - float(r.b0_G2) / (2 * math.pi) * (x["S0"] - x["C2"])
        inv_alpha = 16.0 * lever * (1.0 - 1.0 / (2 * math.pi))        # the cubic's root, by iteration
        for _ in range(60):
            inv_alpha = 16.0 * lever * ((1.0 - 1.0 / (2 * math.pi)) - 1.0 / (2 * math.pi**2 * inv_alpha**2))
        return {"M_b": M_b, "1/alpha": inv_alpha,
                "G": r.G_CODATA * (r.M_PL_CODATA_GEV * 1e3)**2 * count / (x["c"] * M_b**2),
                "v_EW": M_b / 1e3 * math.exp(-x["S0"] + x["C2"] - x["dS"] + x["w16"] * e * e),
                "m_mu": M_b * closure(2), "m_tau": M_b * closure(0)}

    def pulls(self, **shift):
        obs, ref = self.observe(**shift), dict(gm.references(), **{"1/alpha": gm.INV_ALPHA[0][1:]})
        return {k: (obs[k] - ref[k][0]) / ref[k][1] for k in OBSERVABLES}

    def jacobian(self):
        """d(pull)/d(entry) for every observable and entry, by central differences."""
        out = {}
        for name, row in self.entries.items():
            up, down = self.pulls(**{name: row.step}), self.pulls(**{name: -row.step})
            out[name] = {k: (up[k] - down[k]) / (2 * row.step) for k in OBSERVABLES}
        return out


def engine():
    """The chain's own values of the four observables' inputs."""
    with contextlib.redirect_stdout(io.StringIO()):
        m = masses.derive(r.derive())
    return {"M_b": r.M_Pl_MeV, "m_mu": m["m_mu"], "m_tau": m["m_tau"], "v_EW": m["v_EW_pred_GeV"],
            "1/alpha": r.inv_alpha_phys}


def run(report=print):
    out = {}
    acc = Account()

    # A1: nothing left out, and the hand written copies agree.
    eng, obs = engine(), acc.observe()
    gaps = {k: abs(obs[k] / eng[k] - 1.0) for k in eng}
    copies = {"g_path_trials.m_pl": gpt.m_pl() * 1e3, "gravity_resolution.m_pl": gr.m_pl() * 1e3}
    copy_gap = max(abs(v / obs["M_b"] - 1.0) for v in copies.values())
    assert max(gaps.values()) < 1e-12 and copy_gap < 1e-12
    out.update(completeness=max(gaps.values()), copies=copy_gap, entries=len(acc.entries))
    report(f"A1 nothing left out: {len(acc.entries)} entries give the engine's M_Pl, m_mu, m_tau, v_EW to "
           f"{max(gaps.values()):.0e}. The hand written inversions agree to {copy_gap:.0e}")

    # A2, A3: readers and tolerances.
    J, base = acc.jacobian(), acc.pulls()
    readers = {n: [k for k in OBSERVABLES if abs(J[n][k] * acc.entries[n].step) > 1e-6] for n in acc.entries}
    tolerance = {n: min((2.0 / abs(J[n][k]) for k in readers[n] if k != "G"), default=math.inf)
                 for n in acc.entries}
    alone = sorted(n for n in acc.entries if readers[n] == ["G"])
    assert readers["c"] == ["G"] and readers["S0"] == ["G", "1/alpha"]
    assert abs(base["G"] + 2.43) < 0.01 and max(abs(base[k]) for k in OBSERVABLES if k != "G") < 0.3
    out.update(readers=readers, tolerance=tolerance, read_by_G_alone=alone, pulls=base)
    report("A2 readers: " + "; ".join(f"{n} <- {', '.join(readers[n]) or 'none'}" for n in acc.entries))
    report("A3 two sigma tolerance without G: " + "; ".join(
        f"{n} {tolerance[n]:.1e}" for n in acc.entries if math.isfinite(tolerance[n]))
        + ". Unbounded: " + ", ".join(n for n in acc.entries if not math.isfinite(tolerance[n])))

    # A4: where G's residual can sit.
    homes, cost = [], {}
    for n in acc.entries:
        if abs(J[n]["G"]) < 1e-12:
            cost[n] = math.inf
            continue
        shift = -base["G"] / J[n]["G"]
        after = acc.pulls(**{n: shift})
        assert abs(after["G"]) < 0.05, n
        cost[n] = (shift, max(abs(after[k]) for k in OBSERVABLES if k != "G"))
        if cost[n][1] < 2.0:
            homes.append(n)
    excluded = {n: v for n, v in cost.items() if v != math.inf and n not in homes}
    assert homes == ["c"] and min(v[1] for v in excluded.values()) > 20 and cost["S0"][1] > 1e4
    out.update(homes=sorted(homes), cost=cost)
    report("A4 homes of G's residual: " + ", ".join(
        f"{n} (shift {cost[n][0]:+.2e}, other pulls below {cost[n][1]:.2f})" for n in sorted(homes))
        + ". Excluded: " + ", ".join(f"{n} ({v[1]:.0f} sigma)" for n, v in excluded.items()))

    # A5: the second order entries, and what G_F reads of them.
    e2 = acc.e**2
    together = acc.pulls(w16=1.0)["v_EW"] - base["v_EW"]              # one unit of w16 alone, in G_F sigmas
    k_common = -cost["c"][0] / (2 * e2)                               # the term, in units of e^2 on each mass
    common = acc.pulls(c=-2 * k_common * e2)                          # a term common to every mass is a factor on M_b
    second_order = sorted(n for n, row in acc.entries.items() if row.depth == 4)
    assert second_order == ["w16", "w4"] and abs(together) > 4 and abs(common["G"]) < 0.05
    assert max(abs(common[k] - base[k]) for k in OBSERVABLES if k != "G") < 1e-6
    w16_band = 2.0 / abs(together)
    out.update(second_order=second_order, w16_band=w16_band, common_shift=common["G"], k_common=k_common)
    report(f"A5 second order: entries {', '.join(second_order)}; G_F holds w16 to {w16_band:.2f} at two sigma; a term "
           f"common to every mass is a factor on M_b, a shift of c, and moves G alone ({common['G']:+.2f} at "
           f"{k_common:.2f} e^2)")

    # A6: comparison constants from primaries.
    exact = math.sqrt(HBAR * C_LIGHT**5 / r.G_CODATA) / E_CHARGE / 1e9
    rounding = 1e6 * (r.M_PL_CODATA_GEV / exact - 1.0)
    G_direct = r.COMMON_MODE * HBAR * C_LIGHT**5 / (r.M_Pl_GeV * 1e9 * E_CHARGE)**2
    pull_direct = (G_direct / r.G_CODATA - 1.0) / r.G_REL_SIGMA
    assert abs(rounding + 0.105) < 0.002 and abs(G_direct - 6.6739371e-11) < 1e-17
    assert abs(pull_direct + 2.419) < 0.002 and abs(r.G_PRED - 6.673936e-11) < 1e-16
    out.update(planck_exact=exact, rounding_ppm=rounding, G_direct=G_direct, pull_direct=pull_direct,
               anchor_uncertainty=2 * REL_M_E)
    report(f"A6 primaries: sqrt(hbar c^5/G) = {exact:.7e} GeV. The printed constant is {rounding:+.3f} ppm from it; "
           f"from primaries G = {G_direct:.7e} ({pull_direct:+.2f}), through the constant {r.G_PRED:.6e} "
           f"({base['G']:+.2f}). The anchor's own uncertainty is {2 * REL_M_E:.0e}")

    report("RESULT: G's residual has one home in the account, " + homes[0] + ", read by G alone. The shared "
           "exponent is held by 1/alpha(0) through the lever arm, and every entry particular to a sector is "
           "excluded by G_F or the lepton ratios.")
    return out


if __name__ == "__main__":
    run()
