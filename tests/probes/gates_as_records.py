"""The six gates as statements about records.

THE QUESTION. The gates of the gravity lane are written in the words
of the theory they recover: energy, weight, cone, tensor. A closure
holds none of these. It holds a clock, the rate at which its phase
turns, and a count, the quanta it keeps. A record is a number two
closures can share: a ratio of two rates, a count, a coincidence of
two arrivals. What does each gate say as a record, and in which case
does the record fail?

PREMISES. P1 the lump of lump_in_empty_field.py, a clock of rate Omega
holding the count Q, in two kinds whose binding differs by a quarter
(two_lumps.py). P2 the exchange of the one list's induced kernel
between two lumps (two_lumps.py W2), V = Phi w, with Phi = -G E_S / r
the bias a source S sets up at the separation r and w the weight the
kernel couples to. P3 the saddle's symmetry (saddle_stabilizer.py), the one
list and the bridge's readings (interference/boundary.py,
horizon_amount.py), the pole's polarizations (quadrupole_radiation.py),
the labelled list's kernel (one_list_kernel.py) and the first
post-Newtonian order (post_newtonian.py).

STATEMENTS (the numbers below check each one).

  U (two clocks keep their ratio). A lump's rate is its weight per
     added quantum, Omega = dE/dQ, so the weight is the rate summed
     over the count. The exchange couples to that weight at every count.
     A bias Phi therefore moves the rate by the fraction Phi in both
     kinds: d ln Omega / d Phi = 1, and the ratio of the two clocks
     records no bias. FAILS for a bias that couples to the count instead
     of the weight (the count is the other amount a lump carries). Each
     rate then moves by Phi m, the fractions are m / Omega, 1.099 and
     1.396, and the ratio of the two clocks moves by 0.297 Phi. It fails between rows when the
     vacuum holds a frame vector (saddle_stabilizer.py Z8).
  C (two signals sent together arrive together). On the saddle every
     row has one symbol, so a scalar, a Weyl and a vector signal cover
     one path in one count of ticks. FAILS with a frame vector in the
     vacuum: four scalar rows then arrive apart by a fraction of the
     travel time of order one.
  S (clocks slow in each other's presence). Phi is negative: the rate
     of a clock falls as a second closure approaches, so closures that
     slow together stay together. On the one list the stiffness is
     positive at every contact. FAILS on the list without the amount:
     the stiffness is -4.5 at the covariant contact, the rate rises,
     and nothing gathers.
  M (the grain is a count of the boundary). The regulator in units of
     the boundary, t0 M_b^2 = Sigma / 2 pi, is a pure number, 4.11,
     6.98 and 9.84 at the three contacts, inside 0.6 to 15. FAILS on
     the list without the amount, where no regulator solves the list
     at the covariant contact, and on the labelled bridge at the
     boundary's mass (horizon_amount.py A5), which keeps b0 and fails
     S and M.
  T (a ring of clocks records two patterns). A bias passing through a
     ring of clocks is recorded in two transverse patterns that keep
     the ring's area.  A sphere of clocks that expands and contracts
     (a breathing mode) sends none. F0/F2 = -2 is not a choice: the kernel's q^2 term is local
     only there (composite_pole.py K9). FAILS for the labelled list,
     the bridge read with its searchlight label at the boundary's mass:
     solved freely, its kernel has a zero at 0.624 M_b, inside the
     boundary (one_list_kernel.py K3).
  P (a local ratio does not record the place). The bias moves rulers
     and clocks by one fraction, gamma = 1, so no ratio taken inside a
     laboratory records where the laboratory is.  A lump's stress
     moves the exchange by less than 1e-6.  The stiffness depends on the
     rows' counts, not their masses (horizon_amount.py A7).  FAILS for
     the amount read as massive below the boundary.  It passes S and M,
     but its stiffness depends on that mass at order one.  Whatever sets
     the mass then sets G, and a local ratio records the place.

RESULT. Each gate is a statement about ratios of rates, counts and
coincidences, with a computed case in which the record fails.
Universality is the statement that a closure's weight is its own rate
summed over its count.  The bias is then a common factor on every
clock, and a common factor is not a record.

Scope: the exchange at linear order between static spherical lumps,
free rows on the saddle, the pole between conserved slow sources.

Usage: python3 tests/probes/gates_as_records.py     (about one minute)
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r_  # noqa: E402
import horizon_amount as ha  # noqa: E402
import one_list_kernel as olk  # noqa: E402
import post_newtonian as pn  # noqa: E402
import quadrupole_radiation as qr  # noqa: E402
import saddle_stabilizer as ss  # noqa: E402
import two_lumps as tl  # noqa: E402

STEP = 0.005                                 # half width in mu/mu0 of the count's variation
SEP = 120.0                                  # separation, carrier lengths
SOURCE = 0.6                                 # mu/mu0 of the lump that lays the bias


def _silent(*args, **kwargs):
    pass


def clock(frac, src, sep=SEP):
    """A kind's rate, and the fraction by which a bias moves it per unit bias, for two readings of the weight."""
    lo, mid, hi = (tl.source((frac + s * STEP)) for s in (-1, 0, 1))
    omega = tl.lf.lift(tl.lf.lump(frac * tl.lf.MU0))["omega"]
    dQ = hi["mQ"] - lo["mQ"]                                          # the count is Q = m Q at m = 1
    phi = -src["E"] / sep                                             # the bias, in units of G
    exchange = lambda s: tl.potential(s, src, sep, weights=(1.0, 1.0)) / sep
    by_kernel = (exchange(hi) - exchange(lo)) / dQ / phi / omega
    return dict(omega=omega, rate_per_quantum=(hi["E"] - lo["E"]) / dQ, kernel=by_kernel, count=1.0 / omega,
                binding=mid["binding"])


def with_amount(n):
    """The one list with n real scalar units of amount in the bridge's place; n = None removes the row."""
    one = B.rows(r_.d10, r_.d11, r_.n7, r_.n26, bridge="amount" if n is not None else "absent")
    return [(name, n if name == "bridge" and n is not None else states, kind, share, reading)
            for name, states, kind, share, reading in one]


def grain(table):
    """The stiffness with no regulator at the covariant contact, the verdict of S, and t0 M_b^2 at each contact."""
    xi = float(r_.xi_bridge)
    roots = [B.anchored(table, c, xi) for c in B.CONTACTS]
    return dict(sigma=B.sigma(table, 0.0, 1.0, xi), inside=[B.inside(table, c, xi) for c in B.CONTACTS],
                t0=[root[0] / (2 * math.pi) if root else float("nan") for root in roots])


def run(report=print):
    out = {}

    # U: two clocks of different kinds in one bias.
    src = tl.source(SOURCE)
    a, b = clock(tl.FRACS[0], src), clock(tl.FRACS[1], src)
    for c in (a, b):
        assert abs(c["rate_per_quantum"] / c["omega"] - 1) < 1e-3
        assert abs(c["kernel"] - 1) < 1e-3                 # the count reading is m / Omega by its definition
    assert a["binding"] > 0 > b["binding"] + 0.2
    ratio_kernel, ratio_count = a["kernel"] - b["kernel"], a["count"] - b["count"]
    assert abs(ratio_kernel) < 1e-3 and abs(ratio_count) > 0.25
    saddle = ss.run(report=_silent)
    assert saddle["weight_over_energy"] < 1e-12 and saddle["frame_eotvos"] > 1e-2
    out.update(rates=(a["omega"], b["omega"]), rate_per_quantum=(a["rate_per_quantum"], b["rate_per_quantum"]),
               shift=(a["kernel"], b["kernel"]), shift_count=(a["count"], b["count"]),
               ratio_moves=ratio_kernel, ratio_moves_count=ratio_count, rows_frame=saddle["frame_eotvos"])
    report(f"U rates {a['omega']:.4f} and {b['omega']:.4f}, weight per added quantum {a['rate_per_quantum']:.4f} "
           f"and {b['rate_per_quantum']:.4f}; a bias moves the rates by the fractions {a['kernel']:.4f} and "
           f"{b['kernel']:.4f} of itself, the ratio of the clocks by {ratio_kernel:+.0e}; a bias reading the "
           f"count gives {a['count']:.3f} and {b['count']:.3f}, the ratio moving by {ratio_count:+.3f}. Between "
           f"rows a frame vector gives {saddle['frame_eotvos']:.2f}")

    # C: signals of different rows over one path.
    counts = saddle["invariant_counts"]
    apart = [abs(s / saddle["shared_speed"] - 1) for s in saddle["frame_speeds"]]
    assert counts["rank 2"] == 1 and counts["rank 1"] == 0 and saddle["weyl_invariants"] == 1
    assert saddle["vector_det"]["null"] < 1e-12 and saddle["weyl_vector_error"] < 1e-12 and min(apart) > 1e-3
    out.update(symbols=1, arrivals_apart=0.0, arrivals_apart_frame=apart)
    report("C one symbol for scalar, Weyl and vector rows, arrivals together. With a frame vector four scalar "
           "rows arrive apart by " + ", ".join(f"{x:.3f}" for x in apart) + " of the travel time")

    # S and M: the sign of the bias and the grain as a count.
    one, bare = (grain(with_amount(n)) for n in (r_.n7 * r_.n26, None))
    band = lambda g: all(0.6 < t < 15.0 for t in g["t0"])
    assert all(one["inside"]) and band(one) and one["sigma"] > 0 > bare["sigma"]
    assert not all(bare["inside"]) and not band(bare)
    labelled = ha.run(report=_silent)["readings"]["boundary"]     # (keeps b0, S, M) for the labelled bridge
    assert labelled == (True, False, False)
    out.update(sigma=one["sigma"], sigma_bare=bare["sigma"], grain=one["t0"], readings_boundary=labelled,
               bias_sign=-1 if one["sigma"] > 0 else 1, bias_sign_bare=-1 if bare["sigma"] > 0 else 1)
    report(f"S stiffness {one['sigma']:.2f}, the rate falls near a second closure. Without the amount "
           f"{bare['sigma']:.2f}, the rate rises")
    report("M t0 M_b^2 = " + ", ".join(f"{t:.2f}" for t in one["t0"]) + ". Without the amount no regulator at "
           "the covariant contact. The labelled bridge at the boundary keeps b0 and fails S and M")

    # T: the patterns a ring of clocks records.
    quad = qr.run(report=_silent)
    zero_one, zero_labelled = ha.zero_q(r_.n7 * r_.n26, 0.0), olk.run(report=_silent)["q_zero_labelled"]
    assert quad["tt_error"] < 1e-12 and abs(quad["breathing"]) < 1e-12
    assert zero_one is None and zero_labelled < 1.0
    out.update(patterns=2, breathing=abs(quad["breathing"]), zero_labelled=zero_labelled)
    report(f"T two transverse patterns to {quad['tt_error']:.0e}, a breathing sphere sends "
           f"{abs(quad['breathing']):.0e}. The labelled list's kernel has a zero at {zero_labelled:.3f} M_b")

    # P: rulers against clocks.
    post = pn.run(report=_silent)
    lumps = tl.run(report=_silent)
    assert post["gamma"] == 1 and post["beta"] == 1
    assert lumps["stress_share"] < 1e-6
    amount = ha.run(report=_silent)
    assert amount["mass_response_one"] == 0.0 and amount["mass_response_massive"] < -0.5
    out.update(gamma=1.0, stress_share=lumps["stress_share"], mass_response_one=amount["mass_response_one"],
               mass_response_massive=amount["mass_response_massive"])
    report(f"P rulers and clocks move by one fraction, gamma = 1, a lump's stress moving the exchange by "
           f"{lumps['stress_share']:.0e}. The amount read as massive responds to its mass by "
           f"{amount['mass_response_massive']:+.2f}, the one list by {amount['mass_response_one']:.0e}")

    report("RESULT: each gate is a statement about ratios of rates, counts and coincidences, with a computed "
           "case in which the record fails. A closure's weight is its rate summed over its count, so the bias "
           "is a common factor on every clock.")
    return out


if __name__ == "__main__":
    run()
