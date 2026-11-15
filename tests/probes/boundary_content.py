"""The labelled bridge on one list: the rows that both the running and the stiffness count.

THE QUESTION. The transmutation of root.py runs alpha_G2 = 1/(24 pi)
down from the boundary M_b with b0 = 32/3. It counts the G2 gauge field
and d11 Dirac 7s. The stiffness entered in interference/gravity.py
counts 24 + 12 vector states, 90 fermion states, 7 G2 scalars and the
182 bridge channels, with every row massless at the regulator. With its
searchlight label, the bridge is n26 real 7s of the same G2 and the
scalars are one real 7, and the running counts neither. What does one
list give if both the running and the stiffness use it with that label?
The engine's one list instead reads the bridge as amount beyond the
searchlight's horizon, with no share of the running
(interference/boundary.py, tests/probes/horizon_amount.py). This probe
computes the labelled reading next to it.

PREMISES. P1 the boundary identification of root.py: M_b, the scale at
which alpha_G2 = 1/(24 pi), is the full Planck mass. With
1/G = Sigma/(2 pi t0) (sakharov_normalisation.py FACT 1), the regulator
is then t0 = Sigma/(2 pi) M_b^-2. P2 one loop running, with each row
leaving it at its mass. The anchor's b0 already takes this step. P3 the
per state weights of ledger_kernel.py R2 (scalar 1/6 - xi, Weyl 1/12,
vector 1/6 - c/2 for vector contact c, c = +1 the covariant value) and
the ledger's own sharp proper time cutoff, under which a row of mass m
enters with E2(m^2 t0) (composite_pole.py K5, one_list_kernel.py K1).
P4 R6 of regulator_free_ratios.py: a composite bound below the
regulator is not a row (the Higgs, and G2's breaking condensate at
v_break = v_EW/sqrt 2, couplings.py). The knots stay rows. A knot is
not a lump, since no carrier mass makes a lump a point particle
(closed_event_scale.py S5). P5 the regulator lies inside the boundary,
t0 >= M_b^-2. The list below M_b is the stiffness's list only up to
M_b, and every placement the repository adopts is inside: the
diagnostic's t0 = d11 l_Pl^2 sits at 0.577 M_b, and the
identification's Sigma/(2 pi) l_Pl^2 sits at 0.482 M_b on the
gravity.py ledger.

STATEMENTS (the numbers below check each one).

  B1 (the label). The G2 index of the E8 adjoint is h(G2) + n26 T(7)
     = 4 + 26 = 30 = h(E8). With this label, the bridge is n26 real 7s
     of the running G2. The 7 scalars are G2's breaking, the
     condensate at v_break = 1.4e-17 M_b, and by R6 they are not a
     row. B2 treats them as elementary light fields. The gravity.py
     ledger counts both in the stiffness and neither in the running.
     For the bridge, that is the amount reading of boundary.py. For the
     scalars, it is a second list.
  B2 (light in both). If they are also counted as light in the
     running, the bridge gives b0 = 19/3, the scalars 21/2, and both
     together 37/6. Running down to Lambda_G2, M_Pl then moves by
     1.6e13, 2.02 and 1.2e14. Through root.py's S_baryo =
     (h(G2)/h(F4)) S_lepton, eta_B moves by 1.4e-6, 0.73 and 5.6e-7.
     If they leave the running at mu* = M_b exp(-(S0 - C2(26))), M_Pl
     moves by 6.0e6, 1.82 and 1.1e7. The G2 endpoint alpha(mu*) =
     2 pi/(b0 C2(26)) moves away from pi/32. The anchor holds only if
     neither is light.
  B3 (absent from both). The light rows give Sigma = 13.5 - 18 c:
     -4.5 at the covariant contact, and 13.5 and 31.5 at c = 0 and -1.
     Attraction needs c < 3/4. The carrier's own six components (+1 at
     most) and the composite metric's own loop (-8/3) do not change
     the sign. A list of the knots alone (7.5) drops the vector rows.
     But b0's 44/3 is the G2 gauge field's own loop, and an elementary
     vector's stiffness is the heat trace row (determinant_vectors.py
     V5). So on one list the vectors count.
  B4 (the regulator inside the boundary). Given P1, t0^(-1/2) lies at
     or below the boundary exactly when Sigma >= 2 pi. That is the
     condition P5 imposes. The 27 of the gravity.py ledger, with the
     bridge read as amount and minimally coupled, puts the regulator
     at 0.48 M_b.
  B5 (the labelled bridge at the boundary). Suppose a mass operator
     exists (the repository does not supply one). A bridge of mass
     m_B >= M_b then leaves the running. If it sat below M_b by more
     than 5.5e-5 in the log, it would move M_Pl past the anchor's
     22 ppm. It enters the stiffness with its heavy row weight. Under
     the sharp cutoff, with the regulator inside the boundary, that
     weight is at most 182 (1/6 - xi) E2(1) = 4.33. So Sigma <= -0.175
     at c = +1, which is 6.46 short of the 2 pi of B4 (9.12 with the
     metric loop). An attractive stiffness whose anchored regulator
     lies inside the boundary needs a boundary row at 0.370 of a
     massless row's weight (0.462 with the loop). The sign alone needs
     0.154. E2 is convex, so among regulators with the same massless
     normalization the sharp cutoff gives a heavy row its smallest
     weight. The regulator exp(-t0/t) gives 2 K1(2) = 0.280, so
     Sigma = +3.65 at t0 = M_b^-2. That is attractive, but 2.64 short
     of the 2 pi that would put its anchored regulator there. Solved
     freely under this regulator, the regulator sits at 1.11 M_b. A
     broader shape can reach 0.370. The coefficient also depends on
     the regulator's field content. Under Pauli-Villars the same
     composite metric gives
     M_Pl^2 = N M^2/(192 pi^2) with three regulator fields and
     N M^2/(288 pi^2) with four (Li and Vaman arXiv:2512.06132,
     Eq. 2.41 and Sec. II). Solved freely under the sharp cutoff,
     Sigma = 3.70 and t0 = 0.589 l_Pl^2. This puts the regulator at
     1.30 M_b, outside the boundary, with 6 pi/Sigma = 5.10 and the
     spin-2 zero at 0.624 M_b (one_list_kernel.py K3). At c = 0,
     Sigma = 14.27 and t0 = 2.27 l_Pl^2, inside.
  B6 (the boundary coupling is not induced). Treat the coupling as a
     compositeness condition with the stiffness's cutoff (Akama,
     Terazawa). The matter rows then induce 1/alpha_G2(mu) =
     (b_m/2 pi) ln(Lambda/mu), with b_m = (4/3) d11 T(7) = 4. Reaching
     1/(24 pi) at the boundary needs ln(Lambda/M_b) = 12 pi^2 = 118.4
     (55.7 with the gravity.py list, b_m = 17/2). That is a cutoff far
     above the boundary, while t0^(-1/2) = 0.48 M_b. The engine's
     coupling at M_b is a boundary datum, not an induced one.
  B7 (what the labelled bridge's sign needs). Under the ledger's own
     cutoff, the 6.46 of B5 equals 39 light G2-singlet real scalars
     beyond the list, or 78 Weyl states. The G2-singlet content that
     the repository names adds 2.34, each item at its most favourable
     weight: the carrier's six components massless (+1), the 26's two
     F4 singlets taken as light G2 singlets rather than bridge
     channels (+1/3), the G2 scalars at the boundary (+0.17), the
     bridge at the anchor's mass edge (+7e-4), the G2-face radial mode
     light (+1/6) and an elementary Higgs doublet light (+2/3, which
     R6 excludes). Sigma at the boundary is then +2.17, still 4.12
     short. The two F4 singlets alone move the anchored regulator to
     t0 = 0.606 l_Pl^2, inside gate M's band, but the sign still fails
     under that cutoff.

RESULT. With its searchlight label on one list, the bridge cannot be
light, and neither can the 7 G2 scalars (the anchor moves by 1.6e13
and 2.02). Without them the covariant contact stiffness is -4.5, which
is repulsive. It is attractive only for a vector contact below 3/4.
Given a mass operator, the labelled bridge at the boundary leaves the
running intact. It enters the stiffness with a heavy row weight set by
the regulator's shape. An attractive list whose anchored regulator
lies inside the boundary needs 0.370 of a massless row's weight there.
The ledger's own sharp cutoff gives 0.148 (Sigma = -0.175, repulsive).
exp(-t0/t) gives 0.280 (Sigma = +3.65, attractive, but its anchored
regulator lies outside the boundary). The G2-singlet content the
repository names supplies 2.34 of the 6.46 the sign needs under the
ledger's own cutoff. So for the labelled bridge the sign depends on
the regulator's shape at the boundary, and the regulator literature
leaves that shape unsettled. Read as amount beyond the horizon, as the
gravity.py ledger counts it and boundary.py's one list reads it, the
bridge takes no share of the running. The sign then holds for any
cutoff shape (tests/probes/horizon_amount.py).

Usage: python3 tests/probes/boundary_content.py
"""

import math
import os
import sys
from fractions import Fraction

from scipy.optimize import brentq
from scipy.special import expn, k1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

XI = float(r.xi_bridge)
ALPHA = r.alpha_G2_Pl
S0 = r.S_lepton
C2 = float(r.C2_26)
CONTACTS = (1.0, 0.0, -1.0)
ANCHOR_PPM = 22.2                       # the anchor's M_Pl residual (g_path_trials.py)
CARRIER = 1.0                           # six real components x 1/6 (composite_pole.py K5)
METRIC_LOOP = -8.0 / 3.0                # on shell (ledger_kernel.py R5)


def rows(bridge="boundary", scalars="composite"):
    """boundary.py's rows. By default the labelled bridge at the boundary."""
    return B.rows(r.d10, r.d11, r.n7, r.n26, bridge=bridge, scalars=scalars)


def a1(kind, c):
    return B.a1(kind, c, XI)


b0 = B.b0


def sigma(table, t0, c, m_over_Mb=1.0, heavy=None):
    """Stiffness at the regulator t0 (in M_b^-2): light 1, boundary its heavy row weight, else 0."""
    return B.sigma(table, t0, c, XI, m_over_Mb, heavy)


def fixed_points(table, c, m_over_Mb=1.0, extra=0.0):
    """Positive solutions of Sigma = sigma(t0 = Sigma/2 pi): the regulator solved from the anchored G."""
    return B.anchored(table, c, XI, m_over_Mb, extra)


def smooth_weight(x):
    """Heavy row weight under rho = exp(-t0/t), relative to a massless row: 2 sqrt(x) K1(2 sqrt(x)), x = m^2 t0."""
    return 2 * math.sqrt(x) * float(k1(2 * math.sqrt(x)))


def light_shift(b0_new):
    """M_Pl and eta_B factors and the G2 endpoint when the running counts b0_new."""
    S_new = 2 * math.pi / (float(b0_new) * ALPHA)
    return dict(b0=b0_new, M_Pl=math.exp(S_new - S0),
                M_Pl_mu_star=math.exp((S0 - C2) * (1 - float(b0_new / r.b0_G2))),
                eta_B=math.exp(-float(Fraction(r.hv_G2, r.hv_F4)) * (S_new - S0)),
                alpha_star=2 * math.pi / (float(b0_new) * C2))


def run(report=print):
    out = {}

    # B1: the bridge is n26 real 7s of the running G2. The 7 scalars are G2's breaking.
    index_E8 = r.hv_G2 + r.n26 * r.T_7
    v_break = math.exp(-r.S_quark) * math.sqrt(float(r.E_v2))
    assert index_E8 == r.hv_E8 == 30
    assert b0(rows()) == r.b0_G2 == Fraction(32, 3)
    booked = rows(bridge="amount", scalars="light")           # the gravity.py ledger: the bridge as amount, xi = 0
    both_light = rows(bridge="light", scalars="light")
    sig_booked = sigma(booked, 0.0, 1.0)
    assert abs(sig_booked - 27.0) < 1e-12 and 1.3e-17 < v_break < 1.5e-17
    out.update(index_E8=index_E8, v_break=v_break, sigma_booked=sig_booked)
    report(f"B1 T(248)|G2 = {r.hv_G2} + {r.n26} T(7) = {index_E8} = h(E8): read with this label the bridge is "
           f"{r.n26} real 7s of the running G2. The 7 scalars are G2's breaking at v_break = {v_break:.1e} M_b, no row by R6. "
           f"The gravity.py ledger's Sigma = {sig_booked:.4f} counts both, b0 = {r.b0_G2} counts neither.")

    # B2: light in both.
    shifts = {"bridge": light_shift(b0(rows(bridge="light"))),
              "G2 scalars": light_shift(b0(rows(scalars="light"))),
              "both": light_shift(b0(both_light))}
    assert shifts["bridge"]["b0"] == Fraction(19, 3) and shifts["G2 scalars"]["b0"] == Fraction(21, 2)
    assert shifts["both"]["b0"] == Fraction(37, 6)
    assert 1e13 < shifts["bridge"]["M_Pl"] < 2e13 and abs(shifts["G2 scalars"]["M_Pl"] - 2.02) < 0.01
    assert min(v["M_Pl_mu_star"] for v in shifts.values()) > 1.8
    assert abs(light_shift(r.b0_G2)["alpha_star"] - math.pi / 32) < 1e-12
    out["light"] = shifts
    for k, v in shifts.items():
        report(f"B2 {k:10s} light in both: b0 = {str(v['b0']):5s} M_Pl x{v['M_Pl']:.3g} (x{v['M_Pl_mu_star']:.3g} "
               f"leaving at mu*), eta_B x{v['eta_B']:.3g}, alpha(mu*) = {v['alpha_star']:.4f} (pi/32 = {math.pi / 32:.4f})")

    # B3: absent from both.
    absent = rows(bridge="absent")
    sig_absent = {c: sigma(absent, 0.0, c) for c in CONTACTS}
    assert all(abs(sig_absent[c] - (13.5 - 18 * c)) < 1e-12 for c in CONTACTS)
    assert sig_absent[1.0] + CARRIER < 0 and sig_absent[1.0] + CARRIER + METRIC_LOOP < 0
    out["sigma_absent"] = sig_absent
    report(f"B3 absent: Sigma = 13.5 - 18 c = {sig_absent[1.0]:+.2f}/{sig_absent[0.0]:+.2f}/{sig_absent[-1.0]:+.2f} "
           f"(c = +1/0/-1). Attractive only for c < 3/4. With the carrier {sig_absent[1.0] + CARRIER:+.2f}, "
           f"and its metric loop {sig_absent[1.0] + CARRIER + METRIC_LOOP:+.2f}")

    # B4: the regulator lies inside the boundary iff Sigma >= 2 pi.
    t0_booked = sig_booked / (2 * math.pi)
    out.update(t0_booked=t0_booked, regulator_booked=t0_booked**-0.5)
    assert abs(t0_booked - 4.297) < 1e-3 and t0_booked**-0.5 < 1 and r.d11**-0.5 < 1
    report(f"B4 t0 = Sigma/(2 pi) M_b^-2: the regulator is inside the boundary iff Sigma >= 2 pi; "
           f"gravity.py ledger t0 = {t0_booked:.3f}, regulator at {t0_booked**-0.5:.3f} M_b (t0 = d11: {r.d11**-0.5:.3f} M_b)")

    # B5: the bridge at the boundary.
    edge = math.exp(-ANCHOR_PPM * 1e-6 * float(r.b0_G2 / (r.b0_G2 - shifts["bridge"]["b0"])))
    at_b = rows(bridge="boundary")
    unit = r.n7 * r.n26 * (1 / 6 - XI)                            # the bridge at full weight
    w_sharp, w_smooth = float(expn(2, 1.0)), smooth_weight(1.0)     # weight at m_B^2 t0 = 1
    sig_inside = sigma(at_b, 1.0, 1.0)
    sig_smooth = sigma(at_b, 1.0, 1.0, heavy=w_smooth)
    need = 2 * math.pi - sig_inside
    w_need = (2 * math.pi - sig_absent[1.0]) / unit
    w_need_loop = (2 * math.pi - sig_absent[1.0] - METRIC_LOOP) / unit
    w_sign = -sig_absent[1.0] / unit
    f_smooth = lambda s: sigma(at_b, s / (2 * math.pi), 1.0, heavy=smooth_weight(s / (2 * math.pi))) - s
    s_smooth = brentq(f_smooth, 1e-6, 2 * math.pi)
    reg_smooth = (s_smooth / (2 * math.pi))**-0.5
    free = {c: fixed_points(at_b, c) for c in CONTACTS}
    assert all(len(v) == 1 for v in free.values())
    s1, s0 = free[1.0][0], free[0.0][0]
    assert abs(edge - 0.999945) < 2e-6 and abs(unit * w_sharp - 4.325) < 1e-3 and -0.2 < sig_inside < 0
    assert w_sharp < w_smooth < w_need < w_need_loop and abs(w_need - 0.3702) < 1e-3 and 0 < sig_smooth < 2 * math.pi
    assert w_sharp < w_sign < w_smooth and f_smooth(2 * math.pi) < 0 and reg_smooth > 1
    assert abs(s1 - 3.698) < 1e-3 and (s1 / (2 * math.pi))**-0.5 > 1.3 and 6 * math.pi / s1 > 5
    assert abs(s0 - 14.271) < 1e-3 and s0 > 2 * math.pi
    assert max(fixed_points(at_b, 1.0, extra=need + 0.01)) > 2 * math.pi
    out.update(anchor_edge=edge, bridge_cap=unit * w_sharp, sigma_inside=sig_inside, needed=need,
               needed_loop=need - METRIC_LOOP, needed_scalars=math.ceil(6 * need), w_sharp=w_sharp,
               w_smooth=w_smooth, sigma_smooth=sig_smooth, w_need=w_need, w_need_loop=w_need_loop,
               w_sign=w_sign, sigma_smooth_free=s_smooth, regulator_smooth_free=reg_smooth,
               free={c: v[0] for c, v in free.items()}, t0_free=s1 / (2 * math.pi),
               regulator_free=(s1 / (2 * math.pi))**-0.5, ratio_free=6 * math.pi / s1)
    report(f"B5 bridge at the boundary: running intact for m_B >= {edge:.6f} M_b. Sharp cutoff, regulator inside: "
           f"bridge <= {unit * w_sharp:.3f}, Sigma <= {sig_inside:+.3f} (c = +1), {need:.2f} short of 2 pi "
           f"({need - METRIC_LOOP:.2f} with the metric loop, {math.ceil(6 * need)} real scalars)")
    report(f"   weight needed {w_need:.3f} of a massless row ({w_need_loop:.3f} with the loop. The sign alone "
           f"{w_sign:.3f}). Sharp {w_sharp:.3f}, exp(-t0/t) {w_smooth:.3f} (Sigma = {sig_smooth:+.2f}, "
           f"{2 * math.pi - sig_smooth:.2f} short. Solved freely, Sigma = {s_smooth:.3f}, regulator at "
           f"{reg_smooth:.3f} M_b)")
    report(f"   solved freely (sharp): Sigma = {s1:.3f}, t0 = {s1 / (2 * math.pi):.3f} l_Pl^2, regulator at "
           f"{(s1 / (2 * math.pi))**-0.5:.3f} M_b, 6 pi/Sigma = {6 * math.pi / s1:.2f} (c = +1). c = 0: Sigma = "
           f"{s0:.3f}, t0 = {s0 / (2 * math.pi):.3f}. c = -1: {free[-1.0][0]:.3f}")

    # B6: the boundary coupling treated as induced (shared cutoff).
    b_m = float(Fraction(4, 3) * r.n_dirac_G2 * r.T_7)
    b_m_booked = b_m - float(shifts["both"]["b0"] - r.b0_G2)
    log_needed = 2 * math.pi / (b_m * ALPHA)
    assert abs(log_needed - 12 * math.pi**2) < 1e-9 and b_m_booked == 8.5
    out.update(b_matter=b_m, induced_log=log_needed, induced_log_booked=2 * math.pi / (b_m_booked * ALPHA))
    report(f"B6 induced 1/alpha_G2(M_b) = (b_m/2 pi) ln(Lambda/M_b), b_m = {b_m:.0f}: 24 pi needs "
           f"ln(Lambda/M_b) = {log_needed:.1f} ({out['induced_log_booked']:.1f} on the gravity.py list), against a regulator "
           f"at {t0_booked**-0.5:.2f} M_b")

    # B7: what the sign needs, against the G2-singlet content the repository names.
    named = {"carrier": CARRIER, "F4 singlets light": 4 / 12, "G2 scalars at M_b": 7 / 6 * w_sharp,
             "bridge at the mass edge": unit * (float(expn(2, edge**2)) - w_sharp),
             "radial mode light": 1 / 6, "Higgs doublet light": 4 / 6}
    stack = sum(named.values())
    t0_f4 = fixed_points(at_b, 1.0, extra=named["F4 singlets light"])[0] / (2 * math.pi)
    assert abs(stack - 2.341) < 2e-3 and abs(need - stack - 4.117) < 2e-3 and 0.6 < t0_f4 < 1
    out.update(named=named, stack=stack, stack_short=need - stack, needed_weyl=math.ceil(12 * need),
               t0_f4_singlets=t0_f4)
    report(f"B7 under the ledger's own cutoff the labelled bridge's sign needs {need:.2f} at the boundary ({math.ceil(6 * need)} light G2-singlet real scalars, "
           f"{math.ceil(12 * need)} Weyl states). The named G2-singlet content adds {stack:.2f} at its most "
           f"favourable weights, {need - stack:.2f} short. The two F4 singlets light alone: t0 = {t0_f4:.3f} l_Pl^2")
    report(f"RESULT: read with its label, the bridge at the boundary needs {w_need:.3f} of a massless row's weight "
           f"for an attractive list with its anchored regulator inside the boundary, which the sharp cutoff "
           f"({w_sharp:.3f}, repulsive) and exp(-t0/t) ({w_smooth:.3f}, attractive, regulator outside) do not give. "
           "Its sign depends on the regulator's shape. Read as amount it takes no share (horizon_amount.py).")
    return out


if __name__ == "__main__":
    run()
