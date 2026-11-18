"""The common mode on the recursive web: gravity's strength as the kernel's fixed point.

THE QUESTION. The web holds the ruler as a node, M_Pl_MeV: the anchor
inversion at web depth four. It holds the venting ledger's rows at
depth one. It holds the common mode's own coefficient G M_b^2 as the
node common_mode (root.py, THE COMMON MODE). That node has base 0 and
two edges at depth four: the relative channels' powers of the lepton
circulant on its depth four rates (common_mode_count.py). On the
searchlight's circulant these edges give the whole that
boundary_whole.py assembles by hand. Register the node alone on the
engine's web, in the engine's grammar, and solve with the kernel
x <- b + W(x). What value does the node take? What does the web then
give for gravity with no dimensional input?

PREMISES. P1 the engine's web after the full chain (root, masses,
mixing, couplings, gravity). P2 the echo law: "The common mode record
of what does not cancel locally is the coherent venting that IS
gravity" (interference/root.py). P3 the vacuum is the empty field, so
the common mode has no base of its own. P4 the forgetting: the first
moment of the harmonic cancels, P v P = 0, and each relative channel
keeps the power |v_q|^2
(interference/searchlight/protected_forgetting.py).

THE NODES (registered for the run and removed after it; the engine's
own common_mode is held and put back).
  common_mode     base 0, two edges at depth four with parent eq_loop,
                  the powers of the relative channels q = 1 and q = 2
                  of the lepton circulant, as the engine's node
  alpha_G(e,e)    the exchange between two electrons through the
                  common mode: common_mode x (m_e / M_Pl_MeV)^2, both
                  rates read from the state
  common_mode(5)  base 0, one edge at depth five whose parent is the
                  formed tau ledger: the second moment of the lepton
                  closure as the web solves it, echoes included

STATEMENTS (the numbers below check each one).

  K1 (the fixed point). The kernel converges with
     common_mode = 3 Q_l - 1 = 1 - 7.50 (alpha/2pi)^2. The node adds
     nothing of its own. On the searchlight's circulant each relative
     channel gives one half, the forgetting's powers.
  K2 (gravity with no dimensional input). alpha_G(e,e) = 1.75171e-45,
     a number fixed by the four integers alone. Moving the anchor m_e
     by a part in a thousand moves M_Pl_MeV by the same fraction and
     leaves alpha_G unchanged to 1e-12. The measured G m_e^2 / hbar c
     is 1.75181e-45. The web's value is 54.6 ppm below it. This is the
     anchor's test of G, in a form that does not depend on the anchor.
  K3 (the node is a reading). One relative channel alone gives 1/2 on
     the searchlight's circulant: the bridge's reading, with c = 2. A
     base of its own (any amplitude of the empty field) adds to the
     node one for one.
  K4 (the web's own second moment). Read from the formed lepton
     ledgers, the closure's second moment is 1 at depths one and two
     and 1 - 1.0e-5 at depth four. The depth four backreaction is
     proportional to each knot's own amplitude, so it does not cancel
     in the moment. The engine's count is that moment. The formed
     ledgers' other factors are common to the three rates and are not
     records (gate U). So the depth five reading off the formed tau
     ledger and the depth four count agree to 1e-15.
  K5 (one reading). The engine's G and the formed closure's are one
     number, 6.673936e-11, 54.6 ppm below the recommended value.

RESULT. On the engine's web the common mode is a node with no base and
two edges: the relative powers of the lepton circulant at depth four.
The kernel returns 3 Q_l - 1 = 1 - 7.50 (alpha/2pi)^2, the formed
closure's own second moment. Gravity's strength between two electrons
is a pure number of the web, 54.6 ppm below the measured one.

Scope: the edges on the searchlight's circulant are the forgetting's
identity. The count is that of common_mode_count.py. That the node
is G M_b^2 is the boundary identification (boundary_whole.py). No exchange at the boundary is
computed.

Usage: python3 tests/probes/common_mode_web.py     (about 20 s, runs the chain)
"""

import contextlib
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402
from interference.searchlight import protected_forgetting as pf  # noqa: E402

LEPTONS = ("m_e", "m_mu", "m_tau")
NODES = ("common_mode", "alpha_G(e,e)", "common_mode(5)")
G_SIGMA = r.G_REL_SIGMA                      # relative uncertainty of the recommended G
CANDIDATE = ("FORCED", "CANDIDATE")


def _chain():
    from interference import couplings, gravity, masses, mixing
    with contextlib.redirect_stdout(io.StringIO()):
        R = r.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
        c = couplings.derive(R, m)
        gravity.derive(R, m, mx, c)


def channel_powers(v):
    """|v_q|^2 for the relative channels q = 1, 2 of a three vector, and its mean, the q = 0 part."""
    f = np.fft.fft(np.asarray(v, dtype=float)) / 3.0
    return float(abs(f[1]) ** 2), float(abs(f[2]) ** 2), float(f[0].real)


def second_moment(rates):
    """P v^2 P of a closure with the three rates given: v = sqrt(rate) over its mean, less one."""
    a = np.sqrt(np.asarray(rates, dtype=float))
    v = a / a.mean() - 1.0
    return float(np.mean(v * v))


_ENGINE = {}                                  # the engine's own nodes of these names, held while a probe runs


def register(web, channels=(1, 2), base=0.0):
    """Put the common mode's nodes on the web, holding the engine's own. Returns the ledgers."""
    for name in NODES:
        if name in web and name not in _ENGINE:
            _ENGINE[name] = (web[name], web.state.get(name))
    p1, p2, first = channel_powers(pf.harmonic())
    assert abs(first) < 1e-12                                     # P v P = 0: nothing vents at first order
    searchlight = dict(web.state, eq_loop=0.0)
    assert all(abs(r._lepton_z3_power(searchlight, q) - p) < 1e-14 for q, p in ((1, p1), (2, p2)))
    cm = r.Ledger("common_mode", base, "add", web_depth=4)
    for q in channels:
        cm.echo([f"Z3(q={q})"], lambda s, q=q: r._lepton_z3_power(s, q), 4, "FORCED",
                "a relative channel's power on the depth four lepton circulant", kind="vent", parent="eq_loop")
    web["common_mode"] = cm
    web["alpha_G(e,e)"] = r.Ledger(
        "alpha_G(e,e)", lambda s: s["common_mode"] * (s["m_e"] / s["M_Pl_MeV"]) ** 2, "mul", web_depth=4)
    formed = r.Ledger("common_mode(5)", 0.0, "add", web_depth=5)
    formed.echo(["read:lepton closure"], lambda s: second_moment([s[k] for k in LEPTONS]), 5, "CANDIDATE",
                "the formed closure's second moment", kind="vent", parent="m_tau")
    web["common_mode(5)"] = formed
    return cm, formed


def remove(web):
    """Take the probe's nodes off the web and put the engine's own back."""
    for name in NODES:
        web.pop(name, None)
        web.state.pop(name, None)
    for name, (led, value) in _ENGINE.items():
        web[name] = led
        web.state[name] = value
    _ENGINE.clear()


def solved(channels=(1, 2), base=0.0):
    """Register, run the kernel, read the three nodes, and remove them."""
    web = r.WEB
    try:
        _, formed = register(web, channels, base)
        web.solve()
        return dict(common_mode=web.state["common_mode"], alpha=web.state["alpha_G(e,e)"],
                    searchlight=web["common_mode"].value_at(dict(web.state, eq_loop=0.0)),
                    formed=formed.value_at(web.state, statuses=CANDIDATE), m_pl=web.state["M_Pl_MeV"],
                    depths={n: (web[n].web_depth, max((t.depth for t in web[n].terms), default=0)) for n in NODES})
    finally:
        remove(web)
        web.solve()


def run(report=print, build=True):
    out = {}
    if build and "Sigma_vent(UV)" not in r.WEB:
        _chain()
    web = r.WEB
    before = dict(web.state)

    # K1: the fixed point.
    one = solved()
    L = web.state["eq_loop"]
    assert abs((one["common_mode"] - 1.0) / L + 7.502) < 0.001 and one["depths"]["common_mode"] == (4, 4)
    assert abs(one["searchlight"] - 1.0) < 1e-12 and abs(one["common_mode"] - r.COMMON_MODE) < 1e-15
    out.update(common_mode=one["common_mode"], vents=one["searchlight"], depths=one["depths"])
    report(f"K1 the kernel's fixed point: common_mode = {one['common_mode']:.12f} = 1 "
           f"{(one['common_mode'] - 1) / L:+.3f} L, base 0, two edges summing to {one['searchlight']:.0f} on the "
           f"searchlight's circulant")

    # K2: gravity with no dimensional input.
    measured = (r.M_E_ANCHOR_MEV / (r.M_PL_CODATA_GEV * 1e3)) ** 2
    ppm = 1e6 * (one["alpha"] / measured - 1)
    anchor = r.M_E_ANCHOR_MEV
    try:
        r.M_E_ANCHOR_MEV = anchor * 1.001
        moved = solved()
    finally:
        r.M_E_ANCHOR_MEV = anchor
        web.solve()
    assert abs(moved["m_pl"] / one["m_pl"] - 1.001) < 1e-12 and abs(moved["alpha"] / one["alpha"] - 1) < 1e-12
    assert abs(ppm + 54.6) < 0.6
    out.update(alpha_G=one["alpha"], alpha_G_measured=measured, ppm=ppm,
               anchor_moves_alpha=moved["alpha"] / one["alpha"] - 1)
    report(f"K2 alpha_G(e,e) = {one['alpha']:.5e} against the measured {measured:.5e}, {ppm:+.1f} ppm. The anchor "
           f"up by 1e-3 moves M_Pl by {moved['m_pl'] / one['m_pl'] - 1:+.1e} and alpha_G by "
           f"{moved['alpha'] / one['alpha'] - 1:+.0e}")

    # K3: the node is a reading.
    half, based = solved(channels=(1,)), solved(base=0.1)
    assert abs(half["searchlight"] - 0.5) < 1e-12 and abs(based["searchlight"] - 1.1) < 1e-12
    out.update(one_channel=half["searchlight"], with_base=based["searchlight"])
    report(f"K3 one relative channel alone gives {half['searchlight']:.2f} on the searchlight's circulant; a base "
           f"of 0.1 gives {based['searchlight']:.2f}")

    # K4: the web's own second moment, by depth.
    by_depth = {d: second_moment([web[k].value(max_depth=d) for k in LEPTONS]) for d in (1, 2, 4)}
    g_ppm = 1e6 * (r.G_PRED / r.G_CODATA - 1)
    formed_ppm = 1e6 * (r.G_DEPTH_ONE * one["formed"] / r.G_CODATA - 1)
    assert abs(by_depth[1] - 1) < 1e-12 and abs(by_depth[2] - 1) < 1e-12
    assert abs(by_depth[4] - one["formed"]) < 1e-14 and abs((by_depth[4] - 1) + 1.0e-5) < 2e-6
    assert abs(one["formed"] - r.COMMON_MODE) < 1e-15 and one["depths"]["common_mode(5)"] == (5, 5)
    out.update(moment_by_depth=by_depth, formed=one["formed"], g_ppm=g_ppm, g_ppm_formed=formed_ppm)
    report("K4 the lepton closure's second moment at depths 1, 2, 4: "
           + ", ".join(f"1{by_depth[d] - 1:+.2e}" for d in (1, 2, 4))
           + f". Read at depth five it puts G at {formed_ppm:+.1f} ppm "
           f"({formed_ppm * 1e-6 / G_SIGMA:+.1f} standard deviations), the count's {g_ppm:+.1f}")

    # K5: the formed closure reading.
    G1, G5 = r.G_PRED, r.G_DEPTH_ONE * one["formed"]
    out.update(G_engine=G1, G_formed=G5, apart_ppm=1e6 * (G1 / G5 - 1))
    assert abs(out["apart_ppm"]) < 1e-6
    report(f"K5 G = {G1:.6e} as the engine reads it and {G5:.6e} read from the formed closure, "
           f"{out['apart_ppm']:.1e} ppm apart")

    # the web is as it was
    assert all(n not in web for n in NODES if n != "common_mode") and web["common_mode"].web_depth == 4
    assert all(abs(web.state[k] - v) <= 1e-12 * max(1.0, abs(v))
                                                    for k, v in before.items())
    report("RESULT: on the engine's web the common mode is a node with no base and two edges, the relative powers "
           "of the lepton circulant at depth four. The kernel returns 3 Q_l - 1 = 1 - 7.50 (alpha/2pi)^2, the "
           "formed closure's own second moment. alpha_G(e,e) is a pure number of the web, 54.6 ppm below the "
           "measured one.")
    return out


if __name__ == "__main__":
    run()
