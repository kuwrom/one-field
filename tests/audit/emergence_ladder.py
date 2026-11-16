"""
emergence_ladder.py: the record of couplings per rung, read off the web.

The one field appears at several levels of description, and each level
closes and persists for its own reason.  This module lists the levels
from the bottom up as RUNGS.  For each rung it records what closes
there, what keeps it, and which edges of the web couple it to the other
rungs.  The word "rung" is used because "layer" already names the four
reading layers (ENGINE, SEARCHLIGHT, SUBSTRATE, AUDIT;
interference/__init__.py).  A rung is what a thing IS at one level of
emergence.  A layer is where it is READ from.

GENERATED from root.WEB after the full chain:
  - every node, with its value and web depth;
  - every echo and vent term, assigned to the rungs its channel path
    runs through;
  - every ratio edge (a live base that reads another node), found
    through the web's own read tracker.
The coupling record and the common mode column are computed, never
written by hand.

CURATED, per rung: what closes, what keeps it, its spins, and its
counterfactuals.  The counterfactuals are ids in candidates.tsv, and
their verdicts are read from that index rather than restated.

Two placement rules keep the record complete: every web node sits on
exactly one rung, and every channel maps to a rung.  A new node or a new
channel fails tests/audit/test_emergence_ladder.py until it is placed.

WHAT THE RECORD SHOWS (each point is frozen by the test):
  1. The common mode.  Sigma is the node that every rung carrying
     degrees of freedom vents into (bridge, forces, leptons, quarks),
     all at depth 1.  This is the ledger form of universality (gate U).
     Without the bridge row the sum is -10/3: the same particles would
     induce gravity of the wrong sign.  This node is the ledger as
     counted, the reference for FACT 7.  The gates read the one list of
     interference/boundary.py, which has the same rows minus the 7 G2
     scalars.  There the bridge's +30.33, read as amount past the
     searchlight's horizon, passes the sign gate (S) at every contact
     (tests/probes/horizon_amount.py).
  2. One path for the deepest backreaction.  Every depth-4 edge has as
     its parent the one depth-3 node, the e<->q charged loop.  The
     depth-4 edges are 1/alpha, the three lepton masses, S_quark, and
     the common mode's two relative powers of the lepton circulant.
  3. Sparse rungs.  Quark masses meet the leptons only by ratio (u and d
     from e, c from mu, t from tau), and meet the algebra only through
     depth-1 WZW vents.  Lepton masses take electromagnetic backreaction
     at depths 2 and 4.  Constants folded into a base (the four
     integers, theta) are not edges.  So an empty cell means there is no
     edge the solver iterates, not that there is no dependence.
  4. Where the web stops.  Five rungs carry no node.  The closure
     algebra enters as channel names and base constants.  The substrate
     is read from the probes.  Hadrons, nuclei and atoms are read from
     QCD and QED with the engine's couplings and masses.  The
     cosmological constant is an integration constant.
  5. The gravity rung's gates.  Each of the lane's six gates appears
     with its verdict, its tier and the number of tests that witness
     it, read from the engine (interference/gravity.py).  The nodes of
     other rungs that read the common mode are listed beside them.

Usage: python3 tests/audit/emergence_ladder.py
"""

import contextlib
import csv
import io
import os
import sys
import textwrap

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_here), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import (couplings, dark_sector, gravity,  # noqa: E402
                          masses, mixing, root)

CANDIDATES = os.path.join(_here, "candidates.tsv")


def _rung(key, name, closes, keeps, spins, nodes=(), read="", counter=()):
    return dict(key=key, name=name, closes=closes, keeps=keeps, spins=spins,
                nodes=tuple(nodes), read=read, counter=tuple(counter))


RUNGS = (
    _rung("substrate", "substrate",
          "the spinning lump: one quantum of common phase winding around a "
          "ring, with the Z3 clock locked inside (committed self-bound sign)",
          "the Noether charges N and L_z = mN. The m = 1 ring holds above a "
          "threshold charge, and more winding or less charge splits it. No "
          "pi_3 knot fits inside the lump",
          "0: three complex components in cyclic Z3 exchange",
          read="tests/probes/spinning_lump.py, knot_in_lump.py",
          counter=("sub-uniform-bdg", "sub-complex-alpha", "sub-long-root-v")),
    _rung("algebra", "closure algebra",
          "the finite closure menu: the SU(3)_3 data (d10 = 2, d11 = 3, "
          "h10 = 2/9) and the conformal embedding E8(1) > G2(1) x F4(1)",
          "integrality: the braid image is finite, so the menu is discrete",
          "none. It supplies the WZW vents",
          read="interference/root.py",
          counter=("theta-menu", "theta-2pih", "theta-extremum")),
    _rung("bridge", "bridge",
          "182 = 7 x 26 coherence channels between the G2 and F4 faces",
          "it stays delocalized: 1/(2pi) of its current response localizes "
          "as baryons, and the rest remains coherence",
          "0 (h_bridge = 1 is integer spin; the Hurwitz ceiling forces the "
          "scalar reading)",
          nodes=("bridge_sq", "a1_bridge"),
          counter=("g-bridge-fermion",)),
    _rung("scales", "scale hierarchy",
          "the closure actions that place each scale at e^-S, as in "
          "v_EW = M_Pl e^-S_quark",
          "decoupling: rungs a factor e^-S apart average each other out "
          "and interact only through echoes tagged with a depth",
          "none",
          nodes=("S_lepton", "S_baryo", "S_quark", "v_EW_GeV"),
          counter=("mode-25", "inst-measure", "as-2loop-g2")),
    _rung("forces", "gauge forces and Higgs",
          "the couplings: 1/alpha from the bridge self-response and the "
          "e<->q loop, sin^2 theta_W, the Higgs quartic at M_Pl",
          "conserved charges: states of different charge sit in different "
          "sectors and never interfere",
          "1 (photon, W, Z, gluons, coset vectors), 0 (Higgs, G2 scalars)",
          nodes=("eq_loop", "inv_alpha", "lepton_EM", "sin2W", "lambda_MPl"),
          counter=("H-veltman", "H-mass-elementary", "kg-threshold-field")),
    _rung("leptons", "leptons (G2 face)",
          "e, mu, tau: three standing waves of the G2 face at theta = "
          "h(fund) = 2/9 with B/A = sqrt 2",
          "species is a conserved node label, since no edge kind changes a "
          "knot's species (root.py, WHY THE PROTON IS STABLE). The "
          "electron is the lightest charged state",
          "1/2",
          nodes=("m_e", "m_mu", "m_tau"),
          counter=("theta-pi3", "wall-triple", "g2-I6-hessian")),
    _rung("quarks", "quarks (F4 face)",
          "u, d, s, c, b, t: word counts on the F4 face, with the crossing "
          "Koide triplets (c, b, t) and (s, c, b)",
          "species and baryon number are conserved node labels. Colour "
          "confines quarks to hadrons",
          "1/2",
          nodes=("m_u", "m_d", "m_s", "m_c", "m_b", "m_t", "Q_cbt", "Q_scb"),
          counter=("ud-bare", "c-ladder", "t-yukawa", "b-q0", "same-charge")),
    _rung("hadrons", "hadrons and nuclei",
          "proton, neutron and nuclei: QCD bound states of the quark rung",
          "the proton is the lightest baryon and no edge kind changes "
          "baryon number: no proton decay at any sensitivity (root.py)",
          "1/2 and composite",
          read="QCD with the engine's alpha_s and quark masses. Most of the "
               "nucleon mass is QCD field energy, not quark mass"),
    _rung("atoms", "atoms",
          "electrons bound to nuclei through the electromagnetic coupling",
          "both constituents are the lightest carriers of a conserved label "
          "(charge, baryon number), and the bound state is gapped",
          "composite",
          read="QED with the engine's alpha and m_e"),
    _rung("gravity", "gravity (common mode)",
          "G from the anchor inversion (M_Pl from m_e) read through the "
          "common mode, G M_b^2 = 3 Q_l - 1. The induced Einstein "
          "term Sigma, every rung's heat kernel weight summed, as its "
          "diagnostic",
          "relative phases between rungs never reach it, so it only adds. "
          "It is universal (gate U), and attractive on the one list, where "
          "the bridge vents as amount (gate S)",
          "2 (the composite metric, a gradient bilinear of the carrier)",
          nodes=("M_Pl_MeV", "common_mode", "Sigma_vent(UV)",
                 "Sigma_vent(broken)", "G_ratio"),
          counter=("g-min-spectrum", "g-screened-poisson", "g-bridge-vector",
                   "g-vortex-tangle", "g-graviton-ondd", "g-one-list")),
    _rung("dark", "dark sector",
          "the bridge coherence that never localizes: dark/baryonic = "
          "2pi - 1",
          "no electromagnetic channel: it couples only through the common "
          "mode",
          "none",
          nodes=("dark_ratio",)),
    _rung("lambda", "cosmological constant",
          "the vacuum term: a constant potential decouples from the "
          "traceless composite metric, so Lambda is an integration constant",
          "nothing to keep: a boundary datum, not a closure",
          "none",
          read="tests/probes/composite_metric.py",
          counter=("L-full-form", "L-ckn", "L-sequester")),
)

# Channel path element (by prefix) -> the rung(s) it runs through.
CHANNELS = (
    ("G₂×WZW", ("leptons",)), ("G₂-face", ("leptons",)),
    ("G₂-scal", ("forces",)), ("WZW", ("algebra",)),
    ("F4-knot", ("quarks",)), ("bridge", ("bridge",)), ("Z3(q=", ("algebra",)),
    ("read:eq_loop", ("forces",)), ("e<->q", ("forces",)),
    ("EM", ("forces",)), ("Higgs", ("forces",)), ("gauge", ("forces",)),
    ("W, Z", ("forces",)), ("coset", ("forces",)),
    ("fermions", ("leptons", "quarks")),
)

KEYS = tuple(r["key"] for r in RUNGS)
NODE_RUNG = {n: r["key"] for r in RUNGS for n in r["nodes"]}


def channel_rungs(element):
    for prefix, rungs in CHANNELS:
        if element.startswith(prefix):
            return rungs
    return None


GATES = {}                                       # the gravity lane's gates, filled by run_chain


def run_chain():
    """The full derivation chain, silently. Returns the solved web and keeps the lane's gates."""
    with contextlib.redirect_stdout(io.StringIO()):
        R = root.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
        c = couplings.derive(R, m)
        g = gravity.derive(R, m, mx, c)
        dark_sector.derive(R, g)
    GATES.clear()
    GATES.update(g["gravity_gates"])
    return root.WEB


def candidate_verdicts(path=CANDIDATES):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and not row[0].startswith("#") and row[0] != "id":
                out[row[0]] = (row[4], row[5])
    return out


def _reads(fn, state):
    tracker = root.Web._ReadTracker(state)
    try:
        fn(tracker)
    except root.Web._SOFT:
        pass
    return tracker.reads


def edges_of(web):
    """Every coupling into every node: echo and vent terms by channel,
    ratio edges by the nodes a live base or factor reads."""
    edges, unplaced = [], set()
    for node, led in web.items():
        if not isinstance(led, root.Ledger):
            continue
        tgt = NODE_RUNG.get(node)
        for t in led.terms:
            srcs = []
            for element in t.path:
                rungs = channel_rungs(element)
                if rungs is None:
                    unplaced.add(element)
                else:
                    srcs += [s for s in rungs if s not in srcs]
            edges.append(dict(node=node, tgt=tgt, srcs=tuple(srcs),
                              depth=t.depth, kind=t.kind,
                              path="->".join(t.path), parent=t.parent,
                              factor=t.value(web.state)))
        if callable(led.base):
            for read in sorted(_reads(led.base, web.state) - {node}):
                edges.append(dict(node=node, tgt=tgt,
                                  srcs=(NODE_RUNG.get(read, "?"),),
                                  depth=None, kind="ratio",
                                  path=f"reads {read}", factor=None))
    return edges, sorted(unplaced)


def coupling_record(edges):
    """(source rung, target rung) -> sorted labels: web depths, r = ratio."""
    rec = {}
    for e in edges:
        label = "r" if e["kind"] == "ratio" else str(e["depth"])
        for s in e["srcs"]:
            rec.setdefault((s, e["tgt"]), set()).add(label)
    return {k: sorted(v) for k, v in rec.items()}


def sources_by_node(edges):
    """node -> the set of OTHER rungs that feed it."""
    out = {}
    for e in edges:
        out.setdefault(e["node"], set()).update(
            s for s in e["srcs"] if s != e["tgt"])
    return out


def common_mode(web, node="Sigma_vent(UV)"):
    rows = []
    for t in web[node].terms:
        rows.append((t.path[0], channel_rungs(t.path[0]) or ("?",),
                     t.value(web.state)))
    return rows, web[node].value()


def missing_sources():
    """Repository paths a rung says it is read from, that do not exist."""
    repo, missing = os.path.dirname(os.path.dirname(_here)), []
    for r in RUNGS:
        folder = ""
        for token in (t.strip() for t in r["read"].split(",")):
            if not token.endswith(".py"):
                continue
            folder = os.path.dirname(token) or folder
            path = os.path.join(repo, folder, os.path.basename(token))
            if not os.path.exists(path):
                missing.append(path)
    return missing


def _wrap(label, text, width=78):
    return textwrap.fill(text, width, initial_indent=f"    {label:<9s}",
                         subsequent_indent=" " * 13)


def run(report=print):
    web = run_chain()
    verdicts = candidate_verdicts()
    nodes = [n for n, led in web.items() if isinstance(led, root.Ledger)]
    placed = [n for r in RUNGS for n in r["nodes"]]
    unplaced_nodes = sorted(set(nodes) - set(placed))
    doubly = sorted({n for n in placed if placed.count(n) > 1})
    phantom = sorted(set(placed) - set(nodes))
    edges, unplaced_channels = edges_of(web)
    record = coupling_record(edges)
    feeders = sources_by_node(edges)
    cm_rows, sigma = common_mode(web)
    missing = sorted({c for r in RUNGS for c in r["counter"]} - set(verdicts))

    report(f"EMERGENCE LADDER: {len(RUNGS)} rungs bottom up, couplings read "
           f"off root.WEB ({len(nodes)} nodes)")
    for i, r in enumerate(RUNGS):
        report(f"\n{i:2d} {r['name']}")
        report(_wrap("closes", r["closes"]))
        report(_wrap("keeps", r["keeps"]))
        report(_wrap("spins", r["spins"]))
        if r["read"]:
            report(_wrap("read", r["read"]))
        for n in r["nodes"]:
            led = web[n]
            unit = f" {led.unit}" if led.unit else ""
            report(f"    node     {n} = {led.value():.10g}{unit} "
                   f"[web depth {led.web_depth}]")
        for e in (e for e in edges if e["tgt"] == r["key"]):
            depth = "ratio" if e["depth"] is None else f"d{e['depth']} {e['kind']}"
            size = "" if e["factor"] is None else f"{e['factor']:+.3e}"
            report(f"      {depth:9s} {','.join(e['srcs']):15s} "
                   f"{e['path'][:30]:30s} {size:>11s} -> {e['node']}")
        if r["key"] == "gravity":                         # the lane's gates: verdict, tier, witnesses
            for name, g in GATES.items():
                report(f"    gate     {name:24s} {g['status']} {g['tier']:8s} {len(g['witness'])} witnesses")
            readers = sorted({e["node"] for e in edges if "gravity" in e["srcs"] and e["tgt"] != "gravity"})
            report(_wrap("readers", ", ".join(readers) or "none: no node of another rung reads the common mode"))
        if r["counter"]:
            cf = ", ".join(f"{c} {'/'.join(verdicts.get(c, ('?', '?')))}"
                           for c in r["counter"])
            report(_wrap("counter", cf))

    cols = [k for k in KEYS if any(t == k for (_, t) in record)]
    rows = [k for k in KEYS if any(s == k for (s, _) in record)]
    report("\nCOUPLING RECORD (row feeds column; digits are web depths, "
           "r a ratio read)")
    report("    " + " " * 10 + "".join(f"{c[:8]:>9s}" for c in cols))
    for s in rows:
        cells = ["".join(record.get((s, t), ["."])) for t in cols]
        report(f"    {s[:10]:10s}" + "".join(f"{c:>9s}" for c in cells))

    report("\nCOMMON MODE (Sigma_vent(UV) rows by rung)")
    for label, rungs, w in cm_rows:
        report(f"    {label:14s} {', '.join(rungs):16s} {w:+9.4f}")
    total = sum(w for *_, w in cm_rows)
    no_bridge = sum(w for _, rungs, w in cm_rows if "bridge" not in rungs)
    report(f"    sum = {total:.4f} = Sigma {sigma:.4f}. Without the bridge "
           f"{no_bridge:+.4f}")

    widest = max(((n, s) for n, s in feeders.items()
                  if NODE_RUNG.get(n) != "gravity"),
                 key=lambda ns: len(ns[1]), default=(None, set()))
    sources = missing_sources()
    report(f"\nFEEDERS: Sigma_vent(UV) <- {sorted(feeders['Sigma_vent(UV)'])}; "
           f"widest elsewhere: {widest[0]} <- {sorted(widest[1])}")
    report(f"PLACEMENT: unplaced nodes {unplaced_nodes}, unplaced channels "
           f"{unplaced_channels}, doubly placed {doubly}, phantom {phantom}, "
           f"unknown candidates {missing}, missing sources {sources}")
    return dict(rungs=RUNGS, nodes=nodes, edges=edges, record=record,
                feeders=feeders, common_mode=cm_rows, sigma=sigma,
                sigma_sum=total, sigma_without_bridge=no_bridge,
                unplaced_nodes=unplaced_nodes,
                unplaced_channels=unplaced_channels, doubly_placed=doubly,
                phantom_nodes=phantom, missing_candidates=missing,
                missing_sources=sources,
                gates={name: (g["status"], g["tier"], len(g["witness"])) for name, g in GATES.items()},
                gravity_readers=sorted({e["node"] for e in edges
                                        if "gravity" in e["srcs"] and e["tgt"] != "gravity"}),
                verdicts={c: verdicts.get(c) for r in RUNGS
                          for c in r["counter"]})


if __name__ == "__main__":
    run()
