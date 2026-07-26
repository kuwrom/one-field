"""
test_blindness.py: the reference-data firewall, executable.

THE CLAIM.  The derivation chain touches measured reality through a
declared input surface and nowhere else:

  ANCHOR      m_e (CODATA)                root.M_E_ANCHOR_MEV
  Λ-SECTOR    H₀                          root.H_0_KM_S_MPC
  ν-SECTOR    Δm²₂₁, Δm²₃₁ (NuFit)        root.NUFIT_OSCILLATION
              (recorded reference; no interference/ code reads it:
              the ν block is structural, no Δm² input)
  EW IMPORTS  Δr̂_W, ρ̂ (PDG Eq. 10.26)    root.PDG_EW['dr_hat_W','rho_hat']
  COORDINATE  M_Z (evaluation scale for   root.PDG_EW['M_Z']
              the α_s comparison; shown
              immaterial below)

Everything else in root.{PDG_MASSES, PDG_CKM, NUFIT_PMNS, PDG_EW,
PDG_COSMO} is comparison-side: it may appear in pulls, χ², and error
columns, never inside a prediction.

THE ENFORCEMENT, three layers:

  1. test_reference_quarantine_fast: perturb EVERY comparison-side
     reference value simultaneously, re-run root → masses → mixing →
     couplings, assert every prediction bit-identical to canon.  This
     extends the nine-mass blind holdout (test_coverage.py) to the
     whole fast chain.
  2. test_reference_quarantine_gravity (slow): same, through gravity
     and the dark sector (the 20 000-step RGE makes it slow).
  3. test_reference_access_surface_frozen: an AST audit freezing the
     complete set of source locations that read a reference dict.  A
     new read anywhere in interference/ fails CI and must be admitted
     here consciously, with its classification.

Plus test_alpha_s_evaluation_scale_immaterial: the one derivation-side
read that is a coordinate rather than a steering input (couplings
evaluates the RGE at the scale where PDG quotes α_s) is shown to move
the prediction by < 3e-5 even at a 5σ perturbation of M_Z, i.e. far
inside the declared band.
"""

import ast
import contextlib
import io
import os
import pathlib

import pytest

import root
import masses
import mixing
import couplings
import gravity
import dark_sector


# ── the declared input surface (unperturbed by the quarantine) ─────────
_DECLARED = {
    ("PDG_EW", "dr_hat_W"),   # declared radiative import
    ("PDG_EW", "rho_hat"),    # declared radiative import
    ("PDG_EW", "M_Z"),        # evaluation coordinate, immaterial (below)
}
# NUFIT_OSCILLATION is a recorded reference; no interference/ code
# consumes it (the ν block is structural, no Δm² input).  Left
# unperturbed.
_REF_NAMES = ("PDG_MASSES", "PDG_CKM", "NUFIT_PMNS", "PDG_EW",
              "PDG_COSMO", "PLANCK_2018", "ACT_DR6")
# PLANCK_2018 / ACT_DR6 live in reference.py and sit inside the
# perturbation set: the DM/baryon prediction must not move when they
# do.


def _perturb_all():
    """Scale every comparison-side reference value (alternating
    ×1.2 / ×0.8), returning a saved snapshot for restore."""
    saved = {}
    i = 0
    for name in _REF_NAMES:
        d = getattr(root, name)
        saved[name] = dict(d)
        for k in d:
            if (name, k) in _DECLARED:
                continue
            f = 1.2 if i % 2 == 0 else 0.8
            i += 1
            v = d[k]
            d[k] = tuple(x * f for x in v) if isinstance(v, tuple) else v * f
    return saved


def _restore(saved):
    for name, snap in saved.items():
        getattr(root, name).clear()
        getattr(root, name).update(snap)


def _rerun(through_gravity=False):
    with contextlib.redirect_stdout(io.StringIO()):
        R = root.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
        c = couplings.derive(R, m)
        if not through_gravity:
            return R, m, mx, c, None, None
        g = gravity.derive(R, m, mx, c)
        d = dark_sector.derive(R, g)
    return R, m, mx, c, g, d


# prediction keys asserted bit-identical under full perturbation
_MASS_KEYS = ("m_e", "m_mu", "m_tau", "m_u", "m_d", "m_s",
              "m_c", "m_b", "m_t", "v_EW_pred_GeV",
              "mu_over_md", "ms_over_mud", "Q_ellipse")
_MIX_KEYS = ("lam", "A", "eta", "rho", "sin2_12", "sin2_13",
             "delta_CP_deg", "J_ckm")
_COUP_KEYS = ("alpha_s_MZ_thresh", "alpha_s_mu_star", "inv_alpha_phys",
              "mu_star_GeV")
_GRAV_KEYS = ("G_ratio_UV", "G_ratio_broken", "lambda_mt_transport", "sin2W",
              "M_W_derived", "M_Z_derived", "eta_B", "rho_Lambda")


def test_reference_quarantine_fast(res):
    """Perturb ALL comparison-side reference data at once; every
    prediction of the fast chain must not move by one part in 1e12."""
    saved = _perturb_all()
    try:
        _, m2, mx2, c2, _, _ = _rerun()
        for k in _MASS_KEYS:
            assert m2[k] == pytest.approx(res["m"][k], rel=1e-12), \
                f"masses[{k}] moved: a reference value steers a prediction"
        for k in _MIX_KEYS:
            assert mx2[k] == pytest.approx(res["mx"][k], rel=1e-12), \
                f"mixing[{k}] moved: a reference value steers a prediction"
        for k in _COUP_KEYS:
            assert c2[k] == pytest.approx(res["c"][k], rel=1e-12), \
                f"couplings[{k}] moved: a reference value steers a prediction"
    finally:
        _restore(saved)
        _rerun()   # restore global web state for later consumers


@pytest.mark.slow
def test_reference_quarantine_gravity(res):
    """Same quarantine, through gravity and the dark sector."""
    saved = _perturb_all()
    try:
        _, _, _, _, g2, d2 = _rerun(through_gravity=True)
        for k in _GRAV_KEYS:
            assert g2[k] == pytest.approx(res["g"][k], rel=1e-12), \
                f"gravity[{k}] moved: a reference value steers a prediction"
        assert d2["DM_baryon_ratio"] == pytest.approx(
            res["d"]["DM_baryon_ratio"], rel=1e-12)
    finally:
        _restore(saved)
        # cheap restore: on pass, the gravity web nodes are provably
        # canon already (that is what the test asserts); the fast
        # chain repairs everything upstream.
        _rerun()


def test_alpha_s_evaluation_scale_immaterial(res):
    """The RGE endpoint M_Z is a coordinate (where PDG quotes α_s),
    not a steering input: a 5σ shift of M_Z moves α_s(M_Z) by less
    than 3e-5, two orders below the declared 0.5% band."""
    saved = dict(root.PDG_EW)
    try:
        root.PDG_EW["M_Z"] = saved["M_Z"] + 5 * 0.0021   # +5σ_exp
        with contextlib.redirect_stdout(io.StringIO()):
            c2 = couplings.derive(res["R"], res["m"])
        shift = abs(c2["alpha_s_MZ_thresh"] - res["c"]["alpha_s_MZ_thresh"])
        assert shift < 3e-5, f"α_s endpoint sensitivity too large: {shift}"
        assert shift > 0.0   # the coordinate is really being read
    finally:
        root.PDG_EW.update(saved)
        with contextlib.redirect_stdout(io.StringIO()):
            couplings.derive(res["R"], res["m"])


# audit stage 3: the frozen access surface
#
# Complete inventory of reference-dict reads in interference/*.py, per
# module.  '[*]' = non-constant subscript (iteration / display loops).
# root.py entries are the definitions themselves plus derive()'s
# pass-through packaging.
_FROZEN_SURFACE = {
    "couplings.py": [
        "PDG_EW", "PDG_EW['M_Z']", "PDG_EW['alpha_s_MZ']"],
    # dark_sector: comparison-side density fits, single-homed in
    # reference.py
    "dark_sector.py": [
        "ACT_DR6", "ACT_DR6['Omega_b_h2']", "ACT_DR6['Omega_b_h2_err']",
        "ACT_DR6['Omega_c_h2']", "ACT_DR6['Omega_c_h2_err']",
        "PDG_COSMO", "PDG_COSMO['Omega_Lambda']",
        "PLANCK_2018", "PLANCK_2018['Omega_b_h2']",
        "PLANCK_2018['Omega_b_h2_err']", "PLANCK_2018['Omega_c_h2']",
        "PLANCK_2018['Omega_c_h2_err']"],
    "gravity.py": [
        "PDG_COSMO", "PDG_COSMO['eta_B']"],
    "higgs_ew.py": [
        "PDG_EW", "PDG_EW['M_W']", "PDG_EW['M_Z']",
        "PDG_EW['dr_hat_W']", "PDG_EW['rho_hat']",
        "PDG_EW['sin2W_MSbar']", "PDG_EW['sin2W_err']"],
    "masses.py": [
        "PDG_MASSES", "PDG_MASSES['*']", "PDG_MASSES['b']",
        "PDG_MASSES['c']", "PDG_MASSES['d']", "PDG_MASSES['s']",
        "PDG_MASSES['t']", "PDG_MASSES['u']"],
    "mixing.py": [
        "NUFIT_PMNS", "NUFIT_PMNS['*']", "NUFIT_PMNS['sin2_12']",
        "NUFIT_PMNS['sin2_13']", "NUFIT_PMNS['sin2_23']",
        "PDG_CKM", "PDG_CKM['*']", "PDG_CKM['A']", "PDG_CKM['J']",
        "PDG_CKM['etabar']", "PDG_CKM['gamma_deg']",
        "PDG_CKM['lambda']", "PDG_CKM['rhobar']"],
    # reference.py: the definitions themselves (the single home;
    # root.py re-exports and does not appear on the surface)
    "reference.py": [
        "ACT_DR6", "NUFIT_OSCILLATION", "NUFIT_PMNS", "PDG_CKM",
        "PDG_COSMO", "PDG_EW", "PDG_MASSES", "PLANCK_2018",
        "SUPER_K"],
    "run.py": [
        "PDG_MASSES", "PDG_MASSES['*']"],
    "words.py": [
        "PDG_MASSES", "PDG_MASSES['tau']"],
}

_ALL_REF = set(_REF_NAMES) | {"NUFIT_OSCILLATION", "SUPER_K"}
# (PLANCK_2018 / ACT_DR6 are in _REF_NAMES, hence scanned too;
# SUPER_K is scan-only: probes read it, the engine does not)

# ── the probe-side frozen surface ──────────────────────────────────────
# Probes are cited as executable witnesses, so measured data can
# reach a claim through tests/probes/.  Probe reference values are
# routed through reference.py and the probe directory is scanned and
# frozen exactly like interference/.  (Literals the scanner cannot
# see remain possible; the routed dicts plus this freeze cover the
# known cases.)
_FROZEN_SURFACE_PROBES = {
    # condensate_map: comparison-side only.  The alpha_s(M_Z) pulls
    # of the condensate-map contraction dichotomy and scale menu; no
    # derived number consumes the reads.
    "condensate_map.py": [
        "PDG_EW", "PDG_EW['M_Z']", "PDG_EW['alpha_s_MZ']"],
    "conversion_rigidity.py": ["PDG_MASSES", "PDG_MASSES['tau']"],
    "edge_derivations.py": [
        "PDG_MASSES", "PDG_MASSES['d']", "PDG_MASSES['s']"],
    "ew_internal.py": [
        "PDG_EW", "PDG_EW['dr_hat_W']", "PDG_EW['rho_hat']"],
    "lambda_assignment.py": ["PDG_MASSES", "PDG_MASSES['b']"],
    "substrate_depth.py": [
        "SUPER_K", "SUPER_K['N_baryon']", "SUPER_K['t_Hubble_s']",
        "SUPER_K['tau_p_yr']", "SUPER_K['yr_s']"],
    "surplus_edges.py": [
        "PDG_MASSES", "PDG_MASSES['e']", "PDG_MASSES['mu']",
        "PDG_MASSES['tau']"],
    "theta_menu.py": [
        "PDG_MASSES", "PDG_MASSES['e']", "PDG_MASSES['mu']",
        "PDG_MASSES['tau']"],
}

# The same AST freeze over tests/audit/.  closure_search is the forward
# search's DATA-CONTACT stage: the one place the generator meets
# measurement, applied once, last, by design (see its docstring).
# Comparison-side by construction.  It reads reference data, so it is
# part of the epistemic surface wherever it lives, so it is frozen
# with tests/audit/ rather than with interference/.
# PDG_MASSES['d']: the report-only down-gen-1 stage, same
# classification as the other closure-search stages.
_FROZEN_SURFACE_AUDIT = {
    "closure_search.py": [
        "PDG_MASSES", "PDG_MASSES['b']", "PDG_MASSES['c']",
        "PDG_MASSES['d']", "PDG_MASSES['e']", "PDG_MASSES['mu']",
        "PDG_MASSES['s']", "PDG_MASSES['t']", "PDG_MASSES['tau']",
        "PDG_MASSES['u']"],
}


def _ref_name(node):
    """The reference-dict name this node reads, or None.

    Covers the bare name (`PDG_MASSES`, the repo convention) and the
    attribute form (`reference.PDG_MASSES`, `root.PDG_EW`), which is
    the same data contact by another spelling."""
    if isinstance(node, ast.Name) and node.id in _ALL_REF:
        return node.id
    if isinstance(node, ast.Attribute) and node.attr in _ALL_REF:
        return node.attr
    return None


def _scan_source(src):
    """Every reference-dict read in one module's source."""
    hits = set()
    for node in ast.walk(ast.parse(src)):
        # `PDG_MASSES[...]` and `reference.PDG_MASSES[...]` both count:
        # an attribute read is the same data contact, and scanning only
        # bare names would leave `import reference` as a silent way
        # through the firewall.
        target = _ref_name(node.value) if isinstance(
            node, ast.Subscript) else None
        if target:
            key = (node.slice.value
                   if isinstance(node.slice, ast.Constant) else "*")
            hits.add(f"{target}[{key!r}]")
        elif _ref_name(node):
            hits.add(_ref_name(node))
    return hits


def _scan_surface(subdir="interference"):
    pkg = pathlib.Path(os.path.join(
        os.path.dirname(__file__), "..", subdir))
    surface = {}
    for p in sorted(q for q in pkg.rglob("*.py")
                    if "__pycache__" not in q.parts):
        hits = _scan_source(p.read_text())
        if hits:
            surface[p.name] = sorted(hits)
    return surface


def test_scanner_sees_both_spellings():
    """The scanner itself, tested.  No module in the repo currently
    uses the attribute form, so nothing else would notice if this
    branch were wrong: the freeze would pass while the hole stayed
    open.  Both spellings must land on the same entry."""
    assert _scan_source("x = PDG_MASSES['t']") == {
        "PDG_MASSES", "PDG_MASSES['t']"}
    assert _scan_source("import reference\nx = reference.PDG_MASSES['t']") \
        == {"PDG_MASSES", "PDG_MASSES['t']"}
    assert _scan_source("x = obj.mass") == set()


def test_reference_access_surface_frozen():
    """The set of source locations reading measured reference data is
    FROZEN.  Adding a reference read anywhere in interference/ fails
    here and must be admitted consciously (edit _FROZEN_SURFACE with
    its classification: comparison-side, declared input, or
    coordinate)."""
    assert _scan_surface() == _FROZEN_SURFACE


def test_probe_reference_surface_frozen():
    """The same AST freeze over tests/probes/: probes are cited as
    witnesses, so this directory is part of the epistemic
    surface too.  A new reference read (or a hardcoded value) must
    be admitted here."""
    assert _scan_surface(os.path.join("tests", "probes")) \
        == _FROZEN_SURFACE_PROBES


def test_audit_reference_surface_frozen():
    """The same AST freeze over tests/audit/.  The closure search is
    the look-elsewhere instrument, and it makes data contact by
    design, so its reference reads stay frozen."""
    assert _scan_surface(os.path.join("tests", "audit")) \
        == _FROZEN_SURFACE_AUDIT
