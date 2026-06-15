# Changelog

Every change to a canonical value or its presentation, dated. Papers
are dated snapshots; this file records where the living canon has moved
since the paper versions (divergences are folded back into the papers
at their next Zenodo version).

## 2026-06-14

- Depth-2 echoes: `sin²θ_W`, `m_s`, `M_W`, `M_Z`, `m_H` updated.
  Watches 3 and 4 superseded with new values.
- Three new modules: `embedding_uniqueness.py` (six-gate proof),
  `protected_forgetting.py` (PvP=0, Pv²P=½P), `nls_soliton.py`
  (Z₃-coupled NLS). Suite: 59 → 62 tests.
- `null_model.py` removed (rhetorical, not structural).
- README rewritten end-to-end: corrected all stale numbers, back-reaction
  terminology, 62 tests, added v_EW at full precision.

## 2026-06-15

- Collaboration structure added: `CONTRIBUTING.md`,
  `watches/REGISTRY.md`, `PAPERS.md`, `CHANGELOG.md`.

## 2026-06-12

- Mass-coordinate doctrine finalized in `masses.py` and `__init__.py`:
  one rule over three dynamical classes; the transport statement is a
  coordinate change on the reader's side, not a step in the prediction.
- Scheme-invariant ratios added as canonical predictions with freeze
  tests: `m_u/m_d = 38/83` (−0.9σ), `m_s/m_ud = 27.318` (+0.2σ,
  registered watch), `Q_ellipse = 22.383`. Suite: 58 → 59 tests.

## 2026-06-11

- α(0) presentation completed: the canonical value is the depth-3
  cubic root `137.035999050` (in the code since 2026-06-10, commit
  062593c); the depth-1 truncation `137.036439` is an intermediate
  identity and is now labeled as such everywhere it appears.
- Gauge-matching scale correction (the μ\* migration): the exact WZW
  cancellation terminates at `μ* = M_Pl·e^−(9π²/2−6) ≈ 253.5 GeV`,
  not at `v_EW`. Earlier releases mislabeled the scale; the conflation
  entered when the 30-mode vertex echo was added to the v_EW exponent.
  Consequences: `α_s(M_Z)` 0.1177 → 0.1184 (+0.43σ), `m_H` → 124.06.
  The pre-migration registration of the 4/13 watch is preserved in
  `interference/__init__.py` with the degeneracy documented.
- Face-split law promoted to canonical: the bridge self-echo splits
  4/13 (G₂) / 9/13 (F₄) in the h∨ metric; `G_ind/G_N = 0.999999917`;
  predicts ABSENCE of anomalous gauge dressing (watch 1).
- Electroweak sector closed at one loop with two declared, recalculated
  PDG 2024 imports: `M_W` −0.33σ, `M_Z` +0.0088%.

## Earlier

- 2026-06-10 (062593c): depth-3 e↔q echo enters the α(0) chain and the
  S_quark ledger (back-reaction, forced multiplicity 16 and cycle
  orientation 2). v_EW closes to G_F at 0.02 ppm.
- 2026-05-23 (054004c): unified polynomial version; depth-1 α(0) with
  its 0.0003% deviation declared in the output.
- 2026-05-21: v0.1 release (per-layer E8 package, M_Pl anchor).
