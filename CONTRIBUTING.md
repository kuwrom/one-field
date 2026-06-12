# Contributing

The framework has no adjustable parameters, so contributions follow a
few rules. They are stated below.
 
## Claim statuses

Every claim in the repo carries exactly one of these statuses:

- **DERIVED**: a closed-form consequence of the four integers
  `(d₁₀, d₁₁, n₇, n₂₆) = (2, 3, 7, 26)` plus π, with the derivation in
  the module docstring and the value frozen in the test suite.
- **ENUMERATED**: selected by an exhaustive, pre-registered scan over a
  declared candidate menu, with the losing candidates listed (see the
  selection checks in `masses.py` and `words.py`).
- **CROSS-PINNED**: the same dictionary object measured in independent
  observables (e.g. `h₁₀ = 2/9` appears in the Cabibbo angle, the
  Koide phase, and the `m_u` correction).
- **CITED IMPORT**: a declared external input (currently two: the
  PDG 2024 one-loop `Δr̂_W` and `ρ̂`), each decomposed and recalculated
  for compatibility in `gravity.py`.

A contribution states which status it claims.

## Ways to contribute

**Verification.** Clone, run `pytest -q` and `python interference/run.py`,
and report your environment and results in `verification/`.

**Critique.** File an issue. A wrong number, a mislabeled comparison,
or a derivation step that does not follow are all useful reports.

**Tooling.** Scheme-transport calculators, higher-loop machinery
(an SMDR-class NNLO Higgs chain is the standing wish), CI, packaging.
Tooling does not touch canonical values.

**Candidate derivations.** See the lifecycle below.

## The candidate lifecycle

New physics goes through `candidates/`:

1. Copy `candidates/_TEMPLATE/` to `candidates/<your-topic>/`.
2. Fill in `CLAIM.md` before computing anything: the target
   observable, the candidate dictionary menu, the promotion bar, and
   the kill condition.
3. Put the computation in `derivation.py` (standalone, no dependencies
   beyond the `interference/` package and numpy).
4. Record the outcome in `STATUS.md`: `CANDIDATE`, `KILLED`, or
   `PROMOTED`, with the date.

**The promotion bar** (all four required):

- The result is a closed-form consequence of the existing dictionary,
  with zero new parameters and zero choices the dictionary does not
  already fix.
- The trial count is declared: how many candidates could have been
  tried, and the resulting false-positive expectation.
- The freeze test extends: the new value enters
  `tests/test_interference.py` at full precision.
- The change is logged: one dated line in `CHANGELOG.md`.

Promotion is gated by the maintainer. Killed candidates move to
`killed/` with their post-mortem and stay there.

## Out of scope

- Corrections containing anything adjustable.
- Fits.
- Changes to canonical values without a derivation meeting the
  promotion bar.
- Hand edits to the README scorecard (run `python tools/scorecard.py
  --write`; CI checks this).

## Style

No em dashes. Plain ASCII in code comments where possible. Every
number printed by the code carries its comparison value and pull.
Docstrings state derivations completely or point to the module that
does.
