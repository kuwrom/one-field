# Contributing

This is a zero-free-parameter framework. That single fact sets the rules
of the game, because the usual way of improving a physics model
(adjusting something until it fits better) is the one move that does not
exist here. There are no knobs. A proposed correction that contains
anything adjustable is not an improvement, it is a different (and
weaker) kind of theory.

## What every claim in this repo is

Exactly one of:

- **DERIVED**: a closed-form consequence of the four integers
  `(d₁₀, d₁₁, n₇, n₂₆) = (2, 3, 7, 26)` plus π, with the derivation in
  the module docstring and the value frozen in the test suite.
- **ENUMERATED**: selected by an exhaustive, pre-registered scan over a
  declared candidate menu (see the selection checks in `masses.py` and
  `words.py`), with the losing candidates listed.
- **CROSS-PINNED**: the same dictionary object measured in independent
  observables (e.g. `h₁₀ = 2/9` is simultaneously the Cabibbo angle,
  the Koide phase, and the `m_u` correction; it cannot be re-chosen per
  observable).
- **CITED IMPORT**: an explicitly declared external input (currently
  two: the PDG 2024 one-loop `Δr̂_W` and `ρ̂`), each decomposed and
  recalculated for compatibility in `gravity.py`.

If a contribution cannot say which of these it is, it is not ready.

## Ways to contribute

**Verification.** Clone, run `pytest -q` and `python interference/run.py`,
and report your environment and results in `verification/`. Independent
replication is the most valuable contribution there is.

**Critique.** File an issue. Quantitative critique (a number that is
wrong, a comparison that is mislabeled, a derivation step that does not
follow) is treated as a gift regardless of tone. Several of the
sharpest tests in this repo exist because a critic pushed on a soft
spot.

**Tooling.** Scheme-transport calculators, higher-loop machinery
(an SMDR-class NNLO Higgs chain is the standing wish), CI, packaging.
Tooling must not touch canonical values.

**Candidate derivations.** The real game. See the lifecycle below.

## The candidate lifecycle

Nothing enters the canon directly. New physics goes through
`candidates/`:

1. Copy `candidates/_TEMPLATE/` to `candidates/<your-topic>/`.
2. Fill in `CLAIM.md` **before computing anything**: the target
   observable, the candidate dictionary menu, the promotion bar, and
   the kill condition. This pre-registration is what makes a later
   success mean something (and what makes "you adjusted it after the
   fact" impossible to allege).
3. Put the computation in `derivation.py` (standalone, no dependencies
   beyond the `interference/` package and numpy).
4. Record the outcome in `STATUS.md`: `CANDIDATE`, `KILLED`, or
   `PROMOTED`, with the date.

**The promotion bar** (all four required):

- The result is forced: a closed-form consequence of the existing
  dictionary with zero new parameters and zero choices the dictionary
  does not already fix.
- The trial count is declared: how many candidates could have been
  tried, and what the false-positive expectation is (see the
  pre-registered scan pattern in the repo history).
- The freeze test extends: the new value enters
  `tests/test_interference.py` at full precision.
- The change is logged: one dated line in `CHANGELOG.md`.

Promotion is gated by the maintainer. Killed candidates move to
`killed/` with their post-mortem; they are never deleted, because the
kill record is the trial-count ledger that keeps the successes honest.

## What is never merged

- Corrections with adjustable anything (coefficients chosen to fit,
  scales chosen by eye, terms added because they help).
- Fits, even good ones. A fit answers "what value works"; this repo
  only accepts answers to "what value is forced".
- Changes to canonical values without a derivation meeting the
  promotion bar.
- README scorecard edits by hand (run `python tools/scorecard.py
  --write`; CI rejects hand edits).

## Style

No em dashes. Plain ASCII in code comments where possible. Every
number printed by the code carries its comparison value and pull.
Docstrings state derivations completely or point to the module that
does. The repo speaks in declarative sentences and shows its residuals
unprompted.
