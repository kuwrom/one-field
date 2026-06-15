# Contributing

The framework has no adjustable parameters, so contributions follow a
few rules. They are stated below.
 
## Claim statuses

Every claim in the repo carries one of these statuses (matching
the labels used in the code itself):

- **FORCED**: a back-reaction term whose multiplicity is fixed by the
  algebra. This is the default status of every `EchoTerm` in the
  `Ledger`/`Web` machinery (`root.py`). A forced echo has zero free
  choices: the weight, depth, and combinatoric factor are all derived
  from the four integers, Casimirs, or central charges.
- **DERIVED**: a closed-form consequence of the four integers
  `(d₁₀, d₁₁, n₇, n₂₆) = (2, 3, 7, 26)` plus π, with the derivation in
  the module docstring and the value frozen in the test suite.
- **ENUMERATED**: selected by an exhaustive scan over a declared
  candidate menu, with the losing candidates listed (see the selection
  checks in `masses.py` and `words.py`).
- **DECLARED IMPORT**: an external input stated with its source and
  uncertainty band (currently two: the PDG 2024 one-loop `Δr̂_W` and
  `ρ̂`), each decomposed and recalculated for compatibility in
  `gravity.py`.

A contribution states which status it claims.

## Ways to contribute

**Verification.** Clone, run `pytest -q` and `python interference/run.py`,
and report your environment and results as an issue.

**Critique.** File an issue. A wrong number, a mislabeled comparison,
or a derivation step that does not follow are all useful reports.

**Tooling.** Scheme-transport calculators, higher-loop machinery
(an SMDR-class NNLO Higgs chain is the standing wish), CI, packaging.
Tooling does not touch canonical values.

**Candidate derivations.** Follow the review process below.

## Review process

1. **Issue.** Open an issue describing the target observable, the
   dictionary objects you plan to use, and the forcing argument.
   Discuss the approach with existing contributors before computing
   anything.
2. **Green light.** At least two contributors signal approval in
   the issue discussion before work proceeds.
3. **Pull request.** Submit the derivation as a PR. The review
   checks that:
   - The result is a closed-form consequence of the existing
     dictionary, with zero new parameters and zero free choices.
   - The trial count is declared: how many candidates could have
     been tried, and the resulting false-positive expectation.
   - The freeze test extends: the new value enters
     `tests/test_interference.py` at full precision.
   - The change is logged: one dated line in `CHANGELOG.md`.
4. **Merge.** At least two reviewers approve the pull request.

Contributors are welcome to list any associated papers in
[`PAPERS.md`](PAPERS.md) as part of their pull request.

## Out of scope

- Corrections containing anything adjustable.
- Fits.
- Changes to canonical values without a reviewed derivation.

## Style

No em dashes. Plain ASCII in code comments where possible. Every
number printed by the code carries its comparison value and pull.
Docstrings state derivations completely or point to the module that
does.
