# Registered watches 

Every live bet the framework has placed ahead of the data, with its
kill condition and the experiment that decides it. A watch is not a
claim; it is a falsifiable commitment registered before the deciding
measurement exists. Source of truth for the values: the code
(`interference/`), frozen in `tests/test_interference.py`.

| # | Watch | The bet | Falsified by | Decided by | Registered |
|---|-------|---------|--------------|------------|------------|
| 1 | 4/13 absence (α_s) | `α_s(M_Z) = 0.1184` with NO anomalous gauge dressing; the G₂-absorbed share of the bridge self-echo is invisible to gauge observables (no-self-dilution symmetry) | `α_s(M_Z) = 0.1177` or `0.1191` (a ±(4/13)/(2π) shift) | lattice α_s at ~0.03% | 2026-06 (pre-μ\*-migration registration quoted 0.1177 under the v_EW-start convention; see `interference/__init__.py` for the documented degeneracy and why the matching scale is fixed first) |
| 2 | α(0) cesium commitment | `1/α(0) = 137.035999050` (depth-3 cubic), +0.13σ of Berkeley Cs | the Cs/Rb discrepancy (5.5σ) resolving on the Rb side | next-generation recoil measurement | 2026-06 |
| 3 | m_s/m_ud ratio | `m_s/m_ud = 27.318` (RG-invariant, no coordinate) | lattice average tightening at 27.30(8) and staying | FLAG-class lattice ratio | 2026-06 (supersedes 27.130; depth-2 echo update 2026-06-14) |
| 4 | sin²θ_W emission weight | `sin²θ_W(M_Z) = 0.231285` | the 5/13 weight (0.231216) or any other family weight | FCC-ee Z-pole (σ ~ 9 × 10⁻⁶) | 2026-06 (supersedes 0.231234; depth-2 echo update 2026-06-14) |
| 5 | m_b down-type vent | `m_b = 4193.8 MeV` with NO vent; if a vent exists it is ±(9/4)(α/2π) with ε-orientation signs (+d, −s, +b) | persistence at >3σ without the pre-registered vent value | PDG m_b error shrinking | 2026-06 |
| 6 | Ω_DM/Ω_b | `2π − 1 = 5.2832` (sits between Planck −1.3σ and ACT +1.2σ) | both datasets converging away from 5.283 | CMB-S4 class | 2026-06 |
| 7 | Σ structural set | `m_ν₁ = 0`, normal ordering | inverted ordering, or a measured m_ν₁ > 0 | JUNO / cosmology | 2026-06 |

## Rules

- A watch is registered before the deciding data exists, with its kill
  condition stated at registration time.
- Watches are never silently edited. A superseded watch keeps its
  original registration with a dated note (see watch 1 for the
  pattern).
- When a watch is decided, it moves to a DECIDED section here with the
  outcome, win or lose.

## Decided

(none yet)
