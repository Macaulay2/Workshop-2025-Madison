# Strategy Comparison — Benchmark Results

Measured by `benchmarks/bench_problems.m2` with `setRandomSeed 42` on branch
`EliminationTemplates-unified` (post-Increment-3, commit `44173f8`). Template
sizes are `rows × cols`. `build+solve` is `cpuTime` in seconds. `res` is
`max |F_j(x*)|` over recovered solutions (QQ only; skipped for ZZ/p).

Supporting analysis and open TODOs are in `matrixhi_scaling_gap.md`.

## Strategy labels

Under the unified Direction A pipeline, every strategy flows through
`getTemplate` → `[E | R | B]` fold → `extractActionFromTemplate` (RREF +
pivot). The only per-strategy branch is inside `getH0`:

| Label | Option combination | `H0` body |
|---|---|---|
| Default | _no options_ (`Strategy => null`) | Gröbner change-of-basis |
| Larsson | `Strategy => "Larsson"` | Default's H0, then mod the first syzygy module of F |
| MatrixHi-mode | `Strategy => "Greedy", AdjustParams => false` | Martyushev gap polys + `buildHSymbolic`, α := 0 (particular solution) |
| Greedy | `Strategy => "Greedy"` (default `AdjustParams => true`) | Martyushev gap polys + `buildHSymbolic` + `adjustParams`; QQ input pinned to ZZ/p internally (Increment 2) |

There is no distinct "MatrixHi" strategy identifier anymore; the paper's
MatrixHi particular solution is reached by the `AdjustParams => false`
flag on Greedy. Polynomial actions are lifted to `R_s = R[s] / ⟨s − a⟩`
automatically (Increment 3) so every strategy accepts both monomial and
polynomial action variables.

## Summary table

| Problem | Ring | deg I | Default | Larsson | MatrixHi-mode | Greedy | Paper std |
|---|---|---:|---|---|---|---|---|
| circle+line (toy) | QQ | 2 | 4×6 ✗ | 4×6 ✗ | 4×6 ✗ | 4×6 ✗ | — |
| two conics | QQ | 4 | 4×9 ✗ | 4×9 ✗ | 4×9 ✗ | 4×9 ✗ | — |
| paper §3 ex (2-var) | QQ | 12 | **7×19** ✓ | 7×19 ✓ | 7×19 ✓ | 7×19 ✓ | — |
| 3-var docs | QQ | 12 | **16×28** ✓ | 16×28 ✓ | 20×33 ✓ | 17×29 ✓ | — |
| 5pt essential +det | QQ | 10 | 10×20 ✓ | 10×20 ✓ | 10×20 ✓ | **10×20** ✓ | 10×20 |
| 6 Demazure cubics | QQ | 10 | 27×34 ✓ | **24×34** ✓ | 25×35 ✓ | **24×34** ✓ | — |
| #1 F+λ 8pt | ZZ/p | 8 | 11×20 | 11×20 | 11×20 | **11×20** | 11×19 |
| #2 E+f 6pt | ZZ/p | 9 | 63×45 | 33×42 | 52×57 | **21×30** | 11×20 |
| #3 f+E+f 6pt | ZZ/p | 15 | 200×120 | 53×73 | 93×99 | **31×50** | 12×27 |
| 3-var mixed | QQ | 6 | **5×11** ✓ | 5×11 ✓ | 5×11 ✓ | 5×11 ✓ | — |

Legend:
- **bold** = strategy-specific best for the problem (ties broken by dominance).
- ✓ = tier-3 residual `< 1e-6`.
- ✗ = pre-existing recovery failure (every strategy returns same partial
  solution); neither introduced nor fixed by Increments 1–3.

Key observations:

- Monomial-action problems where free α = 0 (paper §3, 5pt+det, F+λ 8pt,
  3-var mixed): all four strategies reach the same size; Direction A
  lifts all of them to paper-reference sizes.
- Free-α problems where `adjustParams` commits parameters (6 Demazure,
  #2 E+f, #3 f+E+f): Greedy strictly wins. The gap between Greedy and
  MatrixHi-mode on #2/#3 comes from `adjustParams` cancelling
  redundant shift support.
- 3-var docs: Default's Gröbner H0 has smaller monomial support than
  `buildHSymbolic`'s — one of the rare cases where Default beats Greedy.

## Timing (seed 42, cpuTime)

| Problem | Default | Larsson | MatrixHi-mode | Greedy |
|---|---:|---:|---:|---:|
| paper §3 ex | 0.005 s | 0.006 s | 0.058 s | 0.078 s |
| 3-var docs | 0.047 s | 0.007 s | 0.169 s | 0.304 s |
| 5pt essential +det | 0.019 s | 0.023 s | 0.148 s | 0.245 s |
| 6 Demazure cubics | 0.025 s | 0.078 s | 0.453 s | **40.1 s** |
| #1 F+λ 8pt | 0.005 s | 0.051 s | 0.122 s | 0.173 s |
| #2 E+f 6pt | 0.018 s | 0.062 s | 0.973 s | **45.9 s** |
| #3 f+E+f 6pt | 0.076 s | 0.023 s | 2.752 s | **149.3 s** |
| 3-var mixed | 0.007 s | 0.007 s | 0.017 s | 0.070 s |

Greedy's cost is dominated by `adjustParams`; 6 Demazure / #2 / #3 spend
most of their time in the `adjustParams` outer loop walking candidate α
commits. ZZ/p pinning inside `adjustParams` (Increment 2) keeps the inner
RREF from blowing up on QQ input but doesn't change the per-iteration
work.

## Polynomial action (Increment 3)

Under Direction A's unified pipeline, polynomial action variables like
`x + 4y` or `random(1, R)` now flow through the same `[E | R | B]` +
RREF + pivot path as monomial actions — the action is lifted to
`s` in `R_s = R[s] / ⟨s − a⟩`, where `s` is itself a ring variable on
the quotient.

| Problem | Action | Template | Residual |
|---|---|---|---|
| `change of ideals` (QQ[x,y], `x²+y²−1`, `x²+y³+xy−2`) | `x + 4y` | 20×26 | ~1e-14 |
| 5pt essential (Demazure cubics, no det) | `random(1, R)` | 37×34 | ~1e-12 (one outlier ~1e-8) |

Both paths still run RREF over QQ natively. A future LU split-solve
fallback (already wired in `getActionMatrix(EliminationTemplate)` via
`try / else`) kicks in if RREF raises a catchable error on adversarial
QQ input.

## Where Greedy strictly shrinks MatrixHi-mode (free α > 0 regime)

| Problem | MatrixHi-mode | Greedy | Savings |
|---|---|---|---|
| 6 Demazure cubics | 25×35 | 24×34 | 1 row, 1 col |
| #2 E+f 6pt | 52×57 | 21×30 | 31 rows, 27 cols |
| #3 f+E+f 6pt | 93×99 | 31×50 | 62 rows, 49 cols |

On every other problem in the table, free α = 0 and Greedy = MatrixHi-mode.

## Reproducing

```bash
cd Workshop-2025-Madison

# Single-problem bench (5pt essential, Larsson / MatrixHi-mode / Greedy)
M2 benchmarks/bench_5pt_all.m2

# adjustParams shrinking demo (6 Demazure cubics, no det)
M2 benchmarks/bench_greedy_effect.m2

# Full 10-problem matrix (source of this document)
M2 benchmarks/bench_problems.m2

# Package regression tests (20 TEST blocks)
M2 -e 'needsPackage "EliminationTemplates"; check "EliminationTemplates"; exit 0'

# Monomial-action size probe (sanity check on paper-reference sizes)
M2 playground/notes/direction_a_size_probe.m2
```
