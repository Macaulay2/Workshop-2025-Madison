# Strategy Comparison — Benchmark Results

Measured by `workshop/tests/bench_problems.m2` with `setRandomSeed 42`. Template sizes are rows × cols. `build+solve` is `cpuTime` in seconds. "res" is max `|F_j(x*)|` over recovered solutions (QQ only; skipped for `ZZ/p`).

Supporting analysis (bottlenecks, comparison with Martyushev CVPR 2022 upstream, roadmap to close the size gap) is in `matrixhi_scaling_gap.md`.

## Summary table

| Problem | Ring | deg I | Default | Larsson | MatrixHi | Greedy | Paper std |
|---|---|---:|---|---|---|---|---|
| paper §3 ex (2-var) | QQ | 12 | 19×19 ✓ | 19×19 ✓ | **7×19** ✓ | 7×19 ✓ | — |
| 3-var docs | QQ | 12 | 28×28 ✓ | 28×28 ✓ | **20×33** ✓ | 20×33 ✓ | — |
| 5pt essential +det | QQ | 10 | 20×20 ✗* | 20×20 ✗* | **10×20** ✓ | 10×20 ✓ | 10×20 |
| 6 Demazure cubics | QQ | 10 | 37×34 ✗* | 34×34 ✗* | 25×35 ✓ | **24×34** ✓ | — |
| #1 F+λ 8pt | ZZ/p | 8 | 19×20 | 19×20 | **11×20** | 11×20 | 11×19 |
| #2 E+f 6pt | ZZ/p | 9 | 72×45 | 42×42 | 52×57 | **49×55** | 11×20 |
| #3 f+E+f 6pt | ZZ/p | 15 | 215×120 | 68×73 | 93×99 | **92×98** | 12×27 |
| 3-var mixed | QQ | 6 | 11×11 ✓ | 11×11 ✓ | **5×11** ✓ | 5×11 ✓ | — |

✓ = tier-3 residual < 1e-6
✗* = pre-existing `recoverSolutions` failure (see `matrixhi_scaling_gap.md §4.4`); MatrixHi/Greedy work correctly on the same problem.

## Timing (seed 42, cpuTime)

| Problem | Default | Larsson | MatrixHi | Greedy |
|---|---:|---:|---:|---:|
| paper §3 ex | 0.007 s | 0.007 s | 0.019 s | 0.071 s |
| 3-var docs | 0.051 s | 0.009 s | 0.172 s | 0.368 s |
| 5pt essential +det | 0.024 s | 0.07 s | 0.102 s | 0.396 s |
| 6 Demazure cubics | 0.029 s | 0.076 s | 0.461 s | **2.012 s** |
| #1 F+λ 8pt | 0.007 s | 0.007 s | 0.161 s | 0.205 s |
| #2 E+f 6pt | 0.021 s | 0.064 s | 0.937 s | **4.060 s** |
| #3 f+E+f 6pt | 0.036 s | 0.075 s | 2.786 s | **163.955 s** |
| 3-var mixed | 0.009 s | 0.008 s | 0.063 s | 0.021 s |

## Where Greedy strictly shrinks MatrixHi (free α > 0 regime)

| Problem | MatrixHi | Greedy | Savings |
|---|---|---|---|
| 6 Demazure cubics | 25×35 | 24×34 | 1 row, 1 col |
| #2 E+f 6pt | 52×57 | 49×55 | 3 rows, 2 cols |
| #3 f+E+f 6pt | 93×99 | 92×98 | 1 row, 1 col |

On every other problem above, free α = 0 and Greedy = MatrixHi.

## Reproducing

```bash
cd workshop

# Single-problem bench (5pt essential, all three active strategies)
M2 tests/bench_5pt_all.m2

# adjustParams shrinking demo (6 Demazure cubics, no det)
M2 tests/bench_greedy_effect.m2

# 10-problem matrix (the source of this document)
M2 tests/bench_problems.m2

# Package regression tests (19 TEST blocks)
M2 -e 'needsPackage "EliminationTemplates"; check "EliminationTemplates"; exit 0'
```
