# Strategy Comparison — Benchmark Results

Measured by `benchmarks/bench_problems.m2` with `setRandomSeed 42`. Template
sizes are `rows × cols`. `build+solve` is `cpuTime` in seconds. `res` is
`max |F_j(x*)|` over recovered solutions (QQ only; skipped for ZZ/p).

## Strategy labels

Every strategy flows through the same `getTemplate → [E | R | B] fold →
extractActionFromTemplate` pipeline. Only the `H0` computation inside
`getH0` differs per strategy:

| Label | Option combination | `H0` body |
|---|---|---|
| Default | _no options_ (`Strategy => null`) | Gröbner change-of-basis |
| Larsson | `Strategy => "Larsson"` | Default's H0, then reduced modulo the first syzygy module of F |
| MatrixHi-mode | `Strategy => "Greedy", AdjustParams => false` | Martyushev gap polynomials + `buildHSymbolic` with α := 0 (particular solution) |
| Greedy | `Strategy => "Greedy"` (default `AdjustParams => true`) | Martyushev gap polynomials + `buildHSymbolic` + `adjustParams`; QQ input is run over ZZ/p internally and the H0 support lifted back |

Polynomial actions are lifted to $R_s = R[s] / \langle s - a \rangle$
automatically, so every strategy accepts both monomial and polynomial
action variables.

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
- ✗ = pre-existing recovery failure (every strategy returns the same partial
  solution).

Key observations:

- On monomial-action problems with 0 free α's (paper §3 ex, 5pt +det, F+λ 8pt,
  3-var mixed), all four strategies reach the same size.
- On free-α problems where `adjustParams` commits parameters (6 Demazure,
  #2 E+f, #3 f+E+f), Greedy strictly wins. The margin between Greedy and
  MatrixHi-mode on #2 / #3 comes from `adjustParams` cancelling redundant
  shift support.
- On 3-var docs, Default's Gröbner H0 has smaller monomial support than
  `buildHSymbolic`'s — one of the few cases where Default beats Greedy.

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
commits. ZZ/p pinning inside `adjustParams` keeps the inner RREF from
blowing up on QQ input but doesn't change the per-iteration work.

## Polynomial action

Polynomial action variables like `x + 4y` or `random(1, R)` flow through
the same `[E | R | B]` + RREF + pivot path as monomial actions. The
action is lifted to `s` in `R_s = R[s] / ⟨s − a⟩`, where `s` is itself a
ring variable on the quotient.

| Problem | Action | Template | Residual |
|---|---|---|---|
| `change of ideals` (QQ[x,y], `x²+y²−1`, `x²+y³+xy−2`) | `x + 4y` | 20×26 | ~1e-14 |
| 5pt essential (Demazure cubics, no det) | `random(1, R)` | 37×34 | ~1e-12 (one outlier ~1e-8) |

Both paths run RREF over QQ natively. An LU split-solve fallback in
`getActionMatrix(EliminationTemplate)` activates (via `try/else`) if RREF
raises on adversarial QQ input.

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

# Single-problem bench (5pt essential)
M2 benchmarks/bench_5pt_all.m2

# adjustParams shrinking demo (6 Demazure cubics, no det)
M2 benchmarks/bench_greedy_effect.m2

# Full 10-problem matrix (source of this document)
M2 benchmarks/bench_problems.m2

# Package regression tests
M2 -e 'needsPackage "EliminationTemplates"; check "EliminationTemplates"; exit 0'
```
