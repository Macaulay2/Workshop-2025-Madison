# MatrixHi / Greedy — Scaling Bottleneck and Gap to CVPR 2022

This document reasons about where our Greedy pipeline matches Martyushev
CVPR 2022 template sizes, where it still falls short, and what would need
to change to close the remaining gap.

Pair with:
- `strategy_comparison.md` — the benchmark numbers this doc reasons about.
- `paper_implementation_report.md` — the line-by-line audit of Martyushev's
  Maple code.

## 0. Snapshot (post-Increment 3)

Every strategy (Default / Larsson / MatrixHi-mode / Greedy) now flows
through Direction A's `[E | R | B]` layout + `extractActionFromTemplate`
(RREF + pivot) downstream; only `getH0` branches per strategy. Polynomial
actions are lifted to `R_s = R[s] / ⟨s − a⟩` so every strategy accepts
both monomial and polynomial actions (Increment 3). ZZ/p pinning inside
Greedy's `adjustParams` keeps QQ templates numerically safe at `getH0`
time (Increment 2).

Resolved since the original draft of this doc:

- §4.2 **Tighten `adjustParams` heuristic** — landed. Commit `cac2990` ports Martyushev 3.3b (DESC sort by `|rm|, np`, early break); commit `19ab0f2` adds 3.3c caching of per-monomial contribution tables across outer iterations.
- §4.3 **Pin α-extension ring to ZZ/p** — landed in commit `e5e4e4d` (Increment 2). Greedy's `adjustParams` runs over ZZ/p on QQ input, lifts the resulting H0 support back to QQ per-coefficient.
- §4.4 **Fix graph-ideal `recoverSolutions` for non-generic actions** — landed. `recoverSolutions` rewritten in commit `62087d9` (Increment 3) to use RREF + pivot on the `[E | R | B]` template; handles monomial- and polynomial-action cases, over-determined shift systems, and degenerate eigenvectors (skipped when the `1` slot is near zero).
- §4.5 **Support polynomial action variables in Greedy** — landed in commit `62087d9` (Increment 3). Action `a ∈ R` lifts to the new variable `s ∈ R_s`; after the lift `s` is a ring variable, `isMonomialAction(s)` is true, and the full Greedy pipeline works on `(s, J_s)` in `R_s`.

What is **not** resolved: template sizes on the hard ZZ/p problems
(`#2 E+f 6pt` at `21×30` vs paper `11×20`, `#3 f+E+f 6pt` at `31×50` vs
paper `12×27`). The gap is not in `adjustParams` itself — it's in the
**basis choice**. See §2 below.

## 1. Pipeline stages and per-stage cost

With `n_G` = #gap polynomials, `n_F` = #generators, `|B|` = quotient basis
size, `|α|` = free α parameters, `|E|` = excessive monomials:

| Stage | File / line | Cost |
|---|---|---|
| `buildGapPolys` | `MartyushevClean.m2:24` | one Gröbner basis of `J` + `n_G` normal-form reductions |
| `buildHSymbolic` | `MartyushevClean.m2:46` | per row: RREF of a `\|target-mons\| × \|params-in-row\|` matrix (degree bumped until feasible) |
| **`adjustParams`** | `MartyushevClean.m2:182` | **outer loop iterates ≤ `MaxIter × \|E\|`; inner step rebuilds a `\|contributing-entries\| × \|α\|` linear system per excessive; ZZ/p-pinned when input is QQ** |
| `[E \| R \| B]` fold | `EliminationTemplates.m2` (`getTemplateMatrix`) | `O(rows × cols)` coefficient lookups |
| `extractActionFromTemplate` | `MartyushevClean.m2:695` | single RREF on the template; pivot lookup per basis monomial |

### 1.1 Where Greedy's time goes

Each outer pass of `adjustParams`:

1. `getShiftMons(Hcurrent)` — walks every nonzero entry of the symbolic `H`
   and calls `coefficients(entry, Variables => baseVars)` to recompute
   shift monomials.
2. `getExcessive(HM)` — expands `m * F_j` for every shift monomial,
   collects all monomials, subtracts `RB` — a set difference over
   potentially thousands of elements.
3. For each remaining excessive monomial `e`:
   - Walks `n_G × n_F` entries of `Hcurrent` (now cached across iterations
     — Increment 3.3c).
   - Calls `coefficient(e, m * F_j)` for every shift monomial `m`.
   - Solves a small linear system over ZZ/p (pinned — Increment 2).
4. After any α commit: re-runs `Hcurrent_(i,j) = substMap(Hcurrent_(i,j))`.

Empirical numbers on seed 42:
- `#3 f+E+f 6pt`: Greedy pass 149 s (down from 164 s pre-caching; the
  extra savings vs size improvement is marginal — see §2).
- `#2 E+f 6pt`: 45.9 s.
- `6 Demazure cubics`: 40 s.

## 2. The remaining gap is basis choice, not `adjustParams`

Martyushev doesn't run `matrixHi + adjustParams` **from scratch on one
standard basis**. The reference pipeline runs it against many candidate
bases per problem and keeps the smallest template:

- `basisFinder.mw` — searches non-standard bases (different monomial
  orders, different "selected" monomials for `B`).
- `monOrdFinder.mw` — searches monomial orders.
- `bases/b_<problem>` — precomputed lists of candidate bases per problem.

The paper's "nstd" column in Table 1 is the result of that search; it can
be substantially smaller than the "std" column (e.g. `#1 F+λ 8pt`: std
11×19 → nstd 7×15).

Our code uses whatever `basis(R/J)` gives — effectively GRevLex. On
`#2 E+f 6pt` that's a ~52-monomial basis, vs the paper's ~11-monomial
basis. A ~5× inflation in `|B|` cascades into `|E|` and `|α|`; `adjustParams`
starts from an already-exploded state and can only shrink it marginally.

This is the dominant source of the residual gap on `#2`, `#3`, `#6 P4P+fr`,
`#12 E+fλ 7pt`, etc. `adjustParams` improvements alone cannot close it.

## 3. Roadmap for closing the remaining gap

In priority order, highest expected-impact first:

### 3.1 Port basis / monomial-order search (biggest remaining win)

**Context — what the paper does.** Martyushev CVPR 2022 §5 (note 1) reports
every template in two columns per problem:
- **std** — smallest template across the entire Gröbner fan (via Gfan) or
  across 1,000 random Gröbner bases when Gfan cannot terminate in
  reasonable time.
- **nstd** — smallest template across 500 non-standard quotient bases
  sampled via the strategy from Larsson CVPR 2018.

To reach those numbers, `matrixHi + adjustParams` is run against every
candidate basis and the smallest template is kept. Our code runs on
exactly one basis (GRevLex), i.e. roughly 1/500 to 1/1,000 of the paper's
search.

Infrastructure we have (post-Increment 4.1 commits `b8c69db`):
- `buildGreedyTemplateWithBasis` — Greedy on a supplied basis.
- `searchBases` — run Greedy over a list of bases, return the smallest.
- `randomStandardBases` — sample `n` random standard bases via random
  weight vectors; per-candidate `cpuTime` budget via `alarm` guards
  against pathological weights (commit `44173f8`).
- `randomNonstandardBases` — sample `n` random non-standard bases per
  Larsson 2018.
- `parseMartyushevBases` — read the cached `bases/b_<prob>` files.

What's missing: a driver that wires these together into the paper's
std / nstd columns for each benchmark. Remaining work is ~200 LOC plus
a validation harness; not in Increment 3's scope.

### 3.2 Polynomial-action support inside `buildGreedyTemplateWithBasis`

The basis-search primitive requires a monomial action (its `buildGapPolys`
call assumes `a * b_i` is a single monomial). Increment 3 lifted actions
inside the main `getTemplate`, but not inside this primitive. Options:
1. Lift the action inside `buildGreedyTemplateWithBasis` mirroring
   Increment 3's code in `getTemplate`.
2. Refactor `buildGreedyTemplateWithBasis` to call through
   `getTemplate` / `copyTemplate`, so the lift happens once.

Low priority — no current user runs polynomial-action basis search.

### 3.3 CRT-lift from ZZ/p for RREF numerical safety

`extractActionFromTemplate` and the rewritten `recoverSolutions` both run
RREF on QQ templates. For adversarial polynomial-action problems on QQ
the rational-coefficient blowup inside `rawLinAlgSolve` can be slow or
hit the SIGSEGV boundary.

We already have two guards:
- `getActionMatrix(EliminationTemplate)` wraps `extractActionFromTemplate`
  in a `try / else` that falls back to origin's LU split-solve on an
  `[E | B]` rebuild if RREF raises a catchable error (commit `44173f8`).
- `recoverSolutions` reduces each excess / residual column via its RREF
  pivot row; no least-squares solve.

The "real" fix for numerical safety across every RREF in the pipeline
would be multi-prime CRT: run RREF in ZZ/p for k distinct primes,
combine via the Chinese Remainder Theorem, rational-reconstruct the QQ
entries. Distinct from Increment 2's ZZ/p pinning (which only reads
support, not QQ values). Cost is substantial (~200 lines, including
prime selection and unlucky-prime detection); defer until an adversarial
problem motivates it.

## 4. What success would look like

- **Today**: matches paper `std` on every problem where the paper's best
  basis equals GRevLex (5pt+det, paper §3, F+λ 8pt std column, 3-var
  mixed).
- **After 3.1 (basis search driver)**: matches paper `std` on the
  remaining table problems. Matches `nstd` on the ~15 problems where
  `nstd < std`.
- **After 3.2 + 3.3**: not about size — about robustness. Polynomial-action
  basis search becomes available; QQ RREF numerical safety is tightened
  beyond the current LU-fallback level.

None of these is a session's work. Basis-search driver is about a day
plus validation; CRT-lift is its own project.
