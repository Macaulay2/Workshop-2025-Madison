# MatrixHi / Greedy — Scaling Bottleneck and Gap to CVPR 2022

This document captures *why* the MatrixHi + `adjustParams` pipeline does not scale to Martyushev CVPR 2022 template sizes on larger problems (#2 E+f 6pt, #3 f+E+f 6pt, #12 E+fλ 7pt, P4P+fr), and what specifically would need to change to close that gap.

Pair with:
- `strategy_comparison.md` — the benchmark numbers this doc reasons about.
- `paper_implementation_report.md` — the line-by-line audit of Martyushev's Maple code.

---

## 1. Pipeline stages and per-stage cost

With `n_G` = #gap polynomials, `n_F` = #generators, `|B|` = quotient basis size, `|α|` = free α parameters, `|E|` = excessive monomials:

| Stage | File / line | Cost |
|---|---|---|
| `buildGapPolys` | `MartyushevClean.m2:24` | one Gröbner basis of `J` + `n_G` normal-form reductions |
| `buildHSymbolic` | `MartyushevClean.m2:46` | per row: RREF of a \|target-mons\| × \|params-in-row\| matrix (degree bumped until feasible) |
| **`adjustParams`** | `MartyushevClean.m2:182` | **outer loop iterates ≤ MaxIter × \|E\|; inner step rebuilds a \|contributing-entries\| × \|α\| linear system per excessive; full `Hcurrent` substitution after every commit** |
| `buildTemplateFromH` | `MartyushevClean.m2:475` | `O(rows × cols)` coefficient lookups |

### 1.1 Where the time goes: `adjustParams` outer loop

Each outer pass of the greedy:

1. **`getShiftMons(Hcurrent)`** — walks every nonzero entry of the symbolic `H` and calls `coefficients(entry, Variables => baseVars)` to recompute shift monomials.
2. **`getExcessive(HM)`** — expands `m * F_j` for every shift monomial, collects all monomials, subtracts `RB` — a set difference over potentially thousands of elements.
3. For each remaining excessive monomial `e`:
   - Walks `n_G × n_F` entries of `Hcurrent`.
   - Calls `coefficient(e, m * F_j)` for every shift monomial `m`.
   - Extracts the linear dependence on the α's.
   - Solves a small linear system over the coefficient field (RREF).
4. After any α commit: re-runs `Hcurrent_(i,j) = substMap(Hcurrent_(i,j))` through the extended α-ring for every entry.

Roughly `O(n_G · n_F · (sum monomials per entry) · |E|)` ring operations per pass, plus the M2 `coefficient` / `coefficients` calls which are not cheap.

### 1.2 Empirical numbers

- `#3 f+E+f 6pt` (ZZ/p): one row/col of savings cost **164 s** for the Greedy pass.
- `P4P+fr` (target 52×68): we expose 1095 free α's but the greedy commits only a fraction (see `paper_implementation_report.md §6`).
- `#12 E+fλ 7pt`: did not complete MatrixHi (no adjustParams) over `ZZ/32749` in 4+ minutes; killed.

---

## 2. Does upstream Maple have the same cost?

**Yes — the core algorithm is identical, so the asymptotic cost is identical.** The reference in `reference_martyushev/_greedyAG/` contains the same `matrixHi` + `adjustParams` procedure we ported.

But Martyushev doesn't run matrixHi + adjustParams **from scratch on a standard basis**. The reference directory includes three upstream accelerators:

- `basisFinder.mw` — searches non-standard bases (different monomial orders, different "selected" monomials for `B`).
- `monOrdFinder.mw` — searches monomial orders.
- `bases/b_<problem>` — precomputed lists of candidate bases per problem (one file per benchmark, many bases each).

For each problem the Maple pipeline runs matrixHi + adjustParams on **many candidate bases**, keeps the smallest template. The paper's "nstd" column in Table 1 is the result of that search; it can be substantially smaller than the "std" column (e.g. #1 F+λ 8pt: std 11×19 → nstd 7×15).

Our code uses whatever `basis(R/J)` gives — effectively GRevLex. On `#2 E+f 6pt` that's a ~52-monomial basis, vs the paper's ~11-monomial basis. A 5× inflation in `|B|` cascades multiplicatively into `|E|` and `|α|`, and `adjustParams` starts from an exploded state.

---

## 3. Three additional gaps, beyond basis choice

Even on the same basis, our implementation is slower and eliminates less than the reference. From `workshop/docs/paper_implementation_report.md §6` and our own profiling:

### 3.1 Weaker greedy heuristic
Our `adjustParams` sorts excessive monomials by "contributing entry count" and tries them one at a time. Martyushev's heuristic additionally considers:
- Degree patterns of the contributing entries.
- Coverage: how many already-fixed basis monomials the candidate α interacts with.
- Lookahead: whether committing α_k *creates* new excessives elsewhere.

On P4P+fr we eliminate roughly a few dozen of 1095 free α's; the reference eliminates nearly all.

### 3.2 No `ZZ/p` specialization inside the greedy
Our code runs in whatever ring the user picks. The benchmarks use `ZZ/32749` because `QQ` explodes during RREF (rational coefficient blowup). Martyushev's code is pinned to `ZZ/p` for all the heavy arithmetic and only lifts back to the intended ground field at the very end.

We could do the same — run `adjustParams` over `ZZ/p` regardless of input ring, then rebuild the final `H` over `QQ` — but we haven't wired it.

### 3.3 Repeated M2 `coefficient` calls inside tight loops
A profiling pass would likely show `adjustParams` spending most of its wallclock time in `coefficient` / `coefficients` / `monomials` — generic M2 ring operations — rather than in useful algebra. Each outer iteration re-extracts everything from scratch.

Caching per-monomial contribution tables across iterations (instead of rebuilding them) would be a substantial constant-factor win.

---

## 4. Roadmap for closing the gap

In priority order, highest expected-impact first:

### 4.1 Port basis / monomial-order search (biggest win)

**Context — what the paper actually does.** Martyushev CVPR 2022 §5 (note
1) reports every template in two columns per problem:
- **std** — smallest template across the *entire Gröbner fan* (via Gfan)
  **or** across 1,000 randomly-selected Gröbner bases when Gfan cannot
  terminate in reasonable time.
- **nstd** — smallest template across 500 *non-standard* quotient bases
  sampled via the random strategy from Larsson CVPR 2018.

To reach those numbers, matrixHi + adjustParams is run against every
candidate basis from the chosen route and the smallest template is
kept. Our current code runs matrixHi on exactly one basis (GRevLex, from
`basis(R/J)`), i.e. roughly 1/500 to 1/1,000 of the search the paper
does. That is the dominant source of the gap on `#2 E+f 6pt`,
`#3 f+E+f 6pt`, `#6 P4P+fr`, `#12 E+fλ 7pt`, etc.

**The reference `_greedyAG/` directory caches artifacts for four
different routes.** The paper uses whichever route is tractable per
problem:

| Route | Cache | # cached problems | What it unlocks |
|---|---|---:|---|
| (a) Gfan + basisFinder (full treatment) | `gfan/gf_<prob>` + `bases/b_<prob>` | 8 | std column when Gfan terminates |
| (b) 1,000 random Gröbner bases | (no cache — generated on demand) | — | std column when Gfan does not terminate |
| (c) 500 random non-standard bases (Larsson 2018) | (no cache — generated on demand) | — | nstd column whenever nstd < std |
| (d) Weighted-degree orderings | `weights/w_<prob>` | 12 | alternate route for specific problems |
| (e) Precomputed change matrices | `cm/cm_<prob>` | 2 | shortcut for 5p and toy |

**A candidate basis is "good"** only empirically: it must be a set of
`|B| = degree(I)` monomials that spans `K[X]/J`. For std bases that is
automatic (every Gröbner basis produces one); for nstd bases you pick
random monomial subsets and keep those with full-rank coefficient
matrix modulo `J`. No closed-form heuristic ranks candidates — you run
matrixHi + adjustParams on each and sort by template size. That is why
`bases/b_5p` has 200 candidates rather than 1.

**The 8 problems with a cached `bases/b_<prob>` file** are mostly the
smaller benchmarks:

| Cache file | Paper problem | Paper target (std / nstd) |
|---|---|---|
| `b_5p` | 5-point essential (§3 reference) | 10×20 / 10×20 |
| `b_8ptF_radial` | #22 λ+F+λ 8pt (both-sided radial) | 31×47 / 31×47 |
| `b_stitching` | #5 Stitching f+R+f+λ 3pt | 48×66 / 18×36 |
| `b_4pra` | #25 Rel. pose E+angle 4pt v2 | 16×36 / 16×36 |
| `b_3pra_st0` | 3-pt relative pose variant (supplementary) | — |
| `b_toy` | toy 2-var example | — |
| `b_wpnp` | "weighted PnP" variant | — |
| `b_wpnp_2x2sym` | symmetric-weighted PnP variant | — |

None of the hard problems (`#2`, `#3`, `#6`, `#12`) are in `bases/`.
Reaching paper sizes on those requires route (b) or (c), which do not
ship with a cache.

**Realistic port, in priority order:**

1. **Stage A — cache consumer (route a).** Read `bases/b_<prob>`, loop
   over its candidates, run matrixHi + adjustParams, return smallest
   template. ~100 LOC. Covers exactly those 8 problems. Validates the
   mechanism empirically before any bigger investment.
2. **Route (b) — random Gröbner-basis sampling.** Draw 1,000 random
   monomial orders, compute the reduced Gröbner basis for each, treat
   the resulting standard basis as a candidate, run matrixHi +
   adjustParams. ~200 LOC. This is the route that actually closes the
   gap on `#6 P4P+fr`, `#12 E+fλ 7pt`, and every large problem where
   Gfan itself would not finish.
3. **Route (c) — Larsson 2018 non-standard sampling.** Sample random
   `|B|`-monomial subsets, test full-rank against `R/J`, keep 500 valid
   ones, run matrixHi + adjustParams on each. ~200–300 LOC. Needed to
   match paper's strictly-smaller nstd column on the ~15 problems
   where that column is smaller than std.
4. **Route (d) — weighted-degree orderings.** ~100 LOC plus a weight-
   vector parser. Low marginal value; most of the 12 weighted-cached
   problems are also reachable via (b) or (c).
5. **Route (e).** Skip — already subsumed by our graph-ideal `getH0`.

Sweet-spot path: **A → (b) → (c)**. Total ~600 LOC. Would match paper
sizes on most of the `strategy_comparison.md` bench. Stage A alone
unlocks 8 problems but only a subset of those have template sizes we
don't already match (we already reproduce `5p` at 10×20 via MatrixHi /
Greedy without any basis search).

### 4.2 Tighten `adjustParams` heuristic
Port the sort / lookahead criteria from Martyushev's `adjustParams`. Cache per-monomial contribution tables across outer iterations.

Expected win: measurable size improvement on problems with large `|α|` (P4P+fr etc.); speed improvement on everything.

Scope: medium — a rewrite of the inner loop in `MartyushevClean.m2` (~100 lines).

### 4.3 Pin α-extension ring to `ZZ/p`
Make `adjustParams` always run over `ZZ/p`, regardless of input ring. Rebuild the final `H` over the user's ring at the end.

Expected win: **order-of-magnitude speedup on `QQ` problems** by avoiding rational coefficient blowup during RREF.

Scope: small — ~50 lines, mostly plumbing.

### 4.4 Fix graph-ideal `recoverSolutions` for non-generic actions
`EliminationTemplates.m2:372` assumes the action variable is a random linear form. Single-variable actions on 3+-variable systems silently return wrong solutions or error with "cannot coerce CC value to ring type." The pre-existing bug is documented in `tasks.md` as an open failure.

**Status:** not fixed. MatrixHi / Greedy bypass the bug via `recoverSolutionsMatrixHi`. The Default / Larsson path still errors or silently fails. A proper fix rewrites `recoverSolutions` to handle variables that aren't in the quotient basis.

Scope: medium — requires understanding the graph-ideal pipeline's monomial-partition structure. Not on the critical path for MatrixHi/Greedy but matters for users of Default/Larsson on those problem shapes.

### 4.5 Support polynomial action variables in MatrixHi/Greedy
`buildGapPolys` requires each `a·b_i` to be a single monomial. Workaround: do the graph-ideal extension inside the MatrixHi/Greedy pipeline (add `s − a = 0` to the ideal; use `s` as monomial action in the augmented ring; strip `s` when reading variable coordinates back).

Expected win: users can pass arbitrary action linear forms (e.g. `random(1, R)`) to MatrixHi/Greedy, matching the Default / Larsson API surface.

Scope: small — ~50 lines in `MartyushevClean.m2`.

---

## 5. What would success look like?

After Stage A (4.1 route a) alone: we can reproduce paper sizes on the
8 problems with cached `bases/` files. Most impact is on problems where
we currently use GRevLex and Martyushev's best basis differs from it.
Several of those 8 (like `5p`) we already match without any basis
search, so Stage A measures the *upper bound* on what basis-search
buys us on those problems rather than being a guaranteed win.

After Stage A + route (b): we match paper **std** sizes on essentially
all Table 1 / Table 2 problems — including `#6 P4P+fr` and `#12 E+fλ 7pt`
where the current output is ~6× too large.

After A + (b) + (c): we also match paper **nstd** on the ~15 problems
where nstd is strictly smaller than std (e.g. `#1 F+λ 8pt` 11×19 →
7×15; `#5 Stitching` 48×66 → 18×36).

After (4.2) + (4.3) on top: `adjustParams` runtime drops an order of
magnitude through `ZZ/p` pinning and per-iteration caching. Still
slower than Maple on very large problems where `|α| > 1000` — that is
the inherent algorithm, not an implementation gap.

None of these is a session's work. Stage A is roughly a day; each of
routes (b), (c), and items 4.2/4.3 is its own project with its own
validation plan.
