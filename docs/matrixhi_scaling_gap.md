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
Replicate `basisFinder` / `monOrdFinder` logic from Maple. For each benchmark, enumerate candidate bases (possibly from the cached `bases/b_<problem>` files), run matrixHi + adjustParams on each, keep the smallest template.

Expected win: **10× template-size reduction on the problems where we're currently far from paper** (#2, #3, #12, P4P+fr). This is the single intervention that would let us match paper sizes on most Table-1 entries.

Scope: substantial — several hundred lines, likely a new file `BasisFinder.m2`. Requires understanding of Gröbner fan / monomial-order geometry.

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

After (4.1) + (4.2) + (4.3):

- Match paper std sizes on #1–#5, #7, #22 (most of Table 1 row count).
- Close to paper nstd on problems where we've also ported basis search.
- `adjustParams` runtime drops an order of magnitude through `ZZ/p` pinning and caching.
- Still slower than Maple on very large problems (#6, #8, #18–#21) where `|α| > 1000`; that's the inherent algorithm, not an implementation gap.

Not a session's work. Each of 4.1, 4.2, 4.3 is its own project with its own validation plan.
