# Report — How the paper's reference implementation actually works

**Source audited:** `reference_martyushev/_greedyAG/greedyAG.mw` (Martyushev et al., CVPR 2022).
**Purpose:** document exactly what the original Maple code does, for comparison with our M2 ports.

---

## TL;DR

The paper's production pipeline is:

```
matrixHi  →  adjustParams  →  constructTemplate
```

- No $H_0$-plus-$\Theta\cdot H_1$ syzygy decomposition is used in the production path.
- The "parameters" the greedy adjusts are the **undetermined scalar coefficients** introduced inside `matrixHi`, *not* a separate syzygy module.
- A syzygy-based variant (`matrixHsyz`) exists in the same file but is **never called by the driver**. It appears to be a dead-code artifact from an earlier design.

---

## 1. The three H-constructors in the Maple source

A grep of procedure definitions in `greedyAG.mw` yields:

| Procedure | Comment in source | Used by driver? |
|---|---|---|
| `matrixHi` | "iteratively construct matrix H of size nG x nF" | **yes** |
| `matrixHni` | "use known change matrix to non-iteratively construct matrix H" | **yes** |
| `matrixHsyz` | "syzygy-based construction of matrix H" | **no — orphaned** |

### The driver, verbatim

```maple
H := `if`(basis = nstd,
          matrixHi(F, G, vars, pp),
          matrixHni(F, G, C, CM, vars, pp)):
# construct matrix H
shifts, M, V := adjustParams(F, H, RB, vars, pp):
# adjust parameters and return template M and monomial vector V
```

The `if` chooses between `matrixHi` (non-standard basis) and `matrixHni` (standard basis with a pre-computed change matrix). Neither branch touches `matrixHsyz`.

### What `matrixHsyz` would do, if it were called

```maple
matrixHsyz := proc(F, G, C, A, tord, pp) local S, Theta:
    S     := syzModule(F, G, A, tord, pp):                  # syzygy module of F
    Theta := Matrix(..., (i,j) -> t[i,j]):                  # scalar unknowns t[i,j]
    return map(expand, C.A + Theta.S) mod pp:               # H = C·A + Theta · S
end proc:
```

This is **exactly the $H = H_0 + \Theta \cdot H_1$ form** — with `C·A` playing the role of $H_0$ and the syzygy module `S` playing the role of $H_1$.

**Crucial observation:** the paper's author wrote this variant, did not choose it for production, and dropped it silently. Our analysis of why (below) matches what the author must have discovered: scalar $\Theta$ cannot decouple a polynomial $H_1 = S$.

---

## 2. What `matrixHi` does

Full body (de-obfuscated from the Maple):

```maple
matrixHi := proc(F, G, vars, pp)
  undetPoly := proc(vars, pv, d, j) local mons:
    # generic polynomial of degree d with undetermined scalars _v[j,k]
    coeffs(expand(1 + add(pv^p, p=1..d)), vars, 'mons'):
    mons := [mons]:
    add(mons[k] * _v[j,k], k=1..nops(mons))
  end proc:

  for i from 1 to nG do                              # one loop iter per gap polynomial G[i]
    d_min := 0
    dG    := [max(0, degree(G[i]) - degree(F[j])) for j in 1..nF]

    for d from 0 while d - d_min <= 0 do             # degree-bumping loop
      uH   := [undetPoly(vars, pv, dG[j]+d, j) for j in 1..nF]
      leqs := coefficients-of( -G[i] + sum_j uH[j]*F[j],  vars )
      sol  := msolve(leqs, pp)                       # modular linear solve

      if nops({sol}) = 0 then
        d_min := d + 1                               # infeasible → force retry at d+1
      end if
    end do

    h[i] := eval(uH, sol) mod pp                     # substitute sol (free α stay symbolic)
  end do

  return [row-vector h[i] for i in 1..nG]
end proc
```

### Key properties

1. **Row-by-row.** One independent subproblem per gap polynomial $G_i$.
2. **Undetermined scalars per monomial.** Each entry $H_{i,j}$ is a generic polynomial of degree $d_G[j] + d$ with a scalar unknown $v_{j,k}$ for every monomial $\mathrm{mon}_k$.
3. **Linear system in the unknowns.** `leqs` is the equational system obtained by equating coefficients of $-G_i + \sum_j H_{i,j} F_j$ to zero, monomial by monomial.
4. **Modular linear solve.** `msolve` over $\mathbb{F}_{pp}$ returns a *parametric* solution that retains free variables as symbolic `_t1, _t2, ...`.
5. **Degree-bumping termination.** The loop exits at the *first* $d$ for which `msolve` returns a solution. It does not try to overshoot.
6. **Free α's survive.** Because `eval(uH, sol)` substitutes the parametric `sol` (which keeps unbound parameters symbolic), the returned `h[i]` still contains undetermined indeterminates wherever the linear system had null space. Those are the handles `adjustParams` will pull on.

### How the loop actually terminates

| Step | d, d_min | `msolve` | d_min update | Next check | Action |
|---|---|---|---|---|---|
| 1 | d=0, d_min=0 | solvable | unchanged | d=1: 1−0=1>0 | **exit** |
| 1 | d=0, d_min=0 | empty | d_min=1 | d=1: 1−1=0≤0 | retry at d=1 |
| 2 | d=1, d_min=1 | solvable | unchanged | d=2: 2−1=1>0 | **exit** |
| 2 | d=1, d_min=1 | empty | d_min=2 | d=2: 2−2=0≤0 | retry at d=2 |

So `matrixHi` stops at the minimum feasible degree per row.

---

## 3. What `adjustParams` does

Full body (de-obfuscated):

```maple
adjustParams := proc(F, H, RB, vars, pp)
  # H1 starts as a WORKING COPY of H. The name "1" is misleading —
  # it is NOT a syzygy matrix. Just a local variable.
  H1 := H:

  # --- Phase 1: column-decompose H, retain only columns whose coefficients
  #              still contain indeterminate α's ---
  Lm   := seq(poly2matrix(H[..,j], vars), j=1..nF)
  cols := [keep only columns whose entries are not already numeric for each j]
  L, mon := [coefficient matrices on those cols], [monomial vectors]

  # --- Phase 2: compute excessive monomials E = V₀ \ RB ---
  sm := [supports of expand(F[j] * mon[j][k])  for all j, k]
  V0 := union of all of sm
  E  := V0  \  RB

  # --- Phase 3: sort E by "difficulty": fewer source positions first ---
  rm := for each e in E, list of (j, k) positions producing e
  np := downstream impact counts
  E  := sort E by (|rm|, np) ascending

  # --- Phase 4: greedy loop ---
  for e in E do
    Lm  := refresh column-decomposition of H1       # H1 shrinks each iter
    L, mon := ...
    if nops(indets(L)) = 0 then break end if         # no free α left

    rm   := positions (j,k) where mon[j][k]·F[j] contains e
    leqs := { α-expressions at those positions }     # set them all to 0 ⇒ e vanishes
    sol  := msolve(leqs, pp)

    if nops({sol}) = 0 then next end if              # skip: e cannot be killed now
    H1 := map(expand @ eval, H1, sol) mod pp         # commit: substitute α values
  end do

  constructTemplate(H1, F, vars, RB)
end proc
```

### Key properties

1. **Operates on `matrixHi`'s output directly.** The working matrix `H1` is a copy of the H produced upstream; **no separate $H_1$ syzygy module exists in this algorithm**. The "1" in the name is a source of confusion I'll just call out once.
2. **The parameters adjusted are the undetermined α's retained inside H from `matrixHi`.** Setting one α to a value is cheap: it zeros exactly one coefficient of one monomial in one entry. No cascading, because α's are scalars multiplying monomials, not scalars multiplying polynomials.
3. **Excessive monomials drive the loop.** For each $e \in E$, the greedy asks: *"which α-expressions control this monomial's presence? Can I set them all to zero simultaneously?"* — that's the linear system `leqs`. If `msolve` finds a solution, it's committed; if not, `e` is skipped.
4. **Sort order is fixed before the loop.** $E$ is sorted once, by fewest source positions ascending. The sort does not refresh as H1 evolves. This is a potential weakness but also ensures determinism.
5. **No rollback, no backtracking.** Every successful `msolve` is committed unconditionally. There is no "does this make the template bigger?" check.
6. **Termination:** three exits.
   - `E` exhausted (natural end).
   - `indets(L) = 0` → all α pinned, early break.
   - `msolve` empty for one `e` → `next`, skip without side effect.

### Important subtlety about what gets eliminated

The algorithm's **goal is to eliminate excessive monomials**, not α's. α's are spent whenever it buys an elimination. If E exhausts while some α are still free, those α stay symbolic; `constructTemplate` then uses the monomial support of H1 (already frozen) and the remaining symbolic α's appear as coefficients in the final template — they are replaced numerically by test-data evaluation during `isActionMatrixFound`, but for template-size counting they cost nothing more.

---

## 4. Why the paper doesn't use the syzygy form

The contrast comes out clearly on a concrete problem like P4P+fr:

| Parametrization | Per-parameter effect | Free params available |
|---|---|---|
| **matrixHi** (per-monomial α) | Zeros one coefficient of one entry; no cascading | **1095** |
| **matrixHsyz** ($H_0 + \Theta\cdot H_1$) | Changes every monomial in one syzygy row; cascading | **0** on many problems, always coupled |

The `matrixHsyz` form looks algebraically elegant — parametrize the solution space using the syzygy module — but because the syzygy module has **polynomial** entries, a single scalar $\theta_{i,s}$ shifts coefficients of every monomial in that syzygy row in lockstep. To zero one excessive monomial you'd typically need to set $\theta$ to a value that un-zeros a different one, or the system is over-constrained and forces $\theta = 0$.

`matrixHi`'s per-monomial scalars avoid this entirely: the α's are linearly independent by construction (they're literally the null-space basis of a linear system), so zeroing one costs nothing elsewhere.

This is the reason the paper's production path is `matrixHi + adjustParams` and not `matrixHsyz + adjustParams`.

---

## 5. Comparison to the dropped θ-greedy prototype

An earlier `Strategy => "Greedy"` (removed 2026-04-16) did:

```m2
H1 := transpose sub(syz(gens J), ring J);
Theta := genericMatrix(ThetaExt, ..., numcols H0, numrows H1);
H := H0e + Theta * H1e;                   -- H = H0 + Θ·H1
```

This is **exactly `matrixHsyz`** — the variant the paper's author dropped. All the engineering layers (shortlisting, lookahead, transaction rollback, seed restarts) were refinements of a search over $\theta$ for a parameterization that is structurally incapable of matching paper sizes on any problem with nontrivial syzygy coupling. No improvement in the *search* part would have closed the gap. The fix is to switch parameterization — which is what the current `"MatrixHi"` / `"Greedy"` strategies do (both built on `buildHSymbolic`).

---

## 6. What the M2 port has correct, and what's missing

### Correct (matches the paper)

- **`buildHSymbolic` in `MartyushevClean.m2`** ports `matrixHi` with a flat extended ring, the degree-bumping loop, and the per-monomial α linear system. Returns `(Hsym, Rext, alphaVars, perRowData)` with free α's retained as ring variables.
- **`adjustParams` in `MartyushevClean.m2`** commits free α's to zero excessive monomials (Martyushev CVPR 2022 §4). Empirically eliminates ~10× more excessive monomials than the earlier Maple-faithful `adjustParams` port.
- **`buildTemplateFromH` in `MartyushevClean.m2`** assembles the template matrix over $\KK$ from an $H$ matrix and the shift monomials — the "constructTemplate" step of the paper.
- **Solve-side extraction** (RREF pivots per residual column → action matrix → eigendecompose) is in `benchmarks/bench_5pt_all.m2:buildMatrixHiLike`. Verified on 5pt essential over $\mathbb{Q}$.

### Missing / broken

- **`getActionMatrix` / `templateSolve` do not support the MatrixHi / Greedy cached templates.** Those paths bypass the graph-ideal extension that the Default / Larsson pipeline relies on to read the action matrix off the template. Users must extract the action matrix separately (see `bench_5pt_all.m2`).
- **`adjustParams` under-eliminates on large problems.** On P4P+fr (target 52×68), we get 339×332 — the 1095 free α's are correctly exposed but the greedy's sort heuristic or `msolve`-equivalent pivot choice misses most of them. Narrowing this gap is the main remaining algorithmic task.

### Plumbing fix (does not require new math)

To make `templateSolve(E, Strategy => "MatrixHi")` work end-to-end:

1. Add a MatrixHi / Greedy branch to `getActionMatrix(EliminationTemplate)` that calls `getTemplateMatrix` and then runs the RREF extractor in `bench_5pt_all.m2` to produce the action matrix.
2. Add a matching branch to `getEigenMatrix` / `templateSolve` that skips `getH0` / `getTemplate` for this family (those extend with the graph ideal, which MatrixHi deliberately avoids).

---

## 7. The shape/size difference in one picture

On the 5-point essential problem ($|B| = 10$, $\deg J = 10$):

```
Default pipeline:
  getH0 → (10×10 matrix, trivial cols are zero)
      → shifts (from non-trivial cols only)
      → graph-ideal extension: adds |B|=10 rows for (s-a)·b
      → template 20×20

MatrixHi pipeline:
  matrixHi → (1 row per gap polynomial, here 10 × nGens = 10×3)
      → shifts (direct from H)
      → NO graph-ideal extension
      → template 10×20
```

The $|B|$-row difference is exactly the graph-ideal extension. It exists in Default because `getActionMatrix` reads the action from the $(s-a)\cdot b$ rows via one linear solve; MatrixHi requires a different action-matrix extraction (`extractAction`), which is why `templateSolve` can't be routed through the same pipeline.

---

## 8. One-sentence summary

> **The paper's production algorithm is `matrixHi + adjustParams` — undetermined scalar α per (entry, monomial) pair, solved row-by-row via modular linear algebra, then a per-excessive-monomial greedy that commits α assignments via `msolve`; the syzygy-based variant exists in the same source file but is never called, presumably because scalar Θ cannot decouple a polynomial syzygy module. Our M2 port replicates this faithfully in `MatrixHi.m2` / `MartyushevClean.m2`, but the package wiring routes `templateSolve` through `getH0` — which has no MatrixHi branch — so the working primitive is not yet user-facing via `templateSolve`.**
