# Efficient pipeline for implementing §4.1 (parameter adjustment) without forming large symbolic matrices

This outlines a concrete implementation strategy for the paper’s §4.1 idea: choose parameters **θ** to make many **columns of W** identically zero, thereby shrinking the elimination template — **without** explicitly forming `H = H0 + Θ·H1` via dense multiplication or constructing every `Z_k`.

---

## 0) Notation (matching the paper)

- `H = H0 + Θ H1`, where:
  - `H0 ∈ R^{d×s}` and `H1 ∈ R^{l×s}` are polynomial matrices
  - `Θ ∈ R^{d×l}` contains adjustable parameters `θ_{i,j}`
- `h_k` is column `k` of `H`
- For each column, `h_k = Z_k c_k` and `W = [Z_1 ... Z_s]`
- Goal (§4.1): pick θ to maximize the number of **zero columns of W**, i.e. solve affine constraints `w_k(θ)=0` for as many `k` as possible.

Key structural property used for speed:
- **Row-local parameters**: row `i` only uses parameters `θ_{i,*}`.  
  Therefore, `w_k = 0` splits into **d independent scalar equations** (one per row) in the **l parameters of that row**.

---

## 1) Data representations (use sparse everywhere)

### 1.1 Polynomial representation
Use a sparse dictionary / hash map:
- `poly = { monomial_id : coeff }`
- monomial_id can be:
  - packed exponent vector, or
  - integer ID from a monomial ordering table.

### 1.2 Sparse row representation for `H1`
Store each row `r_j` of `H1` as:
- `H1_rows[j] = [(k, poly), (k2, poly2), ...]`  (only nonzero columns)

### 1.3 Sparse row representation for `H0`
Similarly:
- `H0_rows[i] = [(k, poly), ...]`

This alone avoids the slow dense multiplication `Θ·H1`.

---

## 2) Avoid forming `H = H0 + Θ·H1` (row-wise accumulation)

Instead of matrix multiplication, compute row `i` of `H` as:

`H_row(i) = H0_row(i) + Σ_{j in active(i)} θ_{i,j} * H1_row(j)`

Implementation:
1. Maintain `active(i)` = indices `j` where θ_{i,j} is currently nonzero.
2. For each active `j`, iterate only the sparse entries of `H1_rows[j]` and accumulate into row `i`.

Complexity becomes proportional to the number of nonzeros touched, not `d*l*s`.

---

## 3) Don’t build `Z_k` explicitly: build **affine constraints** for columns of W

Greedy selection only needs to know whether a given column of `W` can be forced to zero, and which θ-values do it.

Represent each relevant scalar constraint as an affine form:

`w = a + Σ_{j=1..l} θ_{i,j} b_j`

Because of row-locality, for each W-column `k` you store **d independent affine forms**, one per row:

`w_{k,i}(θ_i) = a_{k,i} + Σ_j θ_{i,j} b_{k,i,j}`

### 3.1 Data structure for affine forms
For each `(k, i)` pair (W-column k, row i):

- `a`: constant polynomial/scalar term
- `b`: sparse map `{ j : coeff }` for θ_{i,j}

Store only nonzero entries:
- `Affine[(k,i)] = (a, {j: b_j})`

Also store an incidence list for fast updates:
- `Incidence[(i,j)] = list of (k) where b_{k,i,j} ≠ 0`

---

## 4) Offline precomputation (symbolic once per solver)

### 4.1 Build monomial tables and supports
- Enumerate monomials in your chosen basis/order.
- Precompute mappings needed for “shift” operations if you use them.

### 4.2 Build sparse `H0_rows`, `H1_rows`
- Construct from your polynomial system and chosen monomial basis.
- Keep strictly sparse.

### 4.3 Construct the affine constraint representation for W
Goal: build `Affine[(k,i)] = (a, b_map)` without materializing W densely.

How:
1. Identify which “W-columns” correspond to which shift/template rows in your implementation.
2. For each W-column `k`:
   - For each row `i` determine the symbolic dependence on θ_{i,*}`:
     - the constant part goes into `a_{k,i}`
     - coefficients in front of each θ_{i,j}` go into `b_{k,i,j}`

Store:
- `Affine[(k,i)]`
- `Incidence[(i,j)]`

> Tip: if your W construction is itself based on moving monomials/coefficients around, implement it as transformations on sparse maps, always emitting affine terms instead of full columns.

---

## 5) Greedy selection (fast, row-wise)

### 5.1 What “zeroing a W-column” means computationally
To zero W-column `k` you must satisfy:
- for all rows `i`: `w_{k,i}(θ_i) = 0`

Since each row uses its own θ-block, solve per-row:

For each row i:
- Solve `a_{k,i} + Σ_j θ_{i,j} b_{k,i,j} = 0`

Common cases:
- If `b_map` has a single variable `θ_{i,j}`:
  - set `θ_{i,j} = -a / b`
- If multiple θ’s appear:
  - choose one free variable and solve (or use least-norm / pivot heuristic)
  - or declare “not solvable under current restrictions” if you constrain θ to {0,1} etc.

### 5.2 Scoring a candidate assignment
For each candidate decision (e.g., choose θ_{i,j}=0, or choose value solving a pivot equation):
- estimate how many W-columns become all-zero (i.e., all rows satisfy zero)
- use an incremental score `σ` (like the paper’s greedy strategy).

### 5.3 Incremental updates (avoid rescanning everything)
When you fix θ_{i,j}:
- Only affine forms that contain that θ change.
- Use `Incidence[(i,j)]` to find affected W-columns `k`.
- Update `a_{k,i}` ← `a_{k,i} + θ_{i,j} * b_{k,i,j}` and remove `j` from `b_map` if you fully substitute.

Maintain a boolean/score for each W-column `k`:
- `is_zero_row[(k,i)]` whether row constraint is currently satisfied
- `zero_count[k]` = number of rows satisfied
- W-column `k` is eliminated if `zero_count[k] == d`

All updates become local.

---

## 6) Final instantiation of H and template extraction (only after θ is chosen)

Once greedy finishes and θ is fixed:

1. Instantiate `H` row-by-row:
   - `H_row(i) = H0_row(i) + Σ_{j in active(i)} θ_{i,j} * H1_row(j)`
2. Extract the reduced elimination template:
   - drop rows/columns corresponding to eliminated W-columns / shifts
3. (Optional but recommended) apply §4.3 cleanup:
   - remove dependent rows
   - remove dependent excessive monomial columns

---

## 7) Complexity wins vs. naïve implementation

### Naïve
- Form ΘH1 as dense multiply: ~`O(d*l*s)` operations on large objects (polynomials)
- Build each `Z_k`: often `O(d*s)` per k → huge

### This pipeline
- Works in **sparse support** + **affine forms**
- Greedy updates touch only incidence neighborhoods:
  - `O(#affected_constraints)` per θ decision
- Final H build touches only active θ’s and nonzeros in relevant H1 rows

---

## 8) Minimal checklist for a working implementation

- [ ] Sparse polynomial type with cheap add / scale
- [ ] Sparse storage of H0 rows and H1 rows
- [ ] Construction of affine constraints `Affine[(k,i)]`
- [ ] Incidence lists `Incidence[(i,j)]`
- [ ] Greedy loop that:
  - proposes θ decisions
  - scores them using `zero_count`
  - applies decision with local updates only
- [ ] Final row-wise instantiation of H
- [ ] Optional §4.3 dependency pruning

---

## 9) Practical tips

- Prefer **integer / rational** arithmetic in greedy if possible (avoid floating noise).
- Cache polynomial simplifications (e.g., remove zero terms after updates).
- If constraints are scalar numeric (not polynomial) after substitution, store them as plain numbers for speed.
- If you restrict θ choices (e.g., θ ∈ {0,1}), encode that in the per-row solver for `w_{k,i}=0`.

---
