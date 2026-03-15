# Greedy Strategy Plan (Simplified Row-Wise Only)

This note describes the **next implementation plan** for `Strategy => "Greedy"`.

## Goals

1. Use only **row-wise greedy** zeroing of `W` columns.
2. Remove the **column-wise greedy** assignment phase for now.
3. Avoid materializing full dense `H = H0e + Theta * H1e`.
4. Avoid scanning all columns each iteration; use a **candidate shortlist**.
5. Prioritize candidate columns with the **fewest free `θ` variables**.

## Keep from current pipeline

1. Compute base `H0` as before from Groebner data.
2. Compute `H1 := transpose sub(syz(gens J), ring J)`.
3. Keep affine constraint solving style used by current row-wise enforcement:
   - affine-linear only
   - single active variable per equation
   - consistency checks for previously assigned `θ`

## Planned structural change

Instead of building dense `Theta`, `H1e`, and `H = H0e + Theta*H1e`:

1. Build a lazy affine representation for each needed equation:
   - constant part from `H0`
   - linear `θ` contributions from `H1`
2. Build `W` constraints directly from this representation.
3. Evaluate/update only touched equations when assignments are added.

This keeps the same math objective while avoiding the expensive dense multiply.

## Row-wise greedy loop (new policy)

At each iteration:

1. Build/refresh a candidate set of nonzero columns.
2. Score columns by:
   - primary: smaller number of free `θ` variables in its active equations
   - secondary: estimated zeroing gain (or current violation count)
3. Keep only top `K` candidates (`shortlist`).
4. Try zeroing columns only from this shortlist.
5. Commit the best improving assignment set.
6. Repeat until no improvement.

## Candidate shortlist details

Suggested defaults:

1. `K = min(50, ceil(0.1 * #activeColumns))` as a starting rule.
2. Recompute shortlist every iteration (or every few iterations if needed).
3. Tie-breakers:
   - fewer equations remaining in the column
   - lower total symbolic complexity
   - deterministic column index

## Temporarily removed step

- No `columnWiseGreedyAssignments` pass.
- No computation of excessive monomials for greedy decision-making.

## Expected benefits

1. Removes the main bottleneck (`Theta * H1e`).
2. Reduces per-iteration search cost by avoiding full-column scans.
3. Keeps implementation simpler while preserving a working greedy baseline.

## Implementation checklist

1. Refactor `getH0(..., Strategy => "Greedy")` to row-wise-only flow.
2. Remove/disable column-wise invocation and scoring in the Greedy branch.
3. Introduce shortlist builder and scoring function (fewest free `θ` first).
4. Add timing logs for:
   - shortlist build
   - row-wise trial phase
   - assignment application
5. Verify output orientation remains compatible with downstream template code.
