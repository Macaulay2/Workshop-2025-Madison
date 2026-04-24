# matrixHi explained with a tiny toy example

## The pieces in plain words

- **Residual monomial** — a monomial you get by multiplying the action variable by a basis element and landing *outside* the basis.
- **Gap polynomial** — the residual monomial minus its normal form. By construction this lies in the ideal, so it can be written as a combination of the generators.
- **Undetermined $\alpha$** — a placeholder scalar coefficient in front of each monomial that appears inside the cofactors $H_{i,j}$. We set up a linear system in these $\alpha$'s, one equation per monomial on both sides, and solve.

Let me unfold this on a baby system you can compute by hand.

---

## Setup

Take the ring $R = \mathbb{Q}[x, y]$ and the ideal

$$
J = \langle F_1, F_2 \rangle, \qquad F_1 = x^2 + y^2 - 1, \quad F_2 = 2x - y.
$$

Geometrically: unit circle intersected with the line $y = 2x$. Two points: $\pm(\tfrac{1}{\sqrt{5}}, \tfrac{2}{\sqrt{5}})$.

### Quotient basis

In lex with $y > x$, the Groebner basis of $J$ is $\{\,2x - y,\; 5x^2 - 1\,\}$.
Leading monomials are $\{y, x^2\}$, so the quotient basis is the monomials *not* divisible by them:

$$
B = \{\,1,\; x\,\}, \qquad \deg J = |B| = 2.
$$

Call them $b_1 = 1$, $b_2 = x$.

### Action variable

Let the action variable be $a = x$. We want to understand multiplication by $x$ on the quotient ring $R/J$ in terms of the basis $B$.

---

## Step 1 — residuals

For each basis element $b$, compute $a \cdot b$ and check: is it still in the basis?

| $b_k$ | $a \cdot b_k$ | in $B$? | classification |
|---|---|---|---|
| $b_1 = 1$ | $x$ | **yes** ($= b_2$) | recyclable |
| $b_2 = x$ | $x^2$ | **no** | **residual** |

**Residual set:** $\{x^2\}$. One residual.

A residual is a "new" monomial: multiplying by $a$ kicked us out of the basis.

---

## Step 2 — gap polynomials

For each residual $r$, the gap polynomial is

$$
G = r - \mathrm{NF}(r, J)
$$

where $\mathrm{NF}(r, J)$ means "reduce $r$ modulo the Groebner basis, express in terms of $B$."

Here $r = x^2$, and $5x^2 - 1 \in J$ tells us $x^2 \equiv \tfrac{1}{5} \pmod J$. So

$$
G_1 = x^2 - \tfrac{1}{5}.
$$

**Why "gap":** $G_1$ lies *in* $J$ (subtracting the normal form always lands you in the ideal), but the term $x^2$ sticks out beyond the basis. The polynomial $G_1$ is literally the gap between the residual $x^2$ and its basis-level representation $\tfrac{1}{5}$.

Since $G_1 \in J$, there must exist polynomials $H_{1,1}, H_{1,2}$ with

$$
G_1 \;=\; H_{1,1}\!\cdot\! F_1 \;+\; H_{1,2}\!\cdot\! F_2.
$$

`matrixHi`'s job is to find such $H_{1,1}, H_{1,2}$ — and ideally to find them with *as few distinct monomials as possible*, because every monomial of $H_{1,j}$ eventually becomes a row of the template.

---

## Step 3 — undetermined $\alpha$'s

We don't guess $H_{1,1}$ and $H_{1,2}$. We **parameterize** them.

### Degree budget

- $\deg G_1 = 2$, $\deg F_1 = 2$, $\deg F_2 = 1$.
- Minimum degree needed for $H_{1,j}$: $d_j^{\min} = \max(0, \deg G_1 - \deg F_j)$.
  - $d_1^{\min} = 0$ → $H_{1,1}$ is a constant.
  - $d_2^{\min} = 1$ → $H_{1,2}$ has monomials $\{1, x, y\}$.

### Write $H$ with $\alpha$ placeholders

$$
H_{1,1} = \alpha_0, \qquad
H_{1,2} = \alpha_1 + \alpha_2\, x + \alpha_3\, y.
$$

The four $\alpha_0, \alpha_1, \alpha_2, \alpha_3 \in \mathbb{Q}$ are *undetermined scalars* — one per (entry, monomial) pair. These are the unknowns.

### The linear system

Expand $H_{1,1} F_1 + H_{1,2} F_2$ and compare monomial-by-monomial with $G_1 = x^2 - \tfrac{1}{5}$.

$$
\alpha_0 (x^2 + y^2 - 1) + (\alpha_1 + \alpha_2 x + \alpha_3 y)(2x - y) \;=\; x^2 - \tfrac{1}{5}.
$$

Multiply out and collect:

| monomial | LHS coefficient | RHS coefficient |
|---|---|---|
| $x^2$ | $\alpha_0 + 2\alpha_2$ | $1$ |
| $y^2$ | $\alpha_0 - \alpha_3$ | $0$ |
| $xy$ | $-\alpha_2 + 2\alpha_3$ | $0$ |
| $x$ | $2\alpha_1$ | $0$ |
| $y$ | $-\alpha_1$ | $0$ |
| $1$ | $-\alpha_0$ | $-\tfrac{1}{5}$ |

Six equations, four unknowns. In matrix form $A\vec{\alpha} = \vec{b}$:

$$
\underbrace{\begin{pmatrix}
1 & 0 & 2 & 0 \\
1 & 0 & 0 & -1 \\
0 & 0 & -1 & 2 \\
0 & 2 & 0 & 0 \\
0 & -1 & 0 & 0 \\
-1 & 0 & 0 & 0
\end{pmatrix}}_{A}
\begin{pmatrix}\alpha_0 \\ \alpha_1 \\ \alpha_2 \\ \alpha_3\end{pmatrix}
=
\begin{pmatrix}1 \\ 0 \\ 0 \\ 0 \\ 0 \\ -\tfrac{1}{5}\end{pmatrix}.
$$

### Solve via RREF

$A$ has rank 4, $(A \mid \vec{b})$ has rank 4 → consistent, unique solution. From the last row: $\alpha_0 = \tfrac{1}{5}$. From rows 4–5: $\alpha_1 = 0$. Then row 3 gives $\alpha_2 = 2\alpha_3$; row 1 gives $\tfrac{1}{5} + 2\alpha_2 = 1$ so $\alpha_2 = \tfrac{2}{5}$, $\alpha_3 = \tfrac{1}{5}$.

**All $\alpha$ determined. Zero free $\alpha$.** The unique cofactors are

$$
H_{1,1} = \tfrac{1}{5}, \qquad H_{1,2} = \tfrac{2}{5} x + \tfrac{1}{5} y.
$$

Check:
$$
\tfrac{1}{5}(x^2 + y^2 - 1) + (\tfrac{2}{5}x + \tfrac{1}{5}y)(2x - y) = \tfrac{1}{5}(5x^2 - 1) = x^2 - \tfrac{1}{5}. \checkmark
$$

---

## Step 4 — from $H$ to the template

The *shift monomials* per generator are the monomials appearing in column $j$ of $H$:

- Shifts for $F_1$: $\{1\}$ (from $H_{1,1} = \tfrac{1}{5}$).
- Shifts for $F_2$: $\{x, y\}$ (from $H_{1,2} = \tfrac{2}{5}x + \tfrac{1}{5}y$).

The *shifted polynomials* are $\{1 \cdot F_1,\; x \cdot F_2,\; y \cdot F_2\}$. Each becomes a row of the template; each distinct monomial appearing in their combined support becomes a column.

**Template:** 3 rows × 5 columns, with each row encoding $F_1$, $xF_2$, or $yF_2$ in terms of the shared monomial basis $\{x^2, xy, y^2, x, y, 1\}$ (before dropping duplicates).

*This is what `buildTemplateFromH` constructs.*

---

## When do free $\alpha$'s appear?

Free $\alpha$'s show up when $A$ has a nontrivial null space. Concretely, this happens when:

1. **There are more generators than needed** (an overdetermined system). Extra generators create syzygies — multiple ways to write the same gap polynomial.
2. **The degree budget is overly generous.** Bumping the degree of $H_{1,j}$ gives more $\alpha$'s than the equations can fix.

### Overdetermined variant

Add a redundant generator $F_3 = 4x^2 - y^2$ (note $F_3 = (2x+y) \cdot F_2$, so $F_3 \in J$ already). Now

$$
J' = \langle F_1, F_2, F_3 \rangle
$$

has the same solutions, same basis, same residual, same $G_1$. But the parameterization grows:

$$
H_{1,1} = \alpha_0, \qquad
H_{1,2} = \alpha_1 + \alpha_2 x + \alpha_3 y, \qquad
H_{1,3} = \alpha_4.
$$

Now 5 unknowns. Redo the linear system — you'll find the equations force

$$
\alpha_0 = \tfrac{1}{5},\quad \alpha_1 = 0,\quad \alpha_2 = 2\alpha_3,\quad \alpha_4 = \tfrac{1}{5} - \alpha_3,
$$

and $\alpha_3$ is **free**.

Different choices of $\alpha_3$ give different valid cofactors:

| $\alpha_3$ | $H_{1,1}$ | $H_{1,2}$ | $H_{1,3}$ | shifts for $F_1, F_2, F_3$ | template rows |
|---|---|---|---|---|---|
| $0$ | $\tfrac{1}{5}$ | $0$ | $\tfrac{1}{5}$ | $\{1\}, \emptyset, \{1\}$ | **2** |
| $\tfrac{1}{5}$ | $\tfrac{1}{5}$ | $\tfrac{2}{5}x + \tfrac{1}{5}y$ | $0$ | $\{1\}, \{x,y\}, \emptyset$ | 3 |
| $\tfrac{1}{10}$ | $\tfrac{1}{5}$ | $\tfrac{1}{5}x + \tfrac{1}{10}y$ | $\tfrac{1}{10}$ | $\{1\}, \{x,y\}, \{1\}$ | 4 |

**Setting $\alpha_3 = 0$ gives the smallest template** (2 rows). This is exactly what the "particular solution" that `matrixHi` returns does — it sets all free $\alpha$ to zero, which in this case is the optimum.

The `adjustParams` greedy is for cases where the particular solution is *not* optimal — for example, when different excessive monomials are controlled by different $\alpha$'s and you need to co-optimize them. That's the hard case on large problems like P4P+fr where 1095 $\alpha$ are free and the search space is enormous.

---

## Summary of vocabulary

| Term | Meaning | In the toy example |
|---|---|---|
| **Residual** $r$ | $a \cdot b \notin B$ | $x^2$ |
| **Gap polynomial** $G$ | $r - \mathrm{NF}(r, J) \in J$ | $x^2 - \tfrac{1}{5}$ |
| **Entry monomial basis** $\mathrm{entMons}[j]$ | monomial basis of $H_{i,j}$ at its degree budget | $\{1\}$ for $F_1$; $\{1,x,y\}$ for $F_2$ |
| **Undetermined $\alpha$** | scalar unknown, one per (generator $j$, monomial $m$) pair | $\alpha_0, \alpha_1, \alpha_2, \alpha_3$ |
| **Linear system** | one equation per monomial of $G$ or expanded support of $\sum \alpha_{j,m} \cdot m \cdot F_j$ | 6 equations above |
| **Particular solution** | solve the system, set all free $\alpha$ to $0$ | $H_{1,1} = \tfrac{1}{5}, H_{1,2} = \tfrac{2x+y}{5}$ |
| **Free $\alpha$** | null-space direction — room to maneuver | $0$ in base example, $1$ in overdetermined variant |
| **Shift monomial** | monomial appearing in column $j$ of $H$ — shifts $F_j$ into the template | $\{1\}$ for $F_1$, $\{x,y\}$ for $F_2$ |
| **Excessive monomial** | a monomial appearing in $\{m \cdot F_j\}$ but not in $R \cup B$ | after shifting, depends on problem |

---

## Back to the big picture

In the *real* problems (5pt essential, PnP, etc.):

- **nG** (number of gap polynomials) is dozens to hundreds.
- The linear system per row is large but solvable; `matrixHi` bumps $d$ automatically when infeasible.
- The particular solution ($\alpha_\text{free} = 0$) already matches the paper's row count on simple cases (5pt: 10×20, stitching: 48×66).
- On harder cases (P4P+fr: 1095 free $\alpha$), you need `adjustParams` to pick non-zero $\alpha$'s that cancel excessive monomials — and that greedy is the remaining research gap.

So in one sentence: **`matrixHi` treats each row of $H$ as a fresh linear-algebra problem in scalar $\alpha$'s, solves it over the field $\mathbb{F}_p$ (or $\mathbb{Q}$), and the free dimensions of its null space are the handles the greedy pulls on.**

---

# What `getH0` does — same toy, different object

Same ring, same ideal, same action variable:

$$
R = \mathbb{Q}[x,y],\quad J = \langle F_1, F_2\rangle,\quad F_1 = x^2+y^2-1,\quad F_2 = 2x-y,\quad B=\{1,x\},\quad a=x.
$$

`getH0` answers a *different* question than `matrixHi`. It packages multiplication by $a$ on $R/J$ — for **every** basis element, trivial ones included — into a single matrix $H_0$.

## What $H_0$ represents

Define the "action-polynomial" row vector $V$ whose $k$-th entry is

$$
V_k \;=\; a \cdot b_k \;-\; \bigl(\text{the expression of } a \cdot b_k \text{ in terms of } B\bigr).
$$

In other words, $V_k$ is **zero if $a \cdot b_k$ is already in the basis** (recyclable), and otherwise equals "residual minus its normal form" — which is exactly the gap polynomial of matrixHi.

Then $H_0$ is chosen so that

$$
V \;=\; \bigl[F_1 \;\; F_2\bigr] \cdot H_0.
$$

That is, each column of $H_0$ gives cofactors for one basis element's action.

## Compute $V$ for the toy

| $k$ | $b_k$ | $a \cdot b_k$ | in $B$? | $\mathrm{NF}$ in $B$ | $V_k$ |
|---|---|---|---|---|---|
| 1 | $1$ | $x$ | yes ($=b_2$) | $x$ | $0$ |
| 2 | $x$ | $x^2$ | no | $\tfrac{1}{5}$ | $x^2 - \tfrac{1}{5}$ |

So $V = \bigl[\,0,\; x^2 - \tfrac{1}{5}\,\bigr]$.

## Compute $H_0$

We want $\bigl[F_1 \;\; F_2\bigr] \cdot H_0 = V$. $H_0$ is $n_\text{gens} \times n_\text{basis} = 2 \times 2$. Column by column:

- **Column 1:** $F_1 H_0[1,1] + F_2 H_0[2,1] = 0$. Trivial choice: $H_0[\,\cdot\,, 1] = (0, 0)$.
- **Column 2:** $F_1 H_0[1,2] + F_2 H_0[2,2] = x^2 - \tfrac{1}{5}$. *This is exactly the linear system matrixHi solved.* One answer: $H_0[1,2] = \tfrac{1}{5}$, $H_0[2,2] = \tfrac{2x+y}{5}$.

$$
H_0 = \begin{pmatrix} 0 & \tfrac{1}{5} \\[2pt] 0 & \tfrac{2x+y}{5} \end{pmatrix}.
$$

## How `getH0` actually produces this in code

Not by solving a linear system in $\alpha$'s. `getH0` uses the **Groebner basis change matrix** machinery *(EliminationTemplates.m2:94-112)*:

```
G := gb(F, ChangeMatrix => true)        -- reduced GB + matrix HGF with F*HGF = G
P := last coefficients(B % F)            -- change of basis: B in quotient coords
V := a*B - lift(a*B * P^{-1}, S/F) * P   -- the action-polynomial vector
HVG := V // gens G                       -- polynomial division of V by GB
H0  := HGF * HVG                         -- compose: F * H0 = V
```

So `H0` is whatever cofactors fall out of GB-reduction. There's **one particular solution**; no $\alpha$'s, no null-space structure, no user-visible freedom.

## Side-by-side comparison

| | `getH0` | `matrixHi` |
|---|---|---|
| **Output shape** | $n_\text{gens} \times n_\text{basis}$ = $2 \times 2$ | $n_\text{gap} \times n_\text{gens}$ = $1 \times 2$ |
| **Covers** | every basis element (trivial cols are $0$) | only non-trivial residuals |
| **Construction** | GB change matrix + polynomial division | linear system in undetermined $\alpha$, RREF |
| **Degrees of entries** | whatever GB reduction produces (often higher than needed) | minimum feasible, auto-bumped only if infeasible |
| **Freedom exposed** | none — one particular solution | pivot / free $\alpha$ structure |
| **Handle for greedy** | only through $H_0 + \Theta \cdot H_1$ (next section) | directly in the free $\alpha$'s |

In the toy, the non-trivial column of $H_0$ agrees with the row of `matrixHi`. On larger problems they generally don't agree — GB reduction tends to produce higher-degree cofactors than the min-degree linear system.

## The three `getH0` strategies, same toy

### Default: `Strategy => null`

Return the $H_0$ above as-is. Shifts per generator come from the column's monomial support:

- Row 1 of $H_0$ (cofactors of $F_1$): $\{\tfrac{1}{5}\}$ → shift $\{1\}$.
- Row 2 of $H_0$ (cofactors of $F_2$): $\{\tfrac{2x+y}{5}\}$ → shifts $\{x, y\}$.

Template: $F_1,\ xF_2,\ yF_2$ plus the graph-ideal row for the action variable = 3 rows + basis rows.

### `Strategy => "Larsson"`

Project $H_0$ onto a smaller representative modulo the syzygy module:

$$
H_0^{\text{Larsson}} = H_0 \bmod \mathrm{image}(\mathrm{syz}(\text{gens}\,J)).
$$

In the toy, $\mathrm{syz}(F_1, F_2)$ is trivial (generated by high-degree relations), so Larsson's output equals default.

On larger problems with many syzygies, this projection strictly shrinks the cofactor support — e.g. it knocks several rows off the P4P template.

### `Strategy => "MatrixHi"`

Build $H$ via the per-row linear system in undetermined scalar $\alpha$'s (one per `(entry, monomial)` pair), then evaluate all free $\alpha$'s at 0. This is the particular solution from `buildHSymbolic`. On problems where the linear system has no free $\alpha$'s it is already optimal.

### `Strategy => "Greedy"` (Martyushev CVPR 2022)

Same as `MatrixHi` through the `buildHSymbolic` step, but keep the free $\alpha$'s symbolic and then run `adjustParams` (Martyushev CVPR 2022 §4): iterate through excessive monomials and commit $\alpha$ values that cancel them. On a problem with 0 free $\alpha$'s this step is a no-op and output matches `MatrixHi` exactly.

### Retired strategies (see git log pre-2026-04-16 if you need them)

- A $\Theta$-space greedy over the Larsson parameterization $H = H_0 + \Theta \cdot H_1$ with $H_1 = \mathrm{syz}(F)^T$. Scalar $\theta$'s but polynomial coupling through syzygy rows; structurally incapable of matching paper sizes on problems with nontrivial syzygy coupling.
- A gap-polynomial $H_0$ construction via Gröbner division (`qk := matrix{{G_k}} // gens J`). Matched paper row counts on small problems but exposed no freedom for further greedy search, so `adjustParams` had nothing to work with.

## The big picture — two parametrizations compared

| | `matrixHi` (Family B) | `getH0` "Greedy" (Family A) |
|---|---|---|
| Parametrization | $H_{i,j} = \sum_m \alpha_{i,j,m} \cdot m$ | $H = H_0 + \Theta H_1$, $\Theta$ scalar |
| What moves when I set one parameter | just the coefficient of one monomial in one entry | coefficients of *every* monomial in one syzygy row |
| Coupling | none — free $\alpha$'s are linearly independent | polynomial (syzygy rows are polynomials) |
| Free parameters on P4P+fr | 1095 | 0 |
| Room for the greedy to actually shrink the template | yes | no |

This is why `matrixHi` is the right primitive for closing the gap to Maple's sizes, and why `getH0` + the syzygy-$\Theta$ greedy (the original pre-port implementation) couldn't get there no matter how clever the search became.

## TL;DR

- `getH0` gives you **the action-polynomial cofactors for every basis element**, packaged as an $n_\text{gens} \times n_\text{basis}$ matrix built from the GB change matrix. It's a single particular solution.
- `matrixHi` gives you **the cofactors only for the non-trivial residuals**, row-by-row, with each row built from a fresh linear system in undetermined scalar $\alpha$'s that exposes the null space.
- Every strategy — Default, Larsson, MatrixHi-mode (`Strategy => "Greedy", AdjustParams => false`), Greedy — skips the graph-ideal extension on monomial actions and builds the template on the `[E | R | B]` layout, saving $|B|$ rows. Polynomial actions lift to $R_s = R[s]/⟨s-a⟩$ and take the same path.
- On the toy, non-trivial columns of $H_0$ coincide with rows of `matrixHi` (both solve the same gap equation). On real problems they diverge: `matrixHi`'s min-degree linear system is strictly more flexible than `getH0`'s GB-dictated particular solution.
- `matrixHi`'s Family B parametrization is exposed via `Strategy => "Greedy"` (adjustParams on) / `AdjustParams => false` (α:=0 particular solution); `getH0`'s Family A is the Default / Larsson path. There is no separate "MatrixHi" strategy identifier anymore — it is a mode of Greedy.
