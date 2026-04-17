-- MartyushevClean.m2: Fresh, minimal port of Martyushev's matrixHi + adjustParams.
-- Rewritten from scratch to avoid accumulated bugs from MartyushevExact.m2.
--
-- Design:
-- - Work in flat ring R = FF[vars..., alphas...] from the start
-- - Use mutable hash table for alpha values (0 = unassigned)
-- - adjustParams substitutes into H directly without intermediate decompositions
-- - All operations use explicit lists (avoid Sequence bugs)
--
-- Loaded two ways:
-- (a) from inside EliminationTemplates.m2 newPackage (no needsPackage — that
--     would be a circular dep).
-- (b) by external tests that `debug needsPackage "EliminationTemplates"`
--     first, then `load "EliminationTemplates/MartyushevClean.m2"`.

-- ============================================================================
-- buildH: main function combining matrixHi + adjustParams + constructTemplate
--
-- Input: action variable, basis B, ideal J
-- Output: (Hoptimized, residualMons, Blist, alphaEnv, extendedRing)
-- where Hoptimized is an nG x nF matrix with symbolic alpha entries that
-- have been partially reduced by the greedy.
-- ============================================================================
buildGapPolys = method()
buildGapPolys(RingElement, Matrix, Ideal) := (aVar, B, J) -> (
    R := ring J;
    aVar2 := sub(aVar, R);
    Blist := flatten entries lift(B, R);
    Bset := set Blist;
    aBlist := apply(Blist, b -> aVar2 * b);
    resMons := select(aBlist, m -> not Bset#?m);
    gapP := apply(resMons, r -> r - (r % J));
    (gapP, toList resMons, Blist)
);

-- Build the symbolic H in a FLAT extended ring.
-- Returns (Hsym, Rext, alphaVars, perRowData) where:
--   Rext = FF[base vars..., alpha_0, ..., alpha_{n-1}]
--   Hsym is nG x nF over Rext
--   alphaVars = list of alpha ring elements
--   perRowData[i] = (entMons, alphaRange) for row i
--
-- IMPORTANT: This handles the degree increment correctly — if the minimum
-- degree H[i,j] doesn't admit a solution, it incrementally raises d until
-- the linear system becomes feasible. This matches Maple's matrixHi.
buildHSymbolic = method(Options => {Verbose => false})
buildHSymbolic(List, List) := o -> (F, gapPolys) -> (
    R := ring first F;
    FF := coefficientRing R;
    nF := #F;
    nG := #gapPolys;
    degF := apply(F, f -> first degree f);

    -- For each row, find the right degree and solve
    rowResults := apply(nG, i -> (
        dG := first degree gapPolys#i;
        -- Try increasing d until feasible (same as Maple's matrixHi)
        d := 0;
        dMinMet := 0;
        chosenEm := null;
        chosenRref := null;
        chosenPivots := null;
        chosenNParams := 0;
        while d - dMinMet <= 0 do (
            em := apply(nF, j -> flatten entries basis(0, max(0, dG - degF#j) + d, R));
            nParams := sum apply(em, m -> #m);
            targetMons := unique flatten (
                flatten apply(nF, j -> flatten apply(em#j, mk -> flatten entries monomials(mk * F#j)))
                | flatten entries monomials(gapPolys#i)
            );
            Amat := matrix(FF, apply(targetMons, m -> (
                flatten apply(nF, j -> apply(em#j, mk -> coefficient(m, mk * F#j)))
            )));
            bvec := matrix(FF, apply(targetMons, m -> {coefficient(m, gapPolys#i)}));
            -- Check feasibility: rank(A) == rank(A|b)
            rkA := rank Amat;
            rkAug := rank(Amat | bvec);
            if rkAug > rkA then (
                -- Infeasible at this degree
                dMinMet = d + 1;
                d = d + 1;
                if d > 30 then error ("matrixHiRow: degree too high");
            ) else (
                chosenEm = em;
                chosenNParams = nParams;
                chosenRref = reducedRowEchelonForm(Amat | bvec);
                break;
            );
        );
        em := chosenEm;
        nParams := chosenNParams;
        rref := chosenRref;
        -- Identify pivots
        pivotCols := {};
        curR := 0;
        for col from 0 to nParams - 1 do (
            if curR >= numRows rref then break;
            if rref_(curR, col) == 1_FF then (
                pivotCols = append(pivotCols, col);
                curR = curR + 1;
            );
        );
        pivotSet := set pivotCols;
        freeCols := toList select(0..nParams-1, c -> not pivotSet#?c);
        -- For each pivot, store: (pivot col, pivot row, list of (freeColIdx, -coef))
        pivotInfo := apply(pivotCols, col -> (
            rowIdx := position(0..numRows rref - 1, r -> rref_(r, col) == 1_FF);
            constant := rref_(rowIdx, nParams);
            deps := apply(freeCols, fc -> (fc, -rref_(rowIdx, fc)));
            (col, constant, deps)
        ));
        (em, nParams, pivotInfo, freeCols)
    ));

    -- Total free params across all rows
    totalFree := sum apply(rowResults, r -> #(r#3));
    if o.Verbose then << "buildHSymbolic: totalFree=" << totalFree << endl;

    -- Build flat extended ring
    baseNames := apply(numgens R, i -> toString R_i);
    alphaNames := apply(totalFree, k -> "aa" | toString k);
    Rext := if totalFree > 0 then
        FF[(baseNames | alphaNames) / getSymbol, MonomialOrder => {numgens R, totalFree}]
    else R;
    toRext := if totalFree > 0 then
        map(Rext, R, apply(numgens R, i -> Rext_i))
    else map(R, R);
    alphaVars := if totalFree > 0 then apply(totalFree, k -> Rext_(numgens R + k)) else {};

    -- Build Hsym row by row
    alphaOffset := 0;
    perRowData := new MutableList from apply(nG, i -> null);
    Hrows := apply(nG, i -> (
        (em, nParams, pivotInfo, freeCols) := rowResults#i;
        nFreeRow := #freeCols;
        rowAlphas := apply(nFreeRow, k -> alphaVars#(alphaOffset + k));

        -- Mapping from free col idx in this row → global alpha var
        freeIdxToAlpha := hashTable apply(nFreeRow, k -> freeCols#k => rowAlphas#k);

        -- Build alpha value expression for each param idx 0..nParams-1
        alphaExpr := new MutableList from apply(nParams, p -> 0_Rext);
        -- Set free cols to their corresponding alpha var
        scan(nFreeRow, k -> (
            alphaExpr#(freeCols#k) = rowAlphas#k;
        ));
        -- Set pivot cols to their expression: constant - sum (coef * alpha)
        scan(pivotInfo, pi -> (
            (col, constant, deps) := pi;
            val := promote(constant, Rext);
            scan(deps, dep -> (
                (fc, coef) := dep;
                if freeIdxToAlpha#?fc then
                    val = val + promote(coef, Rext) * freeIdxToAlpha#fc;
            ));
            alphaExpr#col = val;
        ));

        -- Build the row: for each gen j, sum (alphaExpr#p * entMon) over monomials of entMons#j
        rowPolys := apply(nF, j -> (
            pStart := sum apply(j, jj -> #(em#jj));
            sum apply(#(em#j), k -> alphaExpr#(pStart + k) * toRext(em#j#k))
        ));

        perRowData#i = (em, toList alphaExpr, rowAlphas);
        alphaOffset = alphaOffset + nFreeRow;
        rowPolys
    ));

    Hsym := matrix(Rext, Hrows);
    (Hsym, Rext, alphaVars, toList perRowData)
);

-- ============================================================================
-- adjustParams: a simpler, direct implementation of Martyushev CVPR 2022 §4.
--
-- Strategy: For each excessive monomial, directly collect the per-entry
-- coefficient expressions that contribute to it, build a linear system,
-- solve, and substitute into Hsym.
-- ============================================================================
-- Note: 5-arg method bundled as (F, Hsym, Rext, (alphaVars, RB))
adjustParams = method(Options => {Verbose => false, MaxIter => 1000})
adjustParams(List, Matrix, Ring, Sequence) := o -> (F, Hsym, Rext, pair) -> (
    (alphaVars, RB) := pair;
    R := ring first F;
    FF := coefficientRing R;
    nF := #F;
    nG := numRows Hsym;
    baseVars := apply(numgens R, i -> Rext_i);
    Fext := apply(F, f -> sub(f, Rext));
    RBset := set apply(RB, m -> sub(m, Rext));

    Hcurrent := mutableMatrix Hsym;

    -- Get current shift monomials per generator (in vars only)
    getShiftMons := (HM) -> (
        apply(nF, j -> (
            unique flatten apply(nG, i -> (
                e := HM_(i,j);
                if e == 0 then {} else (
                    (mm, cc) := coefficients(e, Variables => baseVars);
                    flatten entries mm
                )
            ))
        ))
    );

    -- Get excessive monomials
    getExcessive := (HM) -> (
        sh := getShiftMons HM;
        allExp := unique flatten apply(nF, j -> (
            flatten apply(sh#j, m -> flatten entries monomials(m * Fext#j, Variables => baseVars))
        ));
        select(allExp, m -> not RBset#?m)
    );

    initialExcessive := getExcessive(matrix Hcurrent);
    if o.Verbose then
        << "adjustParams: " << #initialExcessive << " initial excessive" << endl;

    if #initialExcessive == 0 then return matrix Hcurrent;

    -- Process excessive monomials in order (smallest producing-count first)
    -- For each, try to zero it out by setting alphas
    iterations := 0;
    remainingExc := initialExcessive;

    while iterations < o.MaxIter and #remainingExc > 0 do (
        iterations = iterations + 1;
        improved := false;

        -- Sort by difficulty
        exWithCount := apply(remainingExc, e -> (
            cnt := sum apply(nF, j -> (
                sh := unique flatten apply(nG, i -> (
                    entry := Hcurrent_(i,j);
                    if entry == 0 then {} else (
                        (mm, cc) := coefficients(entry, Variables => baseVars);
                        flatten entries mm
                    )
                ));
                #select(sh, m -> coefficient(e, m * Fext#j) != 0_Rext)
            ));
            (cnt, e)
        ));
        sortedEx := apply(sort exWithCount, p -> p#1);

        -- Try each in order
        tookStep := false;
        scan(sortedEx, e -> (
            if tookStep then return;

            -- Collect all entries of H that contribute to e (when expanded via F)
            -- For each gen j, find shift monomials m such that coefficient(e, m * F[j]) != 0
            -- Then for each row i, take coefficient(m, H[i,j]) — this is an expression in alphas
            -- that must be zero across the sum over rows

            -- The specific equation: sum_i [(total coeff of e contributed by m * F[j] row i)] = 0
            -- More precisely, for each (j, m) pair, for each row i, the contribution to e is
            --   coefficient(m, H[i,j]) * coefficient(e, m * F[j])
            -- Setting ALL row contributions to zero means: for each i with coefficient(e, m*F[j]) != 0,
            --   coefficient(m, H[i,j]) = 0
            -- (or we'd need the SUM to be zero, but that requires summing rows)

            -- Actually each row of H contributes a separate shifted polynomial to the template,
            -- so each row's contribution to e must be independently zero for e to not appear.
            -- So we get one equation per (i, j, m) tuple where m * F[j] contains e.

            eqs := flatten apply(nG, i -> (
                flatten apply(nF, j -> (
                    entry := Hcurrent_(i,j);
                    if entry == 0 then return {};
                    (mm, cc) := coefficients(entry, Variables => baseVars);
                    monList := flatten entries mm;
                    coefList := flatten entries cc;
                    -- For each monomial m in the entry, check if m * F[j] contains e
                    eqsForEntry := toList apply(#monList, k -> (
                        m := monList#k;
                        if coefficient(e, m * Fext#j) != 0_Rext then
                            coefList#k  -- the alpha-expression coefficient
                        else
                            null
                    ));
                    select(eqsForEntry, x -> x =!= null)
                ))
            ));

            -- Filter to nonzero
            nonzeroEqs := select(eqs, x -> x != 0_Rext);
            if #nonzeroEqs == 0 then return;

            -- Build linear system in alphaVars
            nAlphas := #alphaVars;
            sysRows := apply(nonzeroEqs, eq -> (
                -- Extract linear coefficients wrt alphaVars
                row := apply(nAlphas, k -> sub(coefficient(alphaVars#k, eq), FF));
                -- Constant term = eq with alphas set to 0
                constSub := map(Rext, Rext, baseVars | apply(nAlphas, k -> 0_Rext));
                constVal := sub(constSub eq, FF);
                row | {-constVal}
            ));

            -- Drop trivial all-zero rows
            sysRows = select(sysRows, r -> not all(r, x -> x == 0_FF));
            if #sysRows == 0 then (
                -- e is already zero, remove from list
                remainingExc = select(remainingExc, m -> m != e);
                improved = true;
                tookStep = true;
                return;
            );

            sysMat := matrix(FF, sysRows);
            coeffMat := sysMat_{0..nAlphas-1};
            rhsCol := sysMat_{nAlphas};

            -- Check consistency
            if rank(coeffMat | rhsCol) > rank coeffMat then return;  -- skip, not solvable

            -- Solve and extract assignment
            solRref := reducedRowEchelonForm(coeffMat | rhsCol);
            newAlphaVals := new MutableList from apply(nAlphas, k -> null);  -- null = unassigned
            curR := 0;
            for col from 0 to nAlphas - 1 do (
                if curR >= numRows solRref then break;
                if solRref_(curR, col) == 1_FF then (
                    newAlphaVals#col = sub(solRref_(curR, nAlphas), Rext);
                    curR = curR + 1;
                );
            );
            -- Set unassigned alphas to 0
            substImages := apply(nAlphas, k ->
                if newAlphaVals#k =!= null then newAlphaVals#k else 0_Rext);
            substMap := map(Rext, Rext, baseVars | substImages);

            -- Apply substitution to Hcurrent
            for i from 0 to nG - 1 do
                for j from 0 to nF - 1 do
                    Hcurrent_(i,j) = substMap(Hcurrent_(i,j));

            improved = true;
            tookStep = true;
            remainingExc = select(remainingExc, m -> m != e);
            if o.Verbose and iterations % 5 == 0 then (
                currExc := getExcessive(matrix Hcurrent);
                << "  iter " << iterations << ": eliminated e, " << #currExc << " remain" << endl;
            );
        ));

        if not improved then break;

        -- Recompute remaining excessive (some may have been zeroed as side effects)
        remainingExc = select(remainingExc, e -> (
            all(nF, j -> (
                sh := unique flatten apply(nG, i -> (
                    entry := Hcurrent_(i,j);
                    if entry == 0 then {} else (
                        (mm, cc) := coefficients(entry, Variables => baseVars);
                        flatten entries mm
                    )
                ));
                any(sh, m -> coefficient(e, m * Fext#j) != 0_Rext)
            ))
        ));
    );

    if o.Verbose then (
        finalExc := getExcessive(matrix Hcurrent);
        << "adjustParams: eliminated " << (#initialExcessive - #finalExc)
           << " of " << #initialExcessive << " excessive (in " << iterations << " iters)" << endl;
    );

    -- Project back to R (substituting all remaining alphas to 0)
    if #alphaVars > 0 then (
        nBase := numgens R;
        finalSubst := map(R, Rext, apply(nBase, i -> R_i) | apply(#alphaVars, k -> 0_R));
        matrix apply(nG, i -> apply(nF, j -> finalSubst(Hcurrent_(i,j))))
    ) else matrix Hcurrent
);

-- ============================================================================
-- extractActionFromTemplate: read the action matrix off a MatrixHi-style
-- template by RREF + pivot extraction.
--
-- Inputs:
--   M           — template matrix over the base field, with column blocks
--                 [excess | residual | basis], built by buildTemplateFromH.
--   resMonsList — residual monomials (a*b_i not in B), in the same column
--                 order as the residual block of M.
--   Blist       — quotient basis monomials, in R, in the same column order
--                 as the basis block of M.
--   aVar        — action variable.
--
-- Output: |B| x |B| matrix Ma over coefficientRing(R) such that, modulo I,
--   a * b_i = sum_k Ma_(k,i) * b_k.
--
-- For each b_i in B compute a*b_i. If a*b_i is in B, the i-th column of Ma
-- is the standard basis vector e_k. Otherwise a*b_i is some residual r_j;
-- we find the RREF pivot row for the j-th residual column, and read off the
-- B-block of that row (negated, since Rref row says: 1*r_j + sum b_k * (...)
-- + ... = 0, so r_j = -sum (...) b_k modulo I).
-- ============================================================================
extractActionFromTemplate = method()
extractActionFromTemplate(Matrix, List, List, RingElement) := (M, resMonsList, Blist, aVar) -> (
    R := ring first Blist;
    FF := coefficientRing R;
    a := sub(aVar, R);
    nR := #resMonsList;
    nB := #Blist;
    nE := numColumns M - nR - nB;
    Rref := reducedRowEchelonForm M;
    matrix(FF, apply(nB, i -> (
        bi := Blist#i;
        abi := a * bi;
        kInB := position(Blist, b -> b == abi);
        if kInB =!= null then (
            apply(nB, k -> if k == kInB then 1_FF else 0_FF)
        ) else (
            jInR := position(resMonsList, r -> r == abi);
            if jInR === null then
                error("extractActionFromTemplate: a*b_i not in R or B: " | toString abi);
            colR := nE + jInR;
            pivRow := null;
            for r from 0 to numRows Rref - 1 do (
                leading := position(0..numColumns Rref - 1, c -> Rref_(r,c) != 0_FF);
                if leading === colR then (pivRow = r; break);
            );
            if pivRow === null then
                error("extractActionFromTemplate: no pivot for residual col " | toString colR);
            apply(nB, k -> -Rref_(pivRow, nE + nR + k))
        )
    )))
)

-- ============================================================================
-- recoverSolutionsMatrixHi: given an action matrix Ma (in coefficientRing R)
-- and the basis monomials Blist (in R), eigendecompose over CC and return
-- the list of complex solutions, one per eigenvector.
--
-- Each eigenvector v is a vector of length |B|. Modulo I, v_k corresponds to
-- the value of b_k at one solution point (up to scale). We rescale so that
-- the slot for the monomial 1 equals 1, then read each variable's coordinate
-- from its position in Blist.
--
-- Output shape matches templateSolve / recoverSolutions: a List of solutions,
-- each a List of CC values, one per variable in declaration order of R.
-- ============================================================================
recoverSolutionsMatrixHi = method()
recoverSolutionsMatrixHi(Matrix, List, Ring) := (Ma, Blist, R) -> (
    TC := sub(Ma, CC);
    (evals, P) := eigenvectors TC;
    oneIdx := position(Blist, b -> b == 1_R);
    if oneIdx === null then
        error("recoverSolutionsMatrixHi: 1 not found in basis Blist");
    varPos := apply(numgens R, i -> position(Blist, b -> b == (gens R)#i));
    sols := {};
    for k from 0 to numColumns P - 1 do (
        v := P_{k};
        sc := v_(oneIdx, 0);
        if abs sc < 1e-12 then continue;
        v = (1/sc) * v;
        coords := apply(varPos, p -> if p =!= null then v_(p,0) else 0_CC);
        sols = append(sols, coords);
    );
    sols
)

-- ============================================================================
-- buildTemplateFromH: from H, F, and RB = residualMons | Blist, build the
-- elimination template matrix M over the base field. For each generator F[j],
-- the shift monomials come from column j of H. The template columns are
-- [ excessMons | RB ].
-- ============================================================================
buildTemplateFromH = method(Options => {Verbose => false})
buildTemplateFromH(Matrix, List, List) := o -> (Hmat, F, RB) -> (
    R := ring first F;
    FF := coefficientRing R;
    nF := #F;
    nG := numRows Hmat;

    sh := apply(nF, j ->
        toList set flatten apply(nG, i -> flatten entries monomials Hmat_(i,j))
    );
    xF := flatten apply(nF, j -> apply(sh#j, m -> m * F#j));
    V0 := toList set flatten apply(xF, p -> flatten entries monomials p);
    RBset := set RB;
    Eexcess := rsort select(V0, m -> not RBset#?m);
    V := Eexcess | RB;
    M := matrix(FF, apply(xF, p -> apply(V, m -> coefficient(m, p))));
    if o.Verbose then
        << "buildTemplateFromH: " << numRows M << "x" << numColumns M
           << " (excess=" << #Eexcess << ")" << endl << flush;
    (sh, M, V, Eexcess)
)

end
