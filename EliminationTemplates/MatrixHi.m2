-- MatrixHi.m2: port of the gap-polynomial / matrixHi construction from
-- Martyushev et al., "Optimizing Elimination Templates by Greedy Parameter
-- Search", CVPR 2022 (arXiv:2203.14901), Maple reference at
-- https://github.com/martyushev/eliminationTemplates (_greedyAG/greedyAG.mw).
--
-- Unlike the default/Larsson/Greedy pipeline (which extends the ideal by a
-- graph generator `s - a` to read off the action matrix), this construction
-- returns the template built directly from the gap polynomials
--   G_i = a * B_i - NF(a * B_i, J)         (for a * B_i not in B)
-- factored as G_i = sum_j H[i,j] * F_j.
--
-- The resulting template has fewer rows (no graph-ideal inflation by |B|)
-- and matches Martyushev's paper sizes on standard benchmarks.

-- ----------------------------------------------------------------------------
-- buildGapPolys: for action element a, quotient basis B, ideal J, returns
--   (gapPolys, residualMons, Blist)
-- where residualMons = { a*m : m in B, a*m not in B } and
--       gapPolys[i]   = residualMons[i] - (residualMons[i] % J).
-- ----------------------------------------------------------------------------
buildGapPolys = method()
buildGapPolys(RingElement, Matrix, Ideal) := (a, B, J) -> (
    R := ring J;
    a2 := sub(a, R);
    Blist := flatten entries lift(B, R);
    Bset := set Blist;
    residualMons := select(apply(Blist, m -> a2 * m), m -> not Bset#?m);
    gapPolys := apply(residualMons, r -> r - (r % J));
    (gapPolys, residualMons, Blist)
)

-- ----------------------------------------------------------------------------
-- matrixHiRow: for a single gap polynomial Gi, find polynomial cofactors
-- uH[j] (one per generator in F) with Gi = sum_j uH[j] * F[j], where each
-- uH[j] is drawn from basis(0, dHmin[j] + d, R) for the smallest d making
-- the linear system feasible (matching Maple's undetPoly in matrixHi).
--
-- Returns (uHrow, entMons, alphaExprs, freeCols, paramMap, nParamsPerGen).
-- ----------------------------------------------------------------------------
matrixHiRow = method(Options => {Verbose => false})
matrixHiRow(RingElement, List) := o -> (Gi, F) -> (
    R := ring Gi;
    FF := coefficientRing R;
    nF := #F;
    dF := apply(F, f -> first degree f);
    dG := first degree Gi;
    dHmin := apply(nF, j -> max(0, dG - dF#j));

    d := 0;
    dMinMet := 0;
    result := null;
    while d - dMinMet <= 0 do (
        entMons := apply(nF, j -> flatten entries basis(0, dHmin#j + d, R));
        nParamsPerGen := apply(entMons, em -> #em);
        nParamsTotal := sum nParamsPerGen;

        -- Precompute ek * F[j] once per (j,k) and hash by monomial.
        prodCoeffs := apply(nF, j -> apply(entMons#j, ek -> (
            p := ek * F#j;
            (mons, cfs) := coefficients p;
            hashTable apply(numColumns mons, c -> (mons_(0,c), cfs_(c,0)))
        )));

        targetSet := new MutableHashTable;
        scan(prodCoeffs, pj -> scan(pj, ht -> scan(keys ht, m -> targetSet#m = true)));
        scan(flatten entries monomials Gi, m -> targetSet#m = true);
        targetMons := keys targetSet;

        GiCoeffs := (
            (gm, gc) := coefficients Gi;
            hashTable apply(numColumns gm, c -> (gm_(0,c), gc_(c,0)))
        );
        zeroFF := 0_FF;

        Amat := matrix(FF, apply(targetMons, m ->
            flatten apply(nF, j -> apply(#(entMons#j), k -> (
                ht := prodCoeffs#j#k;
                if ht#?m then sub(ht#m, FF) else zeroFF
            )))
        ));
        bvec := matrix(FF, apply(targetMons, m ->
            {if GiCoeffs#?m then sub(GiCoeffs#m, FF) else zeroFF}));

        rkA := rank Amat;
        rkAug := rank(Amat | bvec);

        if rkAug > rkA then (
            dMinMet = d + 1;
            d = d + 1;
            if o.Verbose then << "  d=" << d-1 << " infeasible, retry d=" << d << endl;
        ) else (
            rref := reducedRowEchelonForm(Amat | bvec);
            pivotCols := {};
            curRow := 0;
            for col from 0 to nParamsTotal - 1 do (
                if curRow >= numRows rref then break;
                if rref_(curRow, col) == 1_FF then (
                    pivotCols = append(pivotCols, (col, curRow));
                    curRow = curRow + 1;
                );
            );
            pivotSet := set apply(pivotCols, p -> p#0);
            freeCols := select(0..nParamsTotal-1, c -> not pivotSet#?c);

            alphaExprs := apply(nParamsTotal, p -> (
                if pivotSet#?p then (
                    pivRow := (select(pivotCols, pc -> pc#0 == p))#0#1;
                    constant := rref_(pivRow, nParamsTotal);
                    deps := apply(freeCols, fc -> (fc, -rref_(pivRow, fc)));
                    (constant, deps)
                ) else (
                    (0_FF, {(p, 1_FF)})
                )
            ));

            initFreeVals := apply(freeCols, fc -> 0_FF);
            evalAlpha := (p, freeVals) -> (
                (constant, deps) := alphaExprs#p;
                val := constant;
                scan(#deps, k -> val = val + deps#k#1 * freeVals#(position(freeCols, c -> c == deps#k#0)));
                val
            );
            uHrow := apply(nF, j -> (
                pStart := sum apply(j, jj -> nParamsPerGen#jj);
                sum apply(#(entMons#j), k -> evalAlpha(pStart + k, initFreeVals) * entMons#j#k)
            ));

            paramMap := flatten apply(nF, j -> apply(#(entMons#j), k -> (j, k, entMons#j#k)));
            result = (uHrow, entMons, alphaExprs, freeCols, paramMap, nParamsPerGen);
            break;
        );
    );

    if result === null then error "matrixHiRow: failed to find solution";
    result
)

-- ----------------------------------------------------------------------------
-- matrixHi: apply matrixHiRow to each gap polynomial. Returns (Hmat, perRow)
-- where Hmat is nG x nF with entries in R (particular solution).
-- ----------------------------------------------------------------------------
matrixHi = method(Options => {Verbose => false})
matrixHi(List, List) := o -> (F, G) -> (
    R := ring first F;
    nG := #G;
    nF := #F;
    if o.Verbose then << "matrixHi: nG=" << nG << " nF=" << nF << endl << flush;
    perRowInfo := apply(nG, i -> (
        t0 := cpuTime();
        if o.Verbose then (<< "  row " << i << " (degG=" << first degree G#i << "): " << flush);
        rowResult := matrixHiRow(G#i, F, Verbose => o.Verbose);
        if o.Verbose then (<< #(rowResult#3) << " free, " << round(2, cpuTime()-t0) << "s" << endl << flush);
        rowResult
    ));
    Hmat := matrix(R, apply(perRowInfo, info -> info#0));
    (Hmat, perRowInfo)
)

-- ----------------------------------------------------------------------------
-- constructTemplate: from H, F, and RB = residualMons | Blist, build the
-- elimination template matrix M over the base field. For each generator
-- F[j], the shift monomials come from column j of H. The template columns
-- are [ excessMons | RB ].
-- ----------------------------------------------------------------------------
buildMatrixHiTemplate = method(Options => {Verbose => false})
buildMatrixHiTemplate(Matrix, List, List) := o -> (Hmat, F, RB) -> (
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
        << "constructTemplate: " << numRows M << "x" << numColumns M
           << " (excess=" << #Eexcess << ")" << endl << flush;
    (sh, M, V, Eexcess)
)
