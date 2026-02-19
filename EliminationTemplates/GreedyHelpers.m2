-- Greedy strategy helper methods used by EliminationTemplates.m2.

-- Input: a matrix H over a polynomial ring
-- Output: (vH, W, Zlist)

monomialVectorAndW = method()

monomialVectorAndWData = method()

monomialVectorAndWData Matrix := (H) -> (
    R := ring H;

    -- 1) collect monomials from the coefficient ring of H's entries
    -- (e.g. if ring H = R[t], collect monomials in R, not in t)
    SH := ring H;
    baseR := coefficientRing SH;
    monsList := flatten apply(flatten entries H, f -> (
        coeffsF := flatten entries(last coefficients f);
        flatten apply(coeffsF, c -> flatten entries monomials sub(c, baseR))
    ));

    mons := toList set monsList;

    -- turn into a column vector vH (as a matrix with 1 column)
    vH := matrix apply(#mons, i -> {mons#i});  -- (#mons) x 1

    -- 2) for each column h_k, find Z_k such that h_k = Z_k * vH
    ncols := numColumns H;
    nrows := numRows H;

    Zlist := for k from 0 to ncols-1 list (
        hk := submatrix(H, (0..nrows-1), {k});  -- nrows x 1 column

        -- Build Zk as nrows x (#mons) matrix:
        -- entry (i,j) is coefficient of monomial mons#j in hk_i
        Zk := matrix apply(nrows, i -> (
            fi := hk_(i,0);
            (tMonsMat, coeffsMat) := coefficients fi;
            tMons := flatten entries tMonsMat;
            coeffs := flatten entries coeffsMat;
            apply(#mons, j -> (
                m := mons#j;
                parts := apply(#tMons, q -> sub(coefficient(m, sub(coeffs#q, baseR)), SH) * tMons#q);
                if #parts == 0 then 0_SH else fold(parts, (a,b) -> a+b)
            ))
        ));

        Zk
    );

    -- 3) stack all Z_k horizontally to get W
    W := fold((A,B) -> A | B, first Zlist, drop(Zlist,1));
    columnInfo := flatten for k from 0 to ncols-1 list (
        apply(#mons, j -> new HashTable from {
            "generatorIndex" => k,
            "monomial" => mons#j
        })
    );

    new HashTable from {
        "vH" => vH,
        "W" => W,
        "Zlist" => Zlist,
        "mons" => mons,
        "columnInfo" => columnInfo
    }
);

monomialVectorAndW Matrix := (H) -> (monomialVectorAndWData H)#"W";

-- Usage:
-- (vH, W, Zlist) = monomialVectorAndW(H);

copyAssignments = (A) -> (
    B := new MutableHashTable from {};
    scan(keys A, k -> B#k = A#k);
    B
);

assignmentIdeal = (A, SH) -> (
    if #keys A == 0 then ideal(0_SH)
    else ideal apply(keys A, v -> v - A#v)
);

reduceWithAssignments = (f, A, SH) -> (
    if #keys A == 0 then f else f % assignmentIdeal(A, SH)
);

isConstantInBaseRing = (c, baseR) -> (
    mons := flatten entries monomials c;
    all(mons, m -> m == 1_baseR)
);

sumInRing = (L, SH) -> if #L == 0 then 0_SH else fold(L, (a,b) -> a+b);

countZeroColumns = (W, A, SH) -> (
    if numColumns W == 0 then 0
    else #select(0..numColumns W - 1, k -> (
        all(0..numRows W - 1, i -> (
            reduceWithAssignments(W_(i,k), A, SH) == 0_SH
        ))
    ))
);

enforceZeroForColumns = (W, colsToZero, thetaVars, thetaToZeroMap, A, SH, baseR) -> (
    scan(colsToZero, k -> (
        scan(0..numRows W - 1, i -> (
            eq0 := reduceWithAssignments(W_(i,k), A, SH);
            if eq0 =!= 0_SH then (
                coeffs := apply(thetaVars, t -> reduceWithAssignments(coefficient(t, eq0), A, SH));
                linPart := sumInRing(apply(#thetaVars, q -> coeffs#q * thetaVars#q), SH);
                constPart := reduceWithAssignments(eq0 - linPart, A, SH);
                -- We only support affine-linear constraints in theta variables.
                if (constPart - thetaToZeroMap constPart) =!= 0_SH then return false;
                active := select(0..#thetaVars - 1, q -> coeffs#q =!= 0_SH);
                if #active == 0 then return false;
                if #active > 1 then return false;
                q := first active;
                t := thetaVars#q;
                aCoeff := sub(coeffs#q, baseR);
                bCoeff := sub(constPart, baseR);
                if not isConstantInBaseRing(aCoeff, baseR) then return false;
                if not isConstantInBaseRing(bCoeff, baseR) then return false;
                aScalar := coefficient(1_baseR, aCoeff);
                bScalar := coefficient(1_baseR, bCoeff);
                if aScalar == 0_(coefficientRing baseR) then return false;
                val := sub(-bScalar / aScalar, SH);
                if A#?t then (
                    if reduceWithAssignments(A#t - val, A, SH) =!= 0_SH then return false
                ) else A#t = val;
            );
        ));
    ));
    true
);

rowWiseGreedyAssignments = (W, thetaVars, thetaToZeroMap, SH, baseR) -> (
    A := new MutableHashTable from {};
    improved := true;
    while improved do (
        improved = false;
        baseZeroCount := countZeroColumns(W, A, SH);
        bestScore := 0;
        bestA := null;
        if numColumns W > 0 then (
            for k from 0 to numColumns W - 1 do (
                trial := copyAssignments A;
                if enforceZeroForColumns(W, {k}, thetaVars, thetaToZeroMap, trial, SH, baseR) then (
                    score := countZeroColumns(W, trial, SH) - baseZeroCount;
                    if score > bestScore then (
                        bestScore = score;
                        bestA = trial;
                    );
                );
            )
        );
        if bestScore > 0 then (
            A = bestA;
            improved = true;
        );
    );
    A
);

computeExcessiveMonomials = (a, B, J, columnInfo, baseR) -> (
    gensJ := toList J_*;
    allMons := set flatten apply(columnInfo, info -> (
        genIdx := info#"generatorIndex";
        m := info#"monomial";
        flatten entries monomials(m * sub(gensJ#genIdx, baseR))
    ));
    monsB := set flatten entries(lift(B, baseR));
    monsR := set flatten entries(a * lift(B, baseR)) - set flatten entries(lift(B, baseR));
    toList(allMons - union(monsR, monsB))
);

columnsForExcessiveMonomial = (e, columnInfo, J, baseR) -> (
    gensJ := toList J_*;
    select(0..#columnInfo - 1, c -> (
        info := columnInfo#c;
        genIdx := info#"generatorIndex";
        m := info#"monomial";
        coefficient(e, m * sub(gensJ#genIdx, baseR)) =!= 0_baseR
    ))
);

columnWiseGreedyAssignments = (W, excessiveMons, columnInfo, J, thetaVars, thetaToZeroMap, SH, baseR) -> (
    A := new MutableHashTable from {};
    improved := true;
    while improved do (
        improved = false;
        baseZeroCount := countZeroColumns(W, A, SH);
        bestScore := 0;
        bestA := null;
        scan(excessiveMons, e -> (
            colsE := columnsForExcessiveMonomial(e, columnInfo, J, baseR);
            if #colsE > 0 then (
                trial := copyAssignments A;
                if enforceZeroForColumns(W, colsE, thetaVars, thetaToZeroMap, trial, SH, baseR) then (
                    score := countZeroColumns(W, trial, SH) - baseZeroCount;
                    if score > bestScore then (
                        bestScore = score;
                        bestA = trial;
                    );
                );
            );
        ));
        if bestScore > 0 then (
            A = bestA;
            improved = true;
        );
    );
    A
);

instantiateHWithAssignments = (H, A, SH, baseR) -> (
    Hred := matrix apply(numRows H, i -> apply(numColumns H, j -> reduceWithAssignments(H_(i,j), A, SH)));
    sub(Hred, baseR)
);
