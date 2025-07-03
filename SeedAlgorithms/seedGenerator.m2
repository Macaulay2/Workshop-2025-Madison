--  File Created: 06/25/25
--  CREDITS  --
--  Advisor: Francesca Gandini
--  Theory: Sumner Strom & Marcus Cassel
--  Code: Gordon Novak, Marcus Cassel
--      * Last Editied 06/26/25
--  Documentation: Sasha Arasha & Gordon Novak
--      * Last Editied 06/26/25
--  Review: Sasha Arasha
--      * Last Reviewed 06/25/25

-- //////////////// --
-- //////////////// --
-- CODE STARTS HERE --
-- //////////////// --
-- //////////////// --


squareSubmatrixGenerator = (M) -> (

    n := numRows M;
    m := numColumns M;

    if (m < n) then (return {});

    subMatrixRep = subsets(toList(0..m-1), n);

    matrixList := {};
    for i in subMatrixRep do (
        finalMatrix := matrix{};
        
        for j to numRows M - 1 do (
            rowMatrix = matrix{};
            for k in i do (
                if (numRows rowMatrix == 0) then (
                    rowMatrix = matrix(M_k)_j;
                ) else (
                    rowMatrix = rowMatrix | matrix(M_k)_j;
                );
            );
            if (j == 0) then (
                finalMatrix = rowMatrix;
            ) else (
                finalMatrix = finalMatrix || rowMatrix;
            );
        );
        matrixList = matrixList | {finalMatrix};
    );

    return matrixList;
    
)


-- METHOD_NAME: elementaryAbelianPGroupSeedGenerator
-- METHOD PURPOSE: generates the seed for an abelian p group
-- INPUTS: (W, p, R)
--      * W: Weight Matrix
--      * p: p-group (cyclic subgroup of order p or cartesian product of cyclic subgroups of order p)
--      * R: Ring of the group action
-- OUTPUT: Polynomial 
--      * Outputs the invariant seed utilized to generate all possible invariants
elementaryAbelianPGroupSeedGenerator = (W, p, R) -> (

    -- First we grab the number of columns and rows
    n = numColumns W;
    m = numRows W;

    -- Then, we store them in the h variable. 
    -- h = factorial(n)/(factorial(m) * factorial(n-m));

    W2 = matrix ({gens R} | toList entries W); -- appends a symbolic row at the top of the weight matrix

    moddedCoefficientList := {};

    for miniMatrix in squareSubmatrixGenerator(W2) do (

        
        --takes the determinant of the new matrix which returns symbolic with coefficients that are the exponent vector of the seed polynomial
        detMiniMatrix = determinant miniMatrix; 

        -- this extracts the coefficients and puts them into a exponent vector

        
        count = 0;
        miniCoef = coefficients detMiniMatrix;
        coef := {};
        for r to #(gens R) - 1 do (
            existingCoef := false;
            for x in (flatten entries miniCoef_0) do (
                if (x == (gens R)_r) then (existingCoef = true);
            );
            
            if (existingCoef) then (
                
                coef = coef | {(flatten entries miniCoef_1)_count};
                count = count + 1;
            )
            else (
                coef = coef | {0};
            );
        );

        -- reduces mod p for our invariant seed
        modCoef = for c in coef list (
            cc = value toString c;
            if (cc < 0) then (
                cc + p
            )
            else (
                cc
            )
        );

        moddedCoefficientList = moddedCoefficientList | {modCoef};
    );

    for i from 0 to ((#moddedCoefficientList) - 1) when (i < (#moddedCoefficientList)) do (
        for k from i+1 to ((#moddedCoefficientList) - 1) when (k < (#moddedCoefficientList)) do (
            if (moddedCoefficientList_i == moddedCoefficientList_k) then (
                moddedCoefficientList = drop(moddedCoefficientList, {k,k});
            );
        );
    );

    return moddedCoefficientList;
)
