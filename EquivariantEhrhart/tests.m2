TEST /// -* isSymmetric hypersimplex(4,2) at 2134 *-
needsPackage "Permutations";
P = convexHull transpose matrix{
    {1,1,0,0},
    {1,0,1,0},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,0,1},
    {0,0,1,1}
    };
assert( isSymmetric(P, matrix permutation {2, 1, 3, 4}) );
///

TEST /// -* simple isEffective tests *-
    -- Triangle in R4
    P = convexHull transpose matrix "1,0,0,0;0,1,0,0;0,0,1,0" 
    g = matrix "0,1,0,0;1,0,0,0;0,0,1,0;0,0,0,1"
    assert(isEffective equivariantEhrhartSeries(P, {g}));

    -- Square in R4
    P = convexHull transpose matrix "1,1,0,0;0,1,1,0;0,0,1,1;1,0,0,1"
    g = matrix "0,0,0,1;1,0,0,0;0,1,0,0;0,0,1,0"
    assert(isEffective equivariantEhrhartSeries(P, {g^2}));

    -- hypersimplex
    hypersimplex = method()
    hypersimplex(ZZ, ZZ) := (n, k) -> (
        convexHull transpose matrix for s in subsets(n, k) list (
        for i from 0 to n-1 list if member(i, s) then 1 else 0
        )
        )
    P = hypersimplex(5,3);
    g = matrix randomPermutation ambDim P;
    assert(isEffective equivariantEhrhartSeries(P,{g}));
///