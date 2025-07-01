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