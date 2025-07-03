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

TEST /// -* representationRing for S_4 *-
(R,T) = representationRing(4, ReturnTable => true);
assert(numgens representationRing 4 == numColumns T);
assert(numgens representationRing 4 == 5)
///

TEST /// -* orbitPolytope of the  point {1/2,0,1} *-
P  = orbitPolytope(transpose matrix{{1/2,0,1}});
expected = transpose matrix{{1,1/2,0},{1/2,1,0},{1,0,1/2}, {0,1,1/2},{1/2,0,1},{0,1/2,1}};
assert(sort vertices P == sort expected)
///

TEST /// -* isSymmetric on a polytope symmetric to matrix{{0,1,0},{1,0,0},{0,0,1}} *-
P = orbitPolytope(transpose matrix{{1/2,0,1}});
assert(isSymmetric(P, matrix{{0,1,0},{1,0,0},{0,0,1}}))
///

TEST /// -* isSymmetric on a polytope with matrice that its not symmetric against*-
P = convexHull transpose matrix "1,0,0;0,2,0;0,0,3";
M = matrix{{1,1,0},{1,0,0},{0,1,1}}
assert(isSymmetric(P, M) == false)
///

TEST /// -* generateGroup on an element of order 4, checks that it returns in order *-
g = matrix{{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,0}}
Gs = generateGroup {g}
assert(length(Gs) == 4)
assert(Gs#0 == g^4)
assert(Gs#1 == g)
assert(Gs#2 == g^2)
///

TEST /// -* generateGroup genererates S4 when given generators for S4 *-
gList = {matrix{{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,0}}, matrix{{0,0,1,0},{1,0,0,0},{0,0,0,1},{0,1,0,0}}}
Gs = generateGroup gList
assert(#Gs == 24)
///

TEST /// -* fixedPolytope of the triangle formed by e_1,e_2,e_3 along e_1 <---> e_2 *-
P = orbitPolytope(transpose matrix{{1,0,0}});
fP = fixedPolytope(P, matrix{{0,1,0},{1,0,0},{0,0,1}});
assert(sort vertices fP == sort transpose matrix{{0,0,1},{1/2,1/2,0}})
///

TEST /// -* cycleTypeRepresentatives on 2 *-
assert(cycleTypeRepresentatives 2 == { matrix{{0,1},{1,0}},matrix{{1,0},{0,1}}})
///

TEST /// -* conjugacyClasses of {{0,1},{1,0}} *-
G = generateGroup {matrix{{0,1},{1,0}}};
ccG = conjugacyClasses G;
assert(ccG == {{matrix{{1,0},{0,1}}},{matrix{{0,1},{1,0}}}})
///
