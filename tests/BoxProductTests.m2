needsPackage "CpMackeyFunctors"
needsPackage "Isomorphism"

-- Verify that the dimensions match up when box product with Burnside
A = makeBurnsideMackeyFunctor(7)
R = makeComplexRepresentationMackeyFunctor(7)
assert(rank getUnderlyingModule R == rank getUnderlyingModule (A**R))
assert(rank getUnderlyingModule R == rank getUnderlyingModule (R**A))
assert(rank getFixedModule R == rank getFixedModule (A**R))
assert(rank getFixedModule R == rank getFixedModule (R**A))

-- Verify that fixed modules of A and A ** R are isomorphic
isIsomorphic(getFixedModule(R), getFixedModule(A ** R))

-- Verifying that induced maps on box products are well-defined
L = complexLinearizationMap(7)
assert(isWellDefined boxProduct(A,L))
assert(isWellDefined boxProduct(L,A))

-- A ** A == A. We can see this with prune!
A = makeBurnsideMackeyFunctor 5
B = A ** A
assert(prune B === A)
f = (prune B).cache.pruningMap;
assert isWellDefined f
assert isWellDefined inverse f

-- Test induced map on box products
for i to 5 do (
    p = (random {2,3,5,7,11,13})#0;
    n = (random toList (0..20))#0;
    M = makeRandomCpMackeyFunctor(p);
    A = makeBurnsideMackeyFunctor(p);
    timesn = map(A,A,n);
    assert (M**timesn == map(M**A,M**A,n));
    assert (timesn**M == map(A**M,A**M,n));
);