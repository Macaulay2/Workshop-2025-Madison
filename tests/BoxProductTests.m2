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