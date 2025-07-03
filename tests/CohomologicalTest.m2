needsPackage "CpMackeyFunctors"

-- verifying some cohomological universal maps are well-defined
f := map(ZZ^3,ZZ^3,matrix({{0,1,0},{-1,-1,0},{0,0,1}}))
M := makeFixedPointMackeyFunctor(3,f);
assert isWellDefined(makeUniversalMapFixedCohomological(M, vector (matrix {{0},{0},{1}})))
assert isWellDefined(makeUniversalMapFixedCohomological(M, gens(getFixedModule M)))

B = makeUnderlyingFreeMackeyFunctor 17
assert isWellDefined(makeUniversalMapCohomological(B, gens(getUnderlyingModule B), gens(getFixedModule B)))