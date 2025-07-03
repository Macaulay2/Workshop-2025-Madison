needsPackage "CpMackeyFunctors"

F = (cokernel matrix {{4}}) 
U = (cokernel matrix {{2}})
r = map(U,F, matrix {{1}})
t = map(F,U, matrix {{2}})
c = map(U,U, matrix {{1}})
MF := makeCpMackeyFunctor(2,r,t,c)


-- verifying some cohomological universal maps are well-defined
f := map(ZZ^3,ZZ^3,matrix({{0,1,0},{-1,-1,0},{0,0,1}}))
M := makeFixedPointMackeyFunctor(3,f);
assert isWellDefined(makeUniversalMapFixedCohomological(M, vector (matrix {{0},{0},{1}})))
assert isWellDefined(makeUniversalMapFixedCohomological(M, gens(getFixedModule M)))

B = makeUnderlyingFreeMackeyFunctor 17
assert isWellDefined(makeUniversalMapCohomological(B, gens(getUnderlyingModule B), gens(getFixedModule B)))

-- verify that constructed surjections are indeed surjective
assert(isTrivialMackeyFunctor (cokernel (makeFreeModuleSurjection MF)))

-- verify free resolutions are in fact complexes
d = resolutionCohomological(MF,3)
for i to (length d) - 2 do (
    comp = d#i * d#(i+1);
    assert(getUnderlyingMap(comp) == 0 and getFixedMap(comp) == 0)
)
