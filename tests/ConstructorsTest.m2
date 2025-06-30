needsPackage "CpMackeyFunctors"

-- Test constructors
assert(isWellDefined (BurnsideMackeyFunctor 17));
assert(isWellDefined (FixedFreeMackeyFunctor 19));
assert(isWellDefined (UnderlyingFreeMackeyFunctor 23));
assert(isWellDefined (ComplexRepresentationMackeyFunctor 29));
assert(isWellDefined (RealRepresentationMackeyFunctor 2));
assert(isWellDefined (RealRepresentationMackeyFunctor 7));

assert(isWellDefined makeZeroMackeyFunctor(3));

assert(isWellDefined makeFixedPointMackeyFunctor(41,matrix{}));
assert(isWellDefined makeFixedPointMackeyFunctor(2,matrix{{0,1,0},{1,0,0},{0,0,1}}));