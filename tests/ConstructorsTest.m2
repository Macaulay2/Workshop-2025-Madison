needsPackage "CpMackeyFunctors"

-- Test constructors
assert(isWellDefinedCpMackeyFunctor (BurnsideMackeyFunctor 17));
assert(isWellDefinedCpMackeyFunctor (FixedFreeMackeyFunctor 19));
assert(isWellDefinedCpMackeyFunctor (UnderlyingFreeMackeyFunctor 23));
assert(isWellDefinedCpMackeyFunctor (ComplexRepresentationMackeyFunctor 29));
assert(isWellDefinedCpMackeyFunctor (RealRepresentationMackeyFunctor 2));
assert(isWellDefinedCpMackeyFunctor (RealRepresentationMackeyFunctor 7));

assert(isWellDefinedCpMackeyFunctor makeZeroMackeyFunctor(3));