needsPackage "CpMackeyFunctors"

-- Test constructors
assert(isWellDefinedCpMackeyFunctor (BurnsideMackeyFunctor 17));
assert(isWellDefinedCpMackeyFunctor (FixedFreeMackeyFunctor 19));
assert(isWellDefinedCpMackeyFunctor (UnderlyingFreeMackeyFunctor 23));
assert(isWellDefinedCpMackeyFunctor (ComplexRepresentationMackeyFunctor 29));