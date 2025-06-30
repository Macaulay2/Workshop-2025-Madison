needsPackage "CpMackeyFunctors"

-- Test constructors
assert(isWellDefinedCpMackeyFunctor (BurnsideMackeyFunctor 2));
assert(isWellDefinedCpMackeyFunctor (UnderlyingFreeMackeyFunctor 2));