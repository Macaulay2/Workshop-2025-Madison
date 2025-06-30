needsPackage "CpMackeyFunctors"

-- Test constructors
assert(isWellDefinedCpMackeyFunctor (BurnsideMackeyFunctor 17));
assert(isWellDefinedCpMackeyFunctor (FixedFreeMackeyFunctor 17));
assert(isWellDefinedCpMackeyFunctor (UnderlyingFreeMackeyFunctor 23));