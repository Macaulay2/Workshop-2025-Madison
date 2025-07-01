needsPackage "CpMackeyFunctors"

A = makeBurnsideMackeyFunctor(5);
K = ker id_A;
assert(getTransfer(K) == 0);
assert(getRestriction(K) == 0);
assert(getFixedModule(K)==0);
assert(getUnderlyingModule(K)==0);
assert(getConjugation(K)==0);


B = makeBurnsideMackeyFunctor(2);
U = makeUnderlyingFreeMackeyFunctor(2);
f = map(U,B,matrix {{2},{2}},matrix {{2,4}});
K = ker f;
assert(getUnderlyingModule(K)==0);
assert(gens(getFixedModule(K))== matrix {{2},{-1}});
assert(getConjugation(K)==0);
assert(getRestriction(K) == 0);
assert(getTransfer(K) == 0);

assert(isWellDefined(B ++ U))