needsPackage "CpMackeyFunctors"

-- not well defined if primes not the same
X := BurnsideMackeyFunctor 7;
Y := BurnsideMackeyFunctor 5;
f := map(getUnderlyingModule Y, getUnderlyingModule X, 0);
g := map(getFixedModule Y, getFixedModule X, 0);
assert not isWellDefined (new MackeyFunctorHomomorphism from {
    Domain => Y,
    Codomain => X,
    UnderlyingMap => f,
    FixedMap => g
});

-- verify the zero map to the zero Mackey functor is well-defined

A := BurnsideMackeyFunctor(7);
Z := makeZeroMackeyFunctor(7);

fixedLevelMap:=map(getFixedModule(Z),getFixedModule(A),0);
underlyingLevelMap:=map(getUnderlyingModule(Z),getUnderlyingModule(A),0);


F := makeMackeyFunctorHomomorphism(Z,A,underlyingLevelMap,fixedLevelMap);
assert isWellDefined(F)