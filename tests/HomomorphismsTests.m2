needsPackage "CpMackeyFunctors"

-- not well defined if primes not the same
X := BurnsideMackeyFunctor 7;
Y := BurnsideMackeyFunctor 5;
f := map(getUnderlyingModule Y, getUnderlyingModule X, 0);
g := map(getFixedModule Y, getFixedModule X, 0);
assert not isWellDefinedCpMackeyFunctorHomomorphism (new MackeyFunctorHomomorphism from {
    Domain => Y,
    Codomain => X,
    UnderlyingMap => f,
    FixedMap => g
});

