needsPackage "CpMackeyFunctors"

-- not well defined if primes not the same
X := makeBurnsideMackeyFunctor 7;
Y := makeBurnsideMackeyFunctor 5;
f := map(getUnderlyingModule Y, getUnderlyingModule X, 0);
g := map(getFixedModule Y, getFixedModule X, 0);
assert not isWellDefined (new MackeyFunctorHomomorphism from {
    Domain => Y,
    Codomain => X,
    UnderlyingMap => f,
    FixedMap => g
});

-- verify the zero map to the zero Mackey functor is well-defined

A := makeBurnsideMackeyFunctor(7);
Z := makeZeroMackeyFunctor(7);

fixedLevelMap:=map(getFixedModule(Z),getFixedModule(A),0);
underlyingLevelMap:=map(getUnderlyingModule(Z),getUnderlyingModule(A),0);
F := map(Z,A,underlyingLevelMap,fixedLevelMap);
assert isWellDefined(F)

-- verify that the linearization maps are well-defined
assert isWellDefined(complexLinearizationMap(5))
assert isWellDefined(realLinearizationMap(5))

-- verifying some universal maps are well-defined

RO := makeRealRepresentationMackeyFunctor 5;
assert isWellDefined(makeUniversalMapFixed(RO, vector (matrix {{4},{1},{-7}})))

B := makeUnderlyingFreeMackeyFunctor 3;
assert isWellDefined(makeUniversalMapUnderlying(B, vector (matrix {{0},{1},{0}})))

-- verify composition does what we want it to do ?
A :=makeBurnsideMackeyFunctor(11);
assert isWellDefined (id_(A) * id_(A))
assert (id_(A) * id_(A) === id_(A))

assert isWellDefined ( complexLinearizationMap(11) * id_(A));

-- Checking direct sum of homomorphisms
B = makeBurnsideMackeyFunctor 2;
U = makeUnderlyingFreeMackeyFunctor 2;
f = map(U, B, matrix {{2},{2}}, matrix {{2,4}});
assert( isWellDefined(directSum({f,id_U,f})))

-- Checking arithmetic of homomorphisms
h = id_B;
assert(h + h + h == 3 * h)
>>>>>>> c4a6adf (added equality, temporary prune)
