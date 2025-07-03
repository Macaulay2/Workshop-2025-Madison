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
assert isWellDefined(makeUniversalMapFixed(RO, gens(getFixedModule RO)))

B := makeUnderlyingFreeMackeyFunctor 3;
assert isWellDefined(makeUniversalMapUnderlying(B, vector (matrix {{0},{1},{0}})))
assert isWellDefined(makeUniversalMapUnderlying(B, gens(getUnderlyingModule B)))

assert isWellDefined(makeUniversalMap(B, gens(getUnderlyingModule B), gens(getFixedModule B)))

-- cursedMackeyFunctor
F = cokernel matrix {{84}}
U = (cokernel matrix {{2}}) ++ module ZZ
r = map(U,F, matrix {{1},{0}})
t = map(F,U, matrix {{42,42}})
c = map(U,U, matrix {{1, 0}, {0,-1}})
cursedMackeyFunctor := makeCpMackeyFunctor(2,r,t,c)
x = gens getFixedModule cursedMackeyFunctor
f1 = makeUniversalMapFixed(cursedMackeyFunctor,x)
x' = (gens getUnderlyingModule cursedMackeyFunctor)
f2 = makeUniversalMapUnderlying(cursedMackeyFunctor, x')
assert(matrix f1.FixedMap == matrix {{1,42}} and matrix f1.UnderlyingMap == matrix{{1},{0}})
assert(matrix f2.FixedMap == matrix {{42,42}} and matrix f2.UnderlyingMap == matrix{{1,1,0,0},{0,0,1,-1}})


-- verify composition does what we want it to do ?
A :=makeBurnsideMackeyFunctor(11);
assert isWellDefined (id_(A) * id_(A))
assert (id_(A) * id_(A) === id_(A))

assert isWellDefined ( complexLinearizationMap(11) * id_(A));


-- test isTrivialMackeyFunctor
assert not (isTrivialMackeyFunctor(makeBurnsideMackeyFunctor(7)));
assert isTrivialMackeyFunctor(makeZeroMackeyFunctor(3));

-- test isIsomorphism
assert not (isIsomorphism(realLinearizationMap(7)));
assert not (isIsomorphism(realLinearizationMap(11)));
assert isIsomorphism(realLinearizationMap(2));
assert isIsomorphism(realLinearizationMap(3));

-- Checking direct sum of homomorphisms
B = makeBurnsideMackeyFunctor 2;
U = makeUnderlyingFreeMackeyFunctor 2;
f = map(U, B, matrix {{2},{2}}, matrix {{2,4}});
assert( isWellDefined(directSum({f,id_U,f})))

-- Checking arithmetic of homomorphisms
h = id_B;
assert(h + h + h == 3 * h)



-- verifying some cohomological universal maps are well-defined
f := map(ZZ^3,ZZ^3,matrix({{0,1,0},{-1,-1,0},{0,0,1}}))
M := makeFixedPointMackeyFunctor(3,f);
assert isWellDefined(makeUniversalMapFixedCohomological(M, vector (matrix {{0},{0},{1}})))
assert isWellDefined(makeUniversalMapFixedCohomological(M, gens(getFixedModule M)))

assert isWellDefined(makeUniversalMapCohomolgical(B, gens(getUnderlyingModule B), gens(getFixedModule B)))