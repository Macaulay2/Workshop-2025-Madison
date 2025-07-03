needsPackage "CpMackeyFunctors"

-- Make sure we can make a random Cp Mackey functor with the numbers of generators and relations specified
assert isWellDefined makeRandomCpMackeyFunctor(5,{4,3,2,6})

-- Make sure we can build one with just a prime as input
assert isWellDefined makeRandomCpMackeyFunctor(2)

-- Make sure we can modify the generator bound
assert isWellDefined makeRandomCpMackeyFunctor(3,GenBound=>6)

assert isWellDefined (makeRandomCpMackeyFunctor(2,{1,2,3,4}) ++ makeRandomCpMackeyFunctor(2,{4,3,2,6}))
assert isWellDefined (makeRandomCpMackeyFunctor(5,GenBound=>6) ++ makeRealRepresentationMackeyFunctor(5))

assert isWellDefined (makeRandomMackeyFunctorHomomorphism(makeRealRepresentationMackeyFunctor(3), makeComplexRepresentationMackeyFunctor(3)))
assert isWellDefined (makeRandomMackeyFunctorHomomorphism(makeRandomCpMackeyFunctor(2,{1,2,3,4}), makeRandomCpMackeyFunctor(2,GenBound=>4) ))
assert isWellDefined (coker (makeRandomMackeyFunctorHomomorphism(makeRandomCpMackeyFunctor(2,{1,2,3,4}), makeRandomCpMackeyFunctor(2,GenBound=>3) )))

M = makeBurnsideMackeyFunctor(5);
N = makeUnderlyingFreeMackeyFunctor(5);
P = makeComplexRepresentationMackeyFunctor(5);
f = makeRandomMackeyFunctorHomomorphism(M,N);
g = makeRandomMackeyFunctorHomomorphism(N,P);
assert isWellDefined(g*f)
assert (source(g*f) == M)
assert (target(g*f) == P)