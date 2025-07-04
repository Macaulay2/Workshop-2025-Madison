restart
debug needsPackage "CpMackeyFunctors"

A = makeBurnsideMackeyFunctor 3
B = makeUnderlyingFreeMackeyFunctor 3
RU = makeComplexRepresentationMackeyFunctor 3
drawVerticalCpMackeyFunctor RU

rand1 = makeRandomCpMackeyFunctor 3
rand2 = makeRandomCpMackeyFunctor 3
makeRandomMackeyFunctorHomomorphism(rand1,rand2)

A**A
prune (A**A)
prune (A**A) == A

-- if (prune (A**A) === A) then (
--     print "so true bestie!!1!"
-- ) else (
--     print "it's over :("
-- )

RU**A
prune(RU**A)
RU

InternalHom(A,A)
prune InternalHom(A,A)

InternalHom(B,B)
prune InternalHom(B,B)

prune(Tor_0(A,A)) == prune(A**A)
isTrivialMackeyFunctor Tor_47(A,A)

-- K-theory example
q = 128 -- order of field 
p = 7 -- degree of extension
C = coker(matrix {{q^p-1}}) -- unit group of field
Ominus1 = makeFixedPointMackeyFunctor(p,inducedMap(C,C,matrix{{q}}))
Circ = makeZeroOnUnderlyingMackeyFunctor(p,coker(matrix{{p}}))

Ext^1(Circ,Ominus1)
isTrivialMackeyFunctor Ext^1(Circ,Ominus1)

-- equivariant cohomology example

Z = makeFixedPointMackeyFunctor(5,id_(ZZ^1))
M = makeZeroOnUnderlyingMackeyFunctor(5,coker(matrix{{5}}))
resolutionCohomological(M,5)

for i to 4 do (
    print prune(ExtCoh(i,M,Z))
)