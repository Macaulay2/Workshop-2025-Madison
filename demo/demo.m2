restart
debug needsPackage "CpMackeyFunctors"

A = makeBurnsideMackeyFunctor 5
B = makeUnderlyingFreeMackeyFunctor 5
RU = makeComplexRepresentationMackeyFunctor 5
drawVerticalCpMackeyFunctor RU

rand1 = makeRandomCpMackeyFunctor 5
rand2 = makeRandomCpMackeyFunctor 5
makeRandomMackeyFunctorHomomorphism(rand1,rand2)

A**A
prune (A**A)
prune (A**A) == A

RU**A
prune(RU**A)
RU

InternalHom(A,A)
prune InternalHom(A,A)

InternalHom(B,B)
prune InternalHom(B,B)

prune(Tor_0(A,A)) == prune(A**A)
isTrivialMackeyFunctor Tor_47(A,A)

q = 16 -- order of field 
p = 5 -- degree of extension
C = coker(matrix {{q^p-1}}) -- unit group of field
Ominus1 = makeFixedPointMackeyFunctor(p,inducedMap(C,C,matrix{{q}}))
Circ = makeZeroOnUnderlyingMackeyFunctor(p,coker(matrix{{p}}))

Ext^1(Circ,Ominus1)
isTrivialMackeyFunctor oo

Z = makeFixedPointMackeyFunctor(p,id_(ZZ^1))
B1 = makeZeroOnUnderlyingMackeyFunctor(p,coker(matrix{{p}}))
resolutionCohomological(B1,5)

for i to 4 do ( print i; print prune(ExtCoh(i,B1,Z)) )
