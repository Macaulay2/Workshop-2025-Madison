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

RU**A
prune(RU**A)
RU

InternalHom(A,A)
prune InternalHom(A,A)

Tor_1(rand1,rand2)
prune oo
Ext^1(rand1,rand2)
prune oo

Tor_0(A,A)
prune(Tor_0(A,A)) == prune(A**A)

q = 16 -- order of field 
p = 5 -- degree of extension
C = coker(matrix {{q^p-1}}) -- unit group of field
K1 = makeFixedPointMackeyFunctor(p,inducedMap(C,C,matrix{{q}}))
C = makeZeroOnUnderlyingMackeyFunctor(p,coker(matrix{{p}}))

Ext^1(K1,C)
prune oo

Z = makeFixedPointMackeyFunctor(p,id_(ZZ^1))
B1 = makeZeroOnUnderlyingMackeyFunctor(p,coker(matrix{{p}}))
resolutionCohomological(B1,5)

for i to 4 do ( print i; print prune(ExtCoh(i,B1,Z)) )
