needsPackage "CpMackeyFunctors"

F = (cokernel matrix {{84}}) ++ (cokernel matrix {{28}}) 
U = (cokernel matrix {{2}}) ++ module ZZ
r = map(U,F, matrix {{1,1}, {0,0}})
t = map(F,U, matrix {{42,42}, {0,14}})
c = map(U,U, matrix {{1, 0}, {0,-1}})
cursedMF := makeCpMackeyFunctor(2,r,t,c)

d = getResolution(cursedMF,10);

-- check that the homology of a free resolution is zero
for i to 9 do (
    assert(isTrivialMackeyFunctor(computeHomology(d_i,d_(i+1))))
)

-- testing Tor
-- check that Tor_i(A,-) and Tor_i(B,-) are zero for i > 0
A = makeBurnsideMackeyFunctor 2;
B = makeBurnsideMackeyFunctor 2;
for i from 1 to 5 do (
    assert(isTrivialMackeyFunctor Tor_i(A,A));
    assert(isTrivialMackeyFunctor Tor_i(A,cursedMF));
    assert(isTrivialMackeyFunctor Tor_i(B,cursedMF))
)

-- texting Ext
-- check that Ext^0 is InternalHom
assert(prune Ext^0(A,A) == prune InternalHom(A,A))
assert(prune Ext^0(A,cursedMF) == prune InternalHom(A,cursedMF))

-- check that Ext^i(A,-) and Ext^i(B,-) are zero for i > 0
for i from 1 to 5 do (
    assert(isTrivialMackeyFunctor Ext^i(A,A));
    assert(isTrivialMackeyFunctor Ext^i(A,cursedMF));
    assert(isTrivialMackeyFunctor Ext^i(B,cursedMF))
)
