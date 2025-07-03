needsPackage "CpMackeyFunctors"

F = (cokernel matrix {{84}}) ++ (cokernel matrix {{28}}) 
U = (cokernel matrix {{2}}) ++ module ZZ
r = map(U,F, matrix {{1,1}, {0,0}})
t = map(F,U, matrix {{42,42}, {0,14}})
c = map(U,U, matrix {{1, 0}, {0,-1}})
cursedMF := makeCpMackeyFunctor(2,r,t,c)

d = getResolution(cursedMF,3);

-- check that the homology of a free resolution is zero
for i to (length d) - 2 do (
    assert(isTrivialMackeyFunctor(computeHomology(d_i,d_(i+1))))
)

-- check that Tor_i(A,-) and Tor_i(B,-) is zero for i > 0
A = makeBurnsideMackeyFunctor 2;
B = makeBurnsideMackeyFunctor 2;
for i from 1 to 3 do (
    assert(isTrivialMackeyFunctor Tor_i(A,A));
    assert(isTrivialMackeyFunctor Tor_i(A,cursedMF));
    assert(isTrivialMackeyFunctor Tor_i(B,cursedMF))
)

