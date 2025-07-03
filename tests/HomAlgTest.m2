needsPackage "CpMackeyFunctors"

F = (cokernel matrix {{84}}) ++ (cokernel matrix {{28}}) 
U = (cokernel matrix {{2}}) ++ module ZZ
r = map(U,F, matrix {{1,1}, {0,0}})
t = map(F,U, matrix {{42,42}, {0,14}})
c = map(U,U, matrix {{1, 0}, {0,-1}})
cursedMF := makeCpMackeyFunctor(2,r,t,c)

d = getResolution(cursedMF,3);

for i to (length d) - 2 do (
    assert(isTrivialMackeyFunctor(computeHomology(d_i,d_(i+1))))
)

