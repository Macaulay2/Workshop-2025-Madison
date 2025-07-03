needsPackage "CpMackeyFunctors"

F = (cokernel matrix {{84}}) ++ (cokernel matrix {{28}}) 
U = (cokernel matrix {{2}}) ++ module ZZ
r = map(U,F, matrix {{1,1}, {0,0}})
t = map(F,U, matrix {{42,42}, {0,14}})
c = map(U,U, matrix {{1, 0}, {0,-1}})
cursedMF := makeCpMackeyFunctor(2,r,t,c)

-- verify that constructed surjections are indeed surjective
assert(isTrivialMackeyFunctor (cokernel (makeFreeModuleSurjection makeRealRepresentationMackeyFunctor 17)))
assert(isTrivialMackeyFunctor (cokernel (makeFreeModuleSurjection cursedMF)))

-- verify free resolutions are in fact complexes
d = getResolution(cursedMF,3)
for i to (length d) - 2 do (
    comp = d#i * d#(i+1);
    assert(getUnderlyingMap(comp) == 0 and getFixedMap(comp) == 0)
)




