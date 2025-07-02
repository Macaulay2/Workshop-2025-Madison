needsPackage "CpMackeyFunctors"

assert isWellDefined makeRandomCpMackeyFunctor(5,{4,3,2,6})

assert isWellDefined makeRandomCpMackeyFunctor(2)
assert isWellDefined makeRandomCpMackeyFunctor(3,GenBound=>6)