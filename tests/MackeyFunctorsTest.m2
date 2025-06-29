needsPackage "CpMackeyFunctors"

-- Test if the Zero Mackey functor can be constructed and is well-defined
R:=matrix({});
C:=R;
T:=R;
M:=makeCpMackeyFunctor(7,R,T,C)
assert(isWellDefinedCpMackeyFunctor M)