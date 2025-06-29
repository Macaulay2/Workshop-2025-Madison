needsPackage "CpMackeyFunctors"

-- Test if the Zero Mackey functor can be constructed and is well-defined
R:=matrix({});
C:=R;
T:=R;
M:=makeCpMackeyFunctor(7,R,T,C)
assert(isWellDefinedCpMackeyFunctor M)

-- Test if constant F_2 is a C_2-Mackey functor
U:=cokernel(matrix({{2}}));
R:=inducedMap(U,U);
T:=inducedMap(U,U,matrix({{0}}));
C:=inducedMap(U,U);
M:=makeCpMackeyFunctor(2,R,T,C);
assert(isWellDefinedCpMackeyFunctor M);

-- Test if F_4 Galois is a C_2-Mackey functor
U:=cokernel(matrix({{2,0},{0,2}}));
C:=inducedMap(U,U,matrix({{1,1},{0,1}}));
F:=kernel(C - id_U);
R:=inducedMap(U,F,matrix({{1,0},{0,1}}));
T:=inducedMap(F,U,matrix({{0,1},{0,0}}));
M:=makeCpMackeyFunctor(2,R,T,C);
assert(isWellDefinedCpMackeyFunctor M);