needsPackage "CpMackeyFunctors"

-- Check that the identity has trivial kernel and cokernel.
A = makeBurnsideMackeyFunctor(5);
K = (ker id_A, coker id_A);
for i to 1 do assert( all( {getTransfer K_i,    getRestriction K_i,
	                    getFixedModule K_i, getUnderlyingModule K_i,
	                    getConjugation K_i}, zero))

-- A non-trivial map between Mackey functors.
B = makeBurnsideMackeyFunctor 2;
U = makeUnderlyingFreeMackeyFunctor 2;
f = map(U, B, matrix {{2},{2}}, matrix {{2,4}});
K' = ker f;
assert( all({getUnderlyingModule K', getConjugation K', getRestriction K', getTransfer K'}, zero) )
assert( gens getFixedModule K'  == matrix {{2},{-1}} )
assert( rank getFixedModule K' == 1 )
C = coker f;
assert( presentation getFixedModule C      == matrix {{2,4}} )
assert( presentation getUnderlyingModule C == matrix {{2},{2}} )

-- Checking that direct sums are well-defined.
assert( isWellDefined(B ++ U) )
assert( isWellDefined(directSum {B,B,B,B}) )


-- Checking direct sum of homomorphisms

assert( isWellDefined(directSum({f,id_U,f})))