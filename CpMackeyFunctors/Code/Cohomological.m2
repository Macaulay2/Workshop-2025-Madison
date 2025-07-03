------additional stuff for cohomological Mackey functors------
makeUniversalMapFixedCohomological = method()
-- Given a Mackey functor M and vector x in fixed module, produce map A -> M
makeUniversalMapFixedCohomological(CpMackeyFunctor,Vector) := MackeyFunctorHomomorphism => (M,x) -> makeUniversalMapFixedCohomological(M,matrix x)
-- Given a Mackey functor M and a matrix of n vectors in fixed module, produce map A^n -> M
makeUniversalMapFixedCohomological(CpMackeyFunctor,Matrix) := MackeyFunctorHomomorphism => (M,x) -> (
    if not isCohomological(M) then error "input is not a cohomological Mackey functor";
    n := numColumns x;
    if n == 0 then return map(M, makeZeroMackeyFunctor(M.PrimeOrder), 0);
    L := {for i to n-1 list (
        X := inducedMap(M.Fixed, , matrix x_i);
        -- TODO: should we error check element containment which is seemingly not implemented?
        p := M.PrimeOrder;
        Z := makeFixedPointMackeyFunctor(p,id_(ZZ^1));
        U := M.Res * X;
        F := X;
        map(M, Z, U, F)
	)};
    matrix L
)

-- Given a Cohomological Mackey functor M,
-- a matrix of n elements X in underlying, and
-- a matrix of m elements Y in fixed,
-- return the universal map B^n ++ Z^m -> M
makeUniversalMapCohomological = method()
makeUniversalMapCohomological(CpMackeyFunctor,Matrix,Matrix) := MackeyFunctorHomomorphism => (M,X,Y) -> (
    return makeUniversalMapUnderlying(M,X) | makeUniversalMapFixedCohomological(M,Y)
)
