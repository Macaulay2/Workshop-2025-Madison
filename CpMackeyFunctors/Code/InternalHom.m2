InternalHom = method()

InternalHom (CpMackeyFunctor,CpMackeyFunctor) := CpMackeyFunctor => (M,N) -> (
    if not M.PrimeOrder === N.PrimeOrder then error "InternalHom: CpMackeyFunctors must have the same prime";
    p := M.PrimeOrder;
    A := makeBurnsideMackeyFunctor p;
    B := makeUnderlyingFreeMackeyFunctor p;
    underlying := Hom(B ** M, N);
    fixed := Hom(A ** M, N);
    tau := makeUniversalMapUnderlying(A, vector(getUnderlyingModule A,{1}));
    restriction := Hom(tau ** M, N);
    rho := makeUniversalMapFixed(B, vector(getFixedModule B,{1}));
    transfer := Hom(rho ** M, N);
    gen := {0,1} | for i to p-3 list 0;
    conj := Hom(makeUniversalMapUnderlying(B, vector(getUnderlyingModule B, gen)) ** M, N);
    makeCpMackeyFunctor(
        p,
        restriction,
        transfer,
        conj
    )
)

InternalHom (CpMackeyFunctor, MackeyFunctorHomomorphism) := MackeyFunctorHomomorphism => (M,N) -> (
    -- TODO
)

InternalHom (MackeyFunctorHomomorphism, CpMackeyFunctor) := MackeyFunctorHomomorphism => (M,N) -> (
    -- TODO
)