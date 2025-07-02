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

InternalHom (CpMackeyFunctor, MackeyFunctorHomomorphism) := MackeyFunctorHomomorphism => (M,f) -> (
    if not M.PrimeOrder === (source f).PrimeOrder then error "InternalHom: Mackey functors are over different primes";
    p := M.PrimeOrder;
    A := makeBurnsideMackeyFunctor p;
    B := makeUnderlyingFreeMackeyFunctor p;
    src := Hom(M, source f);
    tgt := Hom(M, target f);
    underlyingMap := Hom(B ** M, f.UnderlyingMap);
    fixedMap := Hom(A ** M, f.FixedMap);
    map(tgt, src, underlyingMap, fixedMap)
)

InternalHom (MackeyFunctorHomomorphism, CpMackeyFunctor) := MackeyFunctorHomomorphism => (f,N) -> (
    -- TODO
)