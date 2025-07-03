-- The Hom group Hom(M,N) is a subgroup of Hom(M(Cp/e),N(Cp/e)) ++ Hom(M(Cp/Cp),N(Cp/Cp))
Hom(CpMackeyFunctor, CpMackeyFunctor) := Module => opts -> (M,N) -> (
    if not M.PrimeOrder === N.PrimeOrder then error "-- Mackey functors must have the same prime";
    homFixedFixed := Hom(getFixedModule M, getFixedModule N);
    homUnderlyingUnderlying := Hom(getUnderlyingModule M, getUnderlyingModule N);
    homFixedUnderlying := Hom(getFixedModule M, getUnderlyingModule N);
    homUnderlyingFixed := Hom(getUnderlyingModule M, getFixedModule N);
    lhs := directSum(homUnderlyingUnderlying, homFixedFixed);
    rhs := directSum(homFixedUnderlying, homUnderlyingFixed, homUnderlyingUnderlying);
    result := kernel map(rhs, lhs, matrix(
        {
            {-Hom(M.Res, getUnderlyingModule N), Hom(getFixedModule M, N.Res)},
            {-Hom(getUnderlyingModule M, N.Tr), Hom(M.Tr, getFixedModule N)},
            {Hom(M.Conj, getUnderlyingModule N) - Hom(getUnderlyingModule M, N.Conj), 0}
        }
    ));
    result.cache.homomorphism = f -> (
        map(
            N, M,
            homomorphism (inducedMap(homUnderlyingUnderlying, result, lhs^[0]) * f),
            homomorphism (inducedMap(homFixedFixed, result, lhs^[1]) * f)
        )
    );
    result.cache.formation = FunctionApplication { Hom, (M, N) };
    result
)

-- Covariant Hom
Hom(CpMackeyFunctor, MackeyFunctorHomomorphism) := Matrix => opts -> (M,f) -> (
    if not M.PrimeOrder === (source f).PrimeOrder then error "-- Mackey functors must have the same prime";
    T := getFixedModule(M);
    B := getUnderlyingModule(M);
    fT := f.FixedMap;
    fB := f.UnderlyingMap;
    phi := Hom(T, fT);
    psi := Hom(B, fB);
    inducedMap(Hom(M,target(f)), Hom(M, source(f)), psi ++ phi) 
)

-- Contravariant Hom
Hom(MackeyFunctorHomomorphism, CpMackeyFunctor) := Matrix => opts -> (f,N) -> (
    if not (source f).PrimeOrder === N.PrimeOrder then error "-- Mackey functors must have the same prime";
    T := getFixedModule(N);
    B := getUnderlyingModule(N);
    fT := f.FixedMap;
    fB := f.UnderlyingMap;
    phi := Hom(fT, T);
    psi := Hom(fB, B);
    inducedMap(Hom(source(f),N), Hom(target(f),N), psi ++ phi) 
)