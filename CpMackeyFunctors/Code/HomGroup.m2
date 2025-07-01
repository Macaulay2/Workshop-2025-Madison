-- The Hom group Hom(M,N) is a subgroup of Hom(M(Cp/Cp),N(Cp/Cp)) ++ Hom(M(Cp/e),N(Cp/e))
Hom(CpMackeyFunctor, CpMackeyFunctor) := Module => opts -> (M,N) -> (
    homFixedFixed := Hom(getFixedModule M, getFixedModule N);
    homUnderlyingUnderlying := Hom(getUnderlyingModule M, getUnderlyingModule N);
    homFixedUnderlying := Hom(getFixedModule M, getUnderlyingModule N);
    homUnderlyingFixed := Hom(getUnderlyingModule M, getFixedModule N);
    lhs := directSum(homUnderlyingUnderlying, homFixedFixed);
    rhs := directSum(homFixedUnderlying, homUnderlyingFixed, homUnderlyingUnderlying);
    kernel map(rhs, lhs, matrix(
        {
            {-Hom(M.Res, getUnderlyingModule N), Hom(getFixedModule M, N.Res)},
            {-Hom(getUnderlyingModule M, N.Tr), Hom(M.Tr, getFixedModule N)},
            {Hom(M.Conj, getUnderlyingModule N) - Hom(getUnderlyingModule M, N.Conj), 0}
        }
    ))
)

Hom(CpMackeyFunctor, MackeyFunctorHomomorphism) := Matrix => opts -> (M,f) -> (

)

Hom(MackeyFunctorHomomorphism, CpMackeyFunctor) := Matrix => opts -> (f,N) -> (
    
)