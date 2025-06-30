Hom(CpMackeyFunctor, CpMackeyFunctor) := Module => opts -> (N,M) -> (
    homFixedFixed := Hom(getFixedModule N, getFixedModule M);
    homUnderlyingUnderlying := Hom(getUnderlyingModule N, getUnderlyingModule M);
    homFixedUnderlying := Hom(getFixedModule N, getUnderlyingModule M);
    homUnderlyingFixed := Hom(getUnderlyingModule N, getFixedModule M);
    lhs := directSum(homFixedFixed, homUnderlyingUnderlying);
    rhs := directSum(homFixedUnderlying, homUnderlyingFixed, homUnderlyingUnderlying);
    kernel map(rhs, lhs, matrix(
        {
            {Hom(getFixedModule M, N.Res), -Hom(M.Res, getUnderlyingModule N)},
            {Hom(M.Tr, getFixedModule N), -Hom(getUnderlyingModule M, N.Tr)},
            {0, Hom(M.Conj, getUnderlyingModule N) - Hom(getUnderlyingModule M, N.Conj)}
        }
    ))
);