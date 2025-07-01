

directSumCpMackeyFunctors = method()
directSumCpMackeyFunctors(CpMackeyFunctor,CpMackeyFunctor) := CpMackeyFunctor => (M,N) -> (
    if not getPrimeOrder(M) == getPrimeOrder(N) then error("-- Primes not the same (incompatable)");
    T := getTransfer(M) ++ getTransfer(N);
    R := getRestriction(M) ++ getRestriction(N);
    C := getConjugation(M) ++ getConjugation(N);
    p := getPrimeOrder(M);

    makeCpMackeyFunctor(p,R,T,C)
)

-- can we overload the directSum method instead of writing new directSumMackeyFunctors method?

CpMackeyFunctor ++ CpMackeyFunctor := directSumCpMackeyFunctors

ker MackeyFunctorHomomorphism := CpMackeyFunctor => options -> F -> (
    T := ker F.FixedMap;
    B := ker F.UnderlyingMap;
    
    C' := inducedMap(B,B,getConjugation(source(F)));
    T' := inducedMap(T,B,getTransfer(source(F)));
    R' := inducedMap(B,T,getRestriction(source(F)));
    p' := getPrimeOrder(source(F));

    makeCpMackeyFunctor(p',R',T',C')
)

