<<<<<<< HEAD
directSumCpMackeyFunctors = method()
directSumCpMackeyFunctors(CpMackeyFunctor,CpMackeyFunctor) := CpMackeyFunctor => (M,N) -> (
    if not getPrimeOrder(M) == getPrimeOrder(N) then error("-- Primes not the same (incompatable)");
=======
-- Direct sums
CpMackeyFunctor ++ CpMackeyFunctor := CpMackeyFunctor => (M,N) -> (
    if getPrimeOrder(M) != getPrimeOrder(N) then error("-- Primes not the same (incompatable)");
>>>>>>> 6fd72f2 (Cleaning up code and adding coker)
    T := getTransfer(M) ++ getTransfer(N);
    R := getRestriction(M) ++ getRestriction(N);
    C := getConjugation(M) ++ getConjugation(N);
    p := getPrimeOrder(M);

    makeCpMackeyFunctor(p,R,T,C)
)

<<<<<<< HEAD
-- can we overload the directSum method instead of writing new directSumMackeyFunctors method?

CpMackeyFunctor ++ CpMackeyFunctor := directSumCpMackeyFunctors

=======
-- Kernels
>>>>>>> 6fd72f2 (Cleaning up code and adding coker)
ker MackeyFunctorHomomorphism := CpMackeyFunctor => options -> F -> (
    T := ker F.FixedMap;
    B := ker F.UnderlyingMap;
    
    C' := inducedMap(B,B,getConjugation(source(F)));
    T' := inducedMap(T,B,getTransfer(source(F)));
    R' := inducedMap(B,T,getRestriction(source(F)));
    p' := getPrimeOrder(source(F));

    makeCpMackeyFunctor(p',R',T',C')
)

-- Cokernels
coker MackeyFunctorHomomorphism := CpMackeyFunctor => F -> (
    T := coker F.FixedMap;
    B := coker F.UnderlyingMap;
    
    C' := inducedMap(B,B,getConjugation(target(F)));
    T' := inducedMap(T,B,getTransfer(target(F)));
    R' := inducedMap(B,T,getRestriction(target(F)));
    p' := getPrimeOrder(target(F));

    makeCpMackeyFunctor(p',R',T',C')
)
