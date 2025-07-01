-- Direct sums
CpMackeyFunctor.directSum = args -> (
    if not same (args/getPrimeOrder) then error "-- Prime not compatible";
    print "hello";
    T := directSum (args/getTransfer);
    R := directSum (args/getRestriction);
    C := directSum (args/getConjugation);
    p := getPrimeOrder(args_0);
    makeCpMackeyFunctor(p,R,T,C)
    )
CpMackeyFunctor ++ CpMackeyFunctor := CpMackeyFunctor => (F, G) -> CpMackeyFunctor.directSum(F, G)
directSum CpMackeyFunctor := CpMackeyFunctor => F -> CpMackeyFunctor.directSum(1 : F)

-- Kernels
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
