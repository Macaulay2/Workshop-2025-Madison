doc ///
    Key
        makeRandomCpMackeyFunctor
        (makeRandomCpMackeyFunctor,ZZ)
        (makeRandomCpMackeyFunctor,ZZ,List)
        [makeRandomCpMackeyFunctor, GenBound]
    Headline
        generate a random Cp-Mackey functor
    Usage
        makeRandomCpMackeyFunctor(p)
    Inputs
        p : ZZ
            a prime number $p$
    Description
        Text
            Generates a @TO2((CpMackeyFunctor),"Cp-Mackey Functor")@ $M$ of prime $p$ with at most GenBound generators
///



doc ///
    Key
        makeRandomMackeyFunctorHomomorphism
        (makeRandomMackeyFunctorHomomorphism, CpMackeyFunctor, CpMackeyFunctor)
    Headline
        generate a random Cp-Macky functor homomorphism
    Usage
        makeRandomMackeyFunctorHomomorphism(M,N)
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Description
        Text
            Generates a @TO2((MackeyFunctorHomomorphism), "Mackey functor homomorphism")@ $f$ between two compatable @TO2((CpMackeyFunctor),"Cp-Mackey Functors")@ $M$ and $N$
///

