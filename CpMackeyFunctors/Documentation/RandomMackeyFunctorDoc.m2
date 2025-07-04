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
        makeRandomCpMackeyFunctor(p, {n,m,k,l})
    Inputs
        p : ZZ
            a prime number $p$
        {n,m,k,l} : List
            a list of four integers, where $n$ is the number of underlying generators, $m$ is the number of fixed generators, $k$ is the number of underlying relations, and $l$ is the number of fixed relations
    Description
        Text
            Generates a @TO2((CpMackeyFunctor),"Cp-Mackey Functor")@ $M$ of prime $p$ with at most GenBound underlying generators and at most GenBound fixed generators.
///

doc ///
    Key
        makeRandomMackeyFunctorHomomorphism
        (makeRandomMackeyFunctorHomomorphism, CpMackeyFunctor, CpMackeyFunctor)
    Headline
        generate a random Cp-Mackey functor homomorphism
    Usage
        makeRandomMackeyFunctorHomomorphism(M,N)
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Description
        Text
            Generates a @TO2((MackeyFunctorHomomorphism), "Mackey functor homomorphism")@ $f$ between two @TO2((CpMackeyFunctor),"Cp-Mackey Functors")@ $M$ and $N$ (for the same prime $p$)
///


doc ///
    Key
        GenBound
    Headline
        optional input for generating random Mackey functors
    SeeAlso
        makeRandomCpMackeyFunctor
///