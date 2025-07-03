doc ///
    Key
        computeHomology
        (computeHomology,MackeyFunctorHomomorphism,MackeyFunctorHomomorphism)
    Headline
        homology of chain complex defined by two Mackey functor homomorphisms
    Usage
        computeHomology(f,g)
    Inputs
        f : MackeyFunctorHomomorphism
        N : MackeyFunctorHomomorphism
    Outputs
        : CpMackeyFunctor
            the homology of the chain complex defined by $f$ and $g$
    Description
        Text
            Given two chain complexes defined by @TO2((MackeyFunctorHomomorphism),"-Mackey functor homomorphisms")@ $f$ and $g$, we can compute the homology of the underlying complex.
///


doc ///
    Key
        (Tor,ZZ,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        computes Tor of two Cp-Mackey fuctors
    Usage
        Tor_i(M,N)
    Inputs
        i : ZZ
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the $i$th Tor Cp-Mackey Functor of $M$ and $N$
    Description
        Text
            Given two @TO2(CpMackeyFunctor,"Cp-Mackey functors")@ $M$ and $N$, we can compute their $i$th Tor.
///