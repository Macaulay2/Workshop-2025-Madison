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
        g : MackeyFunctorHomomorphism
    Outputs
        : CpMackeyFunctor
            the homology of the chain complex defined by $f$ and $g$
    Description
        Text
            Given two chain complexes defined by @TO2((MackeyFunctorHomomorphism),"-Mackey functor homomorphisms")@ $f$ and $g$, we can compute the homology of the underlying complex.
///

doc ///
    Key
        (Tor,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        computes Tor of two Cp-Mackey fuctors
    Usage
        Tor(M,N)
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the Tor Cp-mackey Functor of the $M$ and $N$
    Description
        Text
            Given two @TO2((CpMackeyFunctor),"-Mackey functors")@ $M$ and $N$, we can compute their Tor.
///