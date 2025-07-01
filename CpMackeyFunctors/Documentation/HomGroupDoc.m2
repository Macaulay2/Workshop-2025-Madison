doc ///
    Key
        (Hom, CpMackeyFunctor,CpMackeyFunctor)
        [Hom,DegreeLimit]
        [Hom,MinimalGenerators]
        [Hom,Strategy]
    Headline
        computes the Hom group between two Cp Mackey functors
    Usage
        Hom(M,N)
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : Module
            the abelian group of homomorphisms $\text{Hom}_{\text{Mack}_{C_p}}(M,N)$
    Description
        Text
            Given two Mackey functors $M$ and $N$, the set of Mackey functor homomorphisms from $M$ to $N$ forms an abelian group under pointwise addition.
        Example
            M = makeBurnsideMackeyFunctor(2);
            N = makeUnderlyingFreeMackeyFunctor(2);
            Hom(M,N)
///