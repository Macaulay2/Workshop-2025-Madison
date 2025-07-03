doc ///
    Key
        (boxProduct,CpMackeyFunctor,CpMackeyFunctor)
        (symbol ⊠, CpMackeyFunctor,CpMackeyFunctor)
        (symbol **,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        box product of Cp-Mackey functors
    Usage 
        M ⊠ N
        boxProduct(M,N)
        M ** N
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the box product of $M$ and $N$
    Description
        Text
            Given two @TO2((CpMackeyFunctor),"Cp-Mackey functors")@ $M$ and $N$, we can form their {\em box product}, defined by taking the levelwise box product of the underlying and fixed modules.
        Example
            boxProduct(makeUnderlyingFreeMackeyFunctor(3),makeUnderlyingFreeMackeyFunctor(3))
///

doc ///
    Key
        (boxProduct, MackeyFunctorHomomorphism,CpMackeyFunctor)
        (symbol ⊠, MackeyFunctorHomomorphism,CpMackeyFunctor)
        (symbol **, MackeyFunctorHomomorphism,CpMackeyFunctor)
        (boxProduct, CpMackeyFunctor,MackeyFunctorHomomorphism)
        (symbol ⊠, CpMackeyFunctor,MackeyFunctorHomomorphism)
        (symbol **, CpMackeyFunctor,MackeyFunctorHomomorphism)
    Headline
        induced maps on box products
    Usage 
        F ⊠ M
        boxProduct(F,M)
        F ** M
        M ⊠ F
        boxProduct(M,F)
        M ** F
    Inputs
        F : MackeyFunctorHomomorphism
        M : CpMackeyFunctor
    Outputs
        : MackeyFunctorHomomorphism
            the map induced by taking the box product of M with a morphism F.
    Description
        Text
            A map $F\colon P\to Q$ of Mackey functors can be tesnored with a Mackey functor $M$ to form a map $M ⊠ F\colon M ⊠ P\to M ⊠ Q$.
        Example
            boxProduct(makeUnderlyingFreeMackeyFunctor(3),complexLinearizationMap(3))
///
